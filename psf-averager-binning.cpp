#include "psf-averager-photoz.hpp"
#include "rebin.hpp"

struct binner_options {
    vec1f bins = {300.0, 1500.0};
    vec1f errors = {0.01, 0.03, 0.1, 0.3, 1.0};
    vec1f biases = {0.00001, 0.00003, 0.0001, 0.0003, 0.001};
    vec1s rebin = {"cst", "trapz", "spline"};
};

class bin_averager : public psf_averager {
public :
    // Individuals
    vec<6,metrics> mb;
    vec<1,metrics> mtrue;

    // Options
    binner_options opts;

    // Internal variables
    vec1d uvj_ngal;

    bin_averager() : psf_averager("bin") {}

    void configure_binner(binner_options topts) {
        opts = topts;
    }

    void process_cached(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {}

    void do_fit(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        uint_t ised = (id_bt == 0 ? id_disk : id_bulge);
        uvj_ngal.safe[ised] += tngal;
    }

    void initialize_redshift_bin(uint_t iz) override {}

    std::string make_cache_hash() override {
        return hash(opts.bins, opts.errors, opts.rebin);
    }

    void initialize_redshift_slice(uint_t itz) override {
        uvj_ngal.resize(use.size());
    }

    void initialize_cache() override {}

    void finalize_redshift_slice(uint_t itz) override {
        double zf = uzf[itz];

        mb.resize(use.size(), opts.bins.size(), opts.rebin.size(), opts.errors.size(), opts.biases.size(), nmc);
        mtrue.resize(use.size());

        uint_t ncomb = opts.bins.size()*opts.errors.size()*opts.biases.size();
        auto pg = progress_start(count(use)*ncomb);
        for (uint_t iuv : range(use.dims[0]))
        for (uint_t ivj : range(use.dims[1])) {
            if (!use.safe(iuv,ivj)) continue;

            vec1d tlam = lam(iuv,ivj,_);
            vec1d tsed = sed(iuv,ivj,_);

            if (!naive_igm) {
                apply_madau_igm(zf, tlam, tsed);
            }

            tsed = lsun2uJy(zf, 1.0, tlam, tsed);
            tlam *= (1.0 + zf);

            uint_t ised = iuv*use.dims[1] + ivj;

            mtrue.safe[ised] = metrics(egg_q11.safe[ised], egg_q12.safe[ised], egg_q22.safe[ised]);

            for (uint_t ib : range(opts.bins)) {
                // Get binned photometry
                double lmin = 0.1, lmax = 2.0;
                double dl = 1e-4*opts.bins[ib];
                uint_t npt = (lmax - lmin)/dl;
                vec1d blam = dl*indgen<double>(npt) + lmin;
                vec1d fobs(npt);
                for (uint_t l : range(npt)) {
                    if (blam[l]-dl/2.0 > tlam.front() && blam[l]+dl/2.0 < tlam.back()) {
                        fobs.safe[l] = integrate(tlam, tsed, blam[l]-dl/2.0, blam[l]+dl/2.0)/dl;
                    }
                }

                uint_t nlam = blam.size();
                vec2d rnd_noise = randomn(seed, nmc, nlam);

                for (uint_t ie : range(opts.errors))
                for (uint_t ibe : range(opts.biases)) {
                    // Normalize errors to zero bias
                    vec2d cerror(rnd_noise.dims);
                    for (uint_t l : range(blam)) {
                        vec1d tmp = e10(0.4*opts.errors.safe[ie]*rnd_noise.safe(_,l));
                        cerror.safe(_,l) = tmp*(e10(0.4*opts.biases.safe[ibe]*(l-nlam/2.0)*dl*1e3)/mean(tmp));
                    }

                    if (ib == 0 && ised == where_first(use)) {
                        vec2d correl(nlam, nlam);
                        for (uint_t i : range(nlam))
                        for (uint_t j : range(nlam)) {
                            if (j >= i) {
                                correl(i,j) = mean((cerror(_,i) - 1.0)*(cerror(_,j) - 1.0));
                            } else {
                                correl(i,j) = correl(j,i);
                            }
                        }

                        fits::write("correl_"+to_string(ie)+"-"+to_string(ibe)+".fits", correl);
                    }

                    for (uint_t i : range(nmc)) {
                        // Generate perturbed photometry
                        vec1f tobs = fobs;
                        for (uint_t l : range(blam)) {
                            tobs.safe[l] *= cerror.safe(i,l);
                        }

                        for (uint_t im : range(opts.rebin)) {
                            // Interpolate
                            vec1d robs;
                            if (opts.rebin[im] == "cst") {
                                robs = rebin_cst(tobs, blam, psf_filter.lam);
                            } else if (opts.rebin[im] == "trapz") {
                                robs = rebin_trapz(tobs, blam, psf_filter.lam);
                            } else if (opts.rebin[im] == "lin") {
                                // No need: actually indistinguishable from trapz
                                robs = interpolate(tobs, blam, psf_filter.lam);
                            } else if (opts.rebin[im] == "spline") {
                                robs = rebin_mcspline(tobs, blam, psf_filter.lam);
                            }

                            // Compute PSF
                            robs *= psf_filter.res;
                            double fvis = integrate(psf_filter.lam, robs);
                            mb.safe(ised,ib,im,ie,ibe,i) = metrics(
                                integrate(psf_filter.lam, robs*mono_q11)/fvis,
                                integrate(psf_filter.lam, robs*mono_q12)/fvis,
                                integrate(psf_filter.lam, robs*mono_q22)/fvis
                            );

                            // if (ib == 0 && ie == 5 && i == 0) {
                            //     fits::write_table("rebin_"+opts.rebin[im]+".fits",
                            //         "blam", blam, "bflx", tobs,
                            //         "elam", psf_filter.lam, "eflx", robs/psf_filter.res,
                            //         "olam", tlam, "oflx", tsed
                            //     );
                            // }
                        }
                    }

                    progress(pg);
                }
            }
        }

        vec1f uv = -2.5*log10(egg_fu/egg_fv);
        vec1f vj = -2.5*log10(egg_fv/egg_fj);

        // Write to disk the individual measurements
        file::mkdir(cache_dir);
        fits::write_table(indiv_filename,
            "e1_true",  get_e1(mtrue),
            "e2_true",  get_e2(mtrue),
            "r2_true",  get_r2(mtrue),
            "e1_obs",   get_e1(mb),
            "e2_obs",   get_e2(mb),
            "r2_obs",   get_r2(mb),
            "bins",     opts.bins,
            "errors",   opts.errors,
            "biases",   opts.biases,
            "rebin",    opts.rebin,
            "ngal",     uvj_ngal,
            "use",      use,
            "uv",       uv,
            "vj",       vj
        );
    }

    void finalize_redshift_bin(uint_t iz, double ntot, double ntot_qu, double ntot_sf) override {}
};

int vif_main(int argc, char* argv[]) {
    // External data
    std::string share_dir = "/home/cschreib/code/egg-analytic/share/";
    std::string filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
    std::string sed_lib = "/home/cschreib/code/egg-analytic/share/opt_lib_fastpp_hd_noigm.fits";
    std::string sed_imf = "chabrier";
    std::string psf_file  = "/home/cschreib/code/euclid_psf/psf-averager/psf-mono.fits";

    // Survey definition
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";

    // Mock photometry
    uint_t nmc = 100;
    double dz = 0.01;
    uint_t seds_step = 1;
    bool write_cache = false;
    bool use_cache = false;
    bool write_individuals = false;
    bool write_averages = false;
    uint_t nthread = 0;
    uint_t iz = 5;
    std::string cache_id = "rebin";
    std::string cache_dir = "cache";
    double prob_limit = 0.1;
    vec1f bins = {300.0, 1500.0};
    vec1f errors = {0.01, 0.03, 0.1, 0.3, 1.0};
    vec1f biases = {0.00001, 0.00003, 0.0001, 0.0003, 0.001};
    vec1s rebin = {"cst", "lin", "trapz", "spline"};

    read_args(argc, argv, arg_list(
        maglim, selection_band, nmc, dz, seds_step, write_cache, use_cache, iz, share_dir,
        filter_db, psf_file, nthread, write_individuals, write_averages, cache_id, sed_lib,
        sed_imf, cache_dir, prob_limit, bins, errors, biases, rebin
    ));

    bin_averager bavg;
    bavg.write_cache = write_cache;
    bavg.use_cache = use_cache;

    // Setup survey
    egg::generator_options opts;
    opts.share_dir = share_dir;
    opts.filter_db = filter_db;
    opts.sed_lib = sed_lib;
    opts.sed_lib_imf = sed_imf;
    opts.filter_flambda = true;
    opts.filter_photons = true;
    opts.trim_filters = true;
    opts.selection_band = selection_band;
    opts.maglim = maglim;
    opts.logmass_steps = 50;
    opts.bt_steps = 2;
    opts.logmass_max = 12.0;
    opts.seds_step = seds_step;
    opts.nthread = nthread;
    bavg.initialize(opts);

    // Setup mock
    mock_options mopts;
    mopts.nmc = nmc;
    mopts.dz = dz;
    mopts.psf_file = psf_file;
    mopts.write_individuals = write_individuals;
    mopts.write_averages = write_averages;
    mopts.force_cache_id = cache_id;
    mopts.cache_dir = cache_dir;
    mopts.prob_limit = prob_limit;
    bavg.configure_mock(mopts);

    // Setup binning
    binner_options bopts;
    bopts.bins = bins;
    bopts.errors = errors;
    bopts.biases = biases;
    bopts.rebin = rebin;
    bavg.configure_binner(bopts);

    // Average PSF metrics
    // for (uint_t iz : range(bavg.zb)) {
    //     if (!bavg.average_redshift_bin(iz)) continue;
    // }

    if (!bavg.average_redshift_bin(iz)) return 1;

    return 0;
}
