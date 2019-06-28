#include "psf-averager-photoz.hpp"
#include "rebin.hpp"

struct binner_options {
    vec1f bins = {5.0, 10.0, 20.0, 50.0, 100.0, 300.0, 500.0, 1500.0};
    vec1s rebin = {"cst", "trapz", "spline", "mcspline"};
};

struct shared_data {
    vec<3,metrics> psf; // [niter,nmeth,nbin]
    vec<3,metrics> lib; // [nuvj,nmeth,nbin]
    binner_options opts;
};

class bin_averager : public psf_averager {
public :
    // Individuals
    shared_data& shared;

    bin_averager(filter_database& db, psf_moments& pm, fitter_base& f, shared_data& sh) :
        psf_averager(db, pm, f), shared(sh) {}

    void on_generated(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        psf_averager::on_generated(iter, id_mass, id_type, id_disk, id_bulge, id_bt, tngal,
            fdisk, fbulge);

        if (just_count) {
            // Just counting number of iterations
            return;
        }

        const double fpsf_bulge = egg_fvis.safe[id_bulge]*bt.safe[id_bt];
        const double fpsf_disk = egg_fvis.safe[id_disk]*(1.0 - bt.safe[id_bt]);
        const double fbtn = fpsf_bulge/(fpsf_disk + fpsf_bulge);
        const double fbti = 1.0 - fbtn;

        for (uint_t ib : range(shared.opts.bins))
        for (uint_t im : range(shared.opts.rebin)) {
            auto& pd = shared.lib(id_disk,im,ib);
            auto& pb = shared.lib(id_bulge,im,ib);
            shared.psf(iter,im,ib) = metrics(
                pd.q11*fbti  + pb.q11*fbtn,
                pd.q12*fbti  + pb.q12*fbtn,
                pd.q22*fbti  + pb.q22*fbtn,
                pd.rlam*fbti + pb.rlam*fbtn
            );
        }
    }
};

class null_fitter : public fitter_base {
public:
    shared_data& shared;
    egg::generator* generator = nullptr;

    explicit null_fitter(filter_database& db, psf_moments& pm, shared_data& sh) :
        fitter_base(db, pm), shared(sh) {}

    void do_prepare_fit(double zf) override {
        note("preparing library at z=", zf, "...");

        uint_t nvj = generator->use.dims[1];
        uint_t nsed = generator->use.size();

        shared.psf.resize(niter, shared.opts.rebin.size(), shared.opts.bins.size());
        shared.lib.resize(nsed, shared.opts.rebin.size(), shared.opts.bins.size());

        for (uint_t ised : range(nsed)) {
            if (!generator->use[ised]) continue;

            uint_t iuv = ised / nvj;
            uint_t ivj = ised % nvj;

            vec1d tlam = generator->lam(iuv,ivj,_);
            vec1d tsed = generator->sed(iuv,ivj,_);

            if (!generator->naive_igm) {
                generator->apply_madau_igm(zf, tlam, tsed);
            }

            double df = lumdist(zf, generator->cosmo);
            tsed = lsun2uJy(zf, df, tlam, tsed);
            tlam *= (1.0 + zf);

            for (uint_t ib : range(shared.opts.bins)) {
                // Get binned photometry
                double lmin = 0.1, lmax = 2.0;
                double dl = 1e-4*shared.opts.bins[ib];
                uint_t npt = (lmax - lmin)/dl;
                vec1d blam = dl*indgen<double>(npt) + lmin;
                vec1d fobs(npt);
                for (uint_t l : range(npt)) {
                    if (blam[l]-dl/2.0 > tlam.front() && blam[l]+dl/2.0 < tlam.back()) {
                        fobs.safe[l] = integrate(tlam, tsed, blam[l]-dl/2.0, blam[l]+dl/2.0)/dl;
                    }
                }

                for (uint_t im : range(shared.opts.rebin)) {
                    // Interpolate
                    vec1d robs;
                    if (shared.opts.rebin[im] == "cst") {
                        robs = rebin_cst(fobs, blam, psf.mono_lam);
                    } else if (shared.opts.rebin[im] == "trapz") {
                        robs = rebin_trapz(fobs, blam, psf.mono_lam);
                    } else if (shared.opts.rebin[im] == "spline") {
                        robs = rebin_spline3(fobs, blam, psf.mono_lam);
                    } else if (shared.opts.rebin[im] == "mcspline") {
                        robs = rebin_mcspline(fobs, blam, psf.mono_lam);
                    }

                    double fvis, q11, q12, q22, rlam;
                    psf.get_moments_same_grid(robs, q11, q12, q22, fvis, rlam);

                    shared.lib(ised,im,ib) = metrics(q11, q12, q22, rlam);
                }
            }
        }
    }

    fit_result do_fit(uint_t iter, const vec1d& ftot) override {
        fit_result fr(nmc);
        return fr;
    }

    void save_individuals(const std::string& filename) override {
        // Write to disk the individual measurements
        fits::table otbl(filename);
        otbl.update_column("e1_bin", get_e1(shared.psf));
        otbl.update_column("e2_bin", get_e2(shared.psf));
        otbl.update_column("r2_bin", get_r2(shared.psf));
        otbl.update_column("rlam_bin", get_rlam(shared.psf));
        otbl.update_column("bins", shared.opts.bins);
        otbl.update_column("rebin", shared.opts.rebin);
    }
};

int vif_main(int argc, char* argv[]) {
    uint_t iz = 5;
    bool allbins = false;

    shared_data shared;
    filter_database db;
    psf_moments psf(db);
    null_fitter fitter(db, psf, shared);
    bin_averager pavg(db, psf, fitter, shared);

    fitter.generator = &pavg;

    // Read setup
    {
        program_arguments opts(argc, argv);

        // Override options
        {
            bool filter_flambda = true; // equivalent to FILTER_FORMAT=1
            bool filter_photons = true; // equivalent to FILTER_FORMAT=1
            bool trim_filters = true;
            opts.write(arg_list(filter_flambda, filter_photons, trim_filters));
        }

        // Default options
        {
            std::string filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
            opts.read(arg_list(filter_db));
        }

        // Setup binning
        opts.read(arg_list(shared.opts.bins, shared.opts.rebin));

        // Setup filter database
        db.read_options(opts);

        // Setup PSF moments
        psf.read_options(opts);

        // Setup survey
        pavg.read_options(opts);

        // Setup fitter
        fitter.read_options(opts);

        opts.read(arg_list(iz, allbins));
    }

    // Average PSF metrics
    if (allbins) {
        for (uint_t tiz : range(pavg.zb.size()-1)) {
            if (!pavg.average_redshift_bin(tiz)) return 1;
        }
    } else {
        if (!pavg.average_redshift_bin(iz)) return 1;
    }

    return 0;
}
