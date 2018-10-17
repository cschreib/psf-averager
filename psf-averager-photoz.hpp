#ifndef PSF_AVERAGER_PHOTOZ_INCLUDED
#define PSF_AVERAGER_PHOTOZ_INCLUDED

#include "egg-analytic.hpp"
#include "fitter-base.hpp"
#include "common.hpp"
#include "metrics.hpp"

struct averager_options {
    double dz = 0.01;
    uint_t nzsplit = 0;
    uint_t nmc = 200;
    vec1f depths;
    bool no_noise = false;
    bool keep_individuals_in_memory = false;
    bool write_individuals = false;
    bool keep_averages_in_memory = false;
    bool write_averages = false;
    bool save_seds = false;
    std::string force_cache_id;
    std::string cache_dir;
    double prob_limit = 0.1;
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    bool global_progress_bar = false;
};

class psf_averager : public egg::generator {
public :
    // Total number of galaxies
    double ngal = 0;

    // Averages
    std::mutex avg_mutex;
    metrics m_tr; // EGG (truth)
    metrics m_ml; // maximum likelihood
    metrics m_ma; // marginalization
    vec1d uzf;
    vec1d dndz;

    // Individuals
    vec1u indiv_id_mass, indiv_id_type, indiv_id_disk, indiv_id_bulge, indiv_id_bt;
    vec1d indiv_ngal;
    vec1f indiv_uv, indiv_vj;
    vec2f indiv_fdisk, indiv_fbulge;
    vec<1,metrics> indiv_tr; // EGG (truth)
    vec<1,metrics> zml, zma;
    vec<2,metrics> indiv_ml, indiv_ma;
    vec2f indiv_chi2, indiv_zml, indiv_zma;

    // Monochromatic PSF library
    psf_moments& psf;

    // Observed SEDs
    vec1f save_sed_lambda;     // [nlam]
    vec2f save_sed_egg_fluxes; // [nsed,nlam]

    // EGG PSF library
    vec1d egg_q11, egg_q12, egg_q22;
    vec1d egg_fu, egg_fv, egg_fj, egg_fvis;

    bool just_count = false;
    uint_t niter = 0;
    double ngal_tot = 0.0;

    std::mutex pg_mutex;
    progress_t pgi;

    // Config
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    double dz = 0.01;
    uint_t nzsplit = 0;
    uint_t nband = npos, nmc = npos;
    double prob_limit = 0.1;

    bool single_pass = false;
    bool global_progress_bar = false;

    fitter_base& fitter;

    bool no_noise = false;

    // Cache
    std::mutex cache_mutex;
    std::string cache_filename;
    fits::table fitter_cache;
    bool write_cache = false;
    std::string force_cache_id;
    std::string cache_dir;

    // Individual measurements
    bool keep_individuals_in_memory = false;
    bool write_individuals = false;
    bool keep_averages_in_memory = false;
    bool write_averages = false;
    bool save_seds = false;

    explicit psf_averager(filter_database& db, psf_moments& pm, fitter_base& f) :
        egg::generator(db), psf(pm), fitter(f) {}

    void read_options(program_arguments& opts) {
        egg::generator_options gopts;
        averager_options aopts;

        // Setup generator
        gopts.bt_steps = 5;
        gopts.logmass_max = 12.0;
        gopts.logmass_steps = 50.0;
        gopts.nthread = 0;

        // Survey definition
        gopts.maglim = 24.5;
        gopts.selection_band = "euclid-vis";
        gopts.filters = {"sdss-u", "sdss-g", "sdss-r", "sdss-i", "sdss-z", "vista-Y", "vista-J", "vista-H"};
        aopts.depths =  {24.5,     24.4,     24.1,     24.1,     23.7,     23.2,      23.2,      23.2}; // 10 sigma
        gopts.seds_step = 5;

        // Averager options
        aopts.nmc = 200;
        aopts.dz = 0.01;
        aopts.nzsplit = 0;
        aopts.no_noise = false;
        write_cache = false;
        aopts.write_individuals = true;
        aopts.write_averages = false;
        aopts.cache_dir = "cache";
        aopts.prob_limit = 0.1;
        aopts.save_seds = false;
        aopts.global_progress_bar = false;

        // External data
        gopts.share_dir = "/home/cschreib/code/egg-analytic/share/";
        gopts.sed_lib = "/home/cschreib/code/egg-analytic/share/opt_lib_fastpp_hd_noigm.fits";
        gopts.sed_lib_imf = "chabrier";

        // Read custom configuration from command line
        opts.read(arg_list(
            gopts.share_dir, gopts.sed_lib, name(gopts.sed_lib_imf, "sed_imf"),
            gopts.selection_band, gopts.filters, gopts.maglim, name(gopts.logmass_steps, "mass_steps"),
            gopts.seds_step, gopts.nthread, write_cache, aopts.depths, aopts.nmc,
            aopts.dz, aopts.nzsplit, aopts.no_noise, aopts.write_individuals,
            aopts.write_averages, name(aopts.force_cache_id, "cache_id"), aopts.cache_dir,
            aopts.prob_limit, aopts.save_seds, aopts.global_progress_bar, aopts.zb
        ));

        // Adjust
        if (aopts.write_individuals) {
            aopts.keep_individuals_in_memory = true;
            opts.write(arg_list(aopts.keep_individuals_in_memory));
        }
        if (aopts.write_averages) {
            aopts.keep_averages_in_memory = true;
            opts.write(arg_list(aopts.keep_averages_in_memory));
        }
        if (aopts.no_noise && aopts.nmc != 1) {
            note("no_noise set, setting nmc=1");
            aopts.nmc = 1;
            opts.write(arg_list(aopts.nmc));
        }

        // Initialize
        egg::generator::initialize(gopts);
        initialize(aopts);
    }

    void initialize(const averager_options& aopts) {
        nband = filters.size();
        dz = aopts.dz;
        nzsplit = aopts.nzsplit;
        zb = aopts.zb;
        no_noise = aopts.no_noise;
        keep_individuals_in_memory = aopts.keep_individuals_in_memory;
        write_individuals = aopts.write_individuals;
        keep_averages_in_memory = aopts.keep_averages_in_memory;
        write_averages = aopts.write_averages;
        force_cache_id = aopts.force_cache_id;
        cache_dir = file::directorize(aopts.cache_dir);
        prob_limit = aopts.prob_limit;
        if (prob_limit == 0.0) prob_limit = dnan;
        save_seds = aopts.save_seds;
        global_progress_bar = aopts.global_progress_bar;
        nmc = aopts.nmc;

        if (save_seds) {
            file::mkdir("seds/");
            save_sed_lambda = rgen_step(0.1,2.0,0.001);
        }

        // Compute actual rest-frame colors of EGG SEDs.
        filter_t rest_filter_u, rest_filter_v, rest_filter_j;
        rest_filter_u = filter_db.read_filter("maiz-U");
        rest_filter_v = filter_db.read_filter("maiz-V");
        rest_filter_j = filter_db.read_filter("2mass-J");

        auto compute_uvj = [&](uint_t ised, const vec1d& tlam, vec1d tsed) {
            tsed = lsun2uJy(0.0, 1.0, tlam, tsed);

            egg_fu[ised] = sed2flux(rest_filter_u.lam, rest_filter_u.res, tlam, tsed);
            egg_fv[ised] = sed2flux(rest_filter_v.lam, rest_filter_v.res, tlam, tsed);
            egg_fj[ised] = sed2flux(rest_filter_j.lam, rest_filter_j.res, tlam, tsed);
        };

        if (single_sed_library) {
            egg_fu.resize(single_use.size());
            egg_fv.resize(single_use.size());
            egg_fj.resize(single_use.size());
            for (uint_t ised : range(single_use)) {
                if (!single_use[ised]) continue;
                compute_uvj(ised, single_lam(ised,_), single_sed(ised,_));
            }
        } else {
            egg_fu.resize(use.size());
            egg_fv.resize(use.size());
            egg_fj.resize(use.size());
            for (uint_t iuv : range(use.dims[0]))
            for (uint_t ivj : range(use.dims[1])) {
                if (!use(iuv, ivj)) continue;
                compute_uvj(iuv*use.dims[1]+ivj, lam(iuv,ivj,_), sed(iuv,ivj,_));
            }
        }
    }

    void on_generated(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        if (just_count) {
            // This is always single-threaded, no need to worry
            ++niter;
            ngal_tot += tngal;
            return;
        }

        // Compute true PSF moments
        // ------------------------

        // Compute flux-weighted B/T
        const double fpsf_bulge = egg_fvis.safe[id_bulge]*bt.safe[id_bt];
        const double fpsf_disk = egg_fvis.safe[id_disk]*(1.0 - bt.safe[id_bt]);
        const double fbtn = fpsf_bulge/(fpsf_disk + fpsf_bulge);
        const double fbti = 1.0 - fbtn;

        // Add up PSFs
        metrics tr(
            egg_q11.safe[id_disk]*fbti + egg_q11.safe[id_bulge]*fbtn,
            egg_q12.safe[id_disk]*fbti + egg_q12.safe[id_bulge]*fbtn,
            egg_q22.safe[id_disk]*fbti + egg_q22.safe[id_bulge]*fbtn
        );

        // Compute actual UVJ colors
        double mbt = bt.safe[id_bt];
        double mbti = 1.0 - mbt;
        double rfu = egg_fu.safe[id_disk]*mbti + mbt*egg_fu.safe[id_bulge];
        double rfv = egg_fv.safe[id_disk]*mbti + mbt*egg_fv.safe[id_bulge];
        double rfj = egg_fj.safe[id_disk]*mbti + mbt*egg_fj.safe[id_bulge];
        double rfuv = -2.5*log10(rfu/rfv);
        double rfvj = -2.5*log10(rfv/rfj);

        // Now generate mocks and fit them
        // -------------------------------

        // Do the fitting
        vec1d ftot(fdisk.size()-1);
        for (uint_t l : range(ftot)) {
            ftot.safe[l] = fdisk.safe[l+1] + fbulge.safe[l+1];
        }

        fit_result fr = fitter.do_fit(iter, ftot);

        bool has_fit = !fr.z_obs.empty();

        // Compute averages
        metrics ma, ma2, ml, ml2;
        if (has_fit) {
            for (uint_t i : range(nmc)) {
                ma  += fr.psf_obs.safe[i];
                ma2 += sqr(fr.psf_obs.safe[i]);
                ml  += fr.psf_obsm.safe[i];
                ml2 += sqr(fr.psf_obsm.safe[i]);
            }

            ma  /= nmc;
            ma2 /= nmc;
            ma2 = sqrt(ma2 - sqr(ma));
            ml  /= nmc;
            ml2 /= nmc;
            ml2 = sqrt(ml2 - sqr(ml));
        }

        // Store individual values if asked
        if (keep_individuals_in_memory) {
            // Could be executed concurrently, but never at the same 'iter'.
            // Therefore this is thread safe.
            indiv_tr.safe[iter]       = tr;
            indiv_id_mass.safe[iter]  = id_mass;
            indiv_id_type.safe[iter]  = id_type;
            indiv_id_disk.safe[iter]  = id_disk;
            indiv_id_bulge.safe[iter] = id_bulge;
            indiv_id_bt.safe[iter]    = id_bt;
            indiv_ngal.safe[iter]     = tngal;
            indiv_uv.safe[iter]       = rfuv;
            indiv_vj.safe[iter]       = rfvj;
            indiv_fdisk.safe(iter,_)  = fdisk;
            indiv_fbulge.safe(iter,_) = fbulge;

            if (has_fit) {
                for (uint_t i : range(nmc)) {
                    indiv_chi2.safe(iter,i) = fr.chi2.safe[i];
                    indiv_zml.safe(iter,i) = fr.z_obs.safe[i];
                    indiv_zma.safe(iter,i) = fr.z_obsm.safe[i];
                    indiv_ml.safe(iter,i) = fr.psf_obs.safe[i];
                    indiv_ma.safe(iter,i) = fr.psf_obsm.safe[i];
                }
            }
        }

        // Save SEDs if asked
        if (save_seds) {
            vec1f tsed = nm.safe[id_mass]*(
                mbti*save_sed_egg_fluxes(id_disk,_) +
                mbt*save_sed_egg_fluxes(id_bulge,_)
            );

            fits::update_table("seds/sed_"+to_string(iter)+".fits",
                "lam", save_sed_lambda, "sed_true", tsed
            );
        }

        // Store average values if asked
        if (keep_averages_in_memory) {
            // Could be executed concurrently, use mutex when in multithreading context
            auto lock = (nthread > 0 ?
                std::unique_lock<std::mutex>(avg_mutex) : std::unique_lock<std::mutex>());

            m_tr += tngal*tr;
            if (has_fit) {
                m_ml += tngal*ml;
                m_ma += tngal*ma;
            }
            ngal += tngal;
        }

        // Save things in the cache
        if (write_cache) {
            // Could be executed concurrently, use mutex when in multithreading context
            auto lock = (nthread > 0 ?
                std::unique_lock<std::mutex>(cache_mutex) : std::unique_lock<std::mutex>());

            fitter_cache.update_elements("im",        id_mass,      fits::at(iter));
            fitter_cache.update_elements("it",        id_type,      fits::at(iter));
            fitter_cache.update_elements("idisk",     id_disk,      fits::at(iter));
            fitter_cache.update_elements("ibulge",    id_bulge,     fits::at(iter));
            fitter_cache.update_elements("ibt",       id_bt,        fits::at(iter));
            fitter_cache.update_elements("ngal",      tngal,        fits::at(iter));
            fitter_cache.update_elements("fbulge",    fbulge,       fits::at(iter,_));
            fitter_cache.update_elements("fdisk",     fdisk,        fits::at(iter,_));
            fitter_cache.update_elements("uv",        rfuv,         fits::at(iter));
            fitter_cache.update_elements("vj",        rfvj,         fits::at(iter));
            fitter_cache.update_elements("e1_true",   tr.e1,        fits::at(iter));
            fitter_cache.update_elements("e2_true",   tr.e2,        fits::at(iter));
            fitter_cache.update_elements("r2_true",   tr.r2,        fits::at(iter));

            if (has_fit) {
                fitter_cache.update_elements("chi2_obs",    fr.chi2,    fits::at(iter,_));
                fitter_cache.update_elements("z_obs",       fr.z_obs,   fits::at(iter,_));
                fitter_cache.update_elements("z_obsm",      fr.z_obsm,  fits::at(iter,_));
                fitter_cache.update_elements("e1_obs",      ml.e1,      fits::at(iter));
                fitter_cache.update_elements("e2_obs",      ml.e2,      fits::at(iter));
                fitter_cache.update_elements("r2_obs",      ml.r2,      fits::at(iter));
                fitter_cache.update_elements("e1_obsm",     ma.e1,      fits::at(iter));
                fitter_cache.update_elements("e2_obsm",     ma.e2,      fits::at(iter));
                fitter_cache.update_elements("r2_obsm",     ma.r2,      fits::at(iter));
                fitter_cache.update_elements("e1_obs_err",  ml2.e1,     fits::at(iter));
                fitter_cache.update_elements("e2_obs_err",  ml2.e2,     fits::at(iter));
                fitter_cache.update_elements("r2_obs_err",  ml2.r2,     fits::at(iter));
                fitter_cache.update_elements("e1_obsm_err", ma2.e1,     fits::at(iter));
                fitter_cache.update_elements("e2_obsm_err", ma2.e2,     fits::at(iter));
                fitter_cache.update_elements("r2_obsm_err", ma2.r2,     fits::at(iter));
            }

            fitter_cache.open(cache_filename);
        }

        if (!global_progress_bar) {
            // Could be executed concurrently, use mutex when in multithreading context
            auto lock = (nthread > 0 ?
                std::unique_lock<std::mutex>(pg_mutex) : std::unique_lock<std::mutex>());

            progress(pgi);
        }
    }

    bool average_redshift_bin(uint_t iz) {
        uint_t ntz;
        if (nzsplit != 0) {
            ntz = nzsplit;
            dz = (zb[iz+1] - zb[iz])/nzsplit;
            uzf = zb[iz] + dz*(0.5 + indgen<double>(nzsplit));
        } else {
            ntz = max(floor((zb[iz+1]-zb[iz])/dz), 1);
            uzf = zb[iz] + dz*indgen<double>(ntz);
        }

        // Start averaging
        vec<1,metrics> ztr;
        if (keep_averages_in_memory) {
            ztr.resize(ntz);
            zml.resize(ntz);
            zma.resize(ntz);
            dndz.resize(ntz);
        }

        if (global_progress_bar) {
            pgi = progress_start(ntz);
        }

        for (uint_t itz : range(ntz)) {
            if (!global_progress_bar) {
                print(itz, "/", ntz);
            }

            double zf = uzf.safe[itz];
            double df = lumdist(zf, cosmo);
            const double ML_cor = e10(-interpolate({0.15,0.15,0.0,0.0,-0.6}, {0.0,0.45,1.3,6.0,8.0}, zf));

            uint_t precision = ceil(max(-log10(dz), 2.0));
            std::string zid = replace(to_string(format::fixed(format::precision(zf, precision))), ".", "p");
            std::string cache_id;
            if (force_cache_id.empty()) {
                cache_id = hash(fitter.make_cache_hash(), niter);
            } else {
                cache_id = force_cache_id;
            }

            if (write_individuals) {
                std::string indiv_filename = cache_dir+"indiv-"+fitter.code_name+"-z"+zid+"-"+cache_id+".fits";
                if (file::exists(indiv_filename)) continue;
            }

            // Compute PSF moments for each template in the library
            auto set_moments = [&](uint_t ised, vec1d tlam, vec1d tsed) {
                if (!naive_igm) {
                    apply_madau_igm(zf, tlam, tsed);
                }

                tsed = lsun2uJy(zf, df, tlam, tsed);
                tlam *= (1.0 + zf);

                if (save_seds) {
                    save_sed_egg_fluxes(ised,_) = ML_cor*resample_sed(tsed, tlam, save_sed_lambda);
                }

                psf.get_moments(
                    tlam, tsed, egg_q11[ised], egg_q12[ised], egg_q22[ised], egg_fvis[ised]
                );
            };

            if (single_sed_library) {
                egg_q11.resize(single_use.size());
                egg_q12.resize(single_use.size());
                egg_q22.resize(single_use.size());
                egg_fvis.resize(single_use.size());

                if (save_seds) {
                    save_sed_egg_fluxes.resize(single_use.size(), save_sed_lambda.size());
                }

                for (uint_t ised : range(single_use)) {
                    if (!single_use[ised]) continue;
                    set_moments(ised, single_lam(ised,_), single_sed(ised,_));
                }
            } else {
                egg_q11.resize(use.size());
                egg_q12.resize(use.size());
                egg_q22.resize(use.size());
                egg_fvis.resize(use.size());

                if (save_seds) {
                    save_sed_egg_fluxes.resize(use.size(), save_sed_lambda.size());
                }

                for (uint_t iuv : range(use.dims[0]))
                for (uint_t ivj : range(use.dims[1])) {
                    if (!use(iuv,ivj)) continue;
                    set_moments(iuv*use.dims[1]+ivj, lam(iuv,ivj,_), sed(iuv,ivj,_));
                }
            }

            // Reset averages
            if (keep_averages_in_memory) {
                m_tr.reset();
                m_ma.reset();
                m_ml.reset();
                ngal = 0;
            }

            // Pre-compute number of iterations
            if (!single_pass || keep_individuals_in_memory || is_finite(prob_limit)) {
                note("compute number of iterations...");

                // Initialize to zero
                niter = 0;
                ngal_tot = 0.0;
                just_count = true;
                use_prob_cut = false;

                // Disable multi-threading for this
                uint_t onthread = nthread;
                nthread = 0;

                // Count number of iterations and total number density
                generate(zf, dz);

                if (is_finite(prob_limit)) {
                    // Count again, now excluding low probability objects
                    use_prob_cut = true;
                    prob_cut = prob_limit*(ngal_tot/niter);
                    niter = 0;
                    ngal_tot = 0.0;
                    generate(zf, dz);
                }

                // Reset multithreading to its original state
                nthread = onthread;
                just_count = false;

                note("done: ", niter);
            }

            // Prepare fitter
            fitter.prepare_fit(niter, zf);

            // Initialize cache
            if (write_cache) {
                cache_filename = cache_dir+"cache-"+fitter.code_name+"-z"+zid+"-"+cache_id+".fits";
                note("cache file: ", cache_filename);
            }

            vec1s tbands(nband+1);
            vec1f lambda(nband+1);
            tbands[0] = selection_band;
            lambda[0] = selection_filter.rlam;
            if (nband > 0) {
                tbands[1-_] = bands;
                for (uint_t l : range(filters)) {
                    lambda[l+1] = filters[l].rlam;
                }
            }

            if (write_cache) {
                note("creating cache on disk...");

                file::mkdir(cache_dir);
                file::remove(cache_filename);
                fitter_cache.open(cache_filename);

                // Create cache arrays
                fitter_cache.allocate_column<uint_t>("im",     niter);
                fitter_cache.allocate_column<uint_t>("it",     niter);
                fitter_cache.allocate_column<uint_t>("idisk",  niter);
                fitter_cache.allocate_column<uint_t>("ibulge", niter);
                fitter_cache.allocate_column<uint_t>("ibt",    niter);
                fitter_cache.allocate_column<float>("fbulge",  niter, nband+1);
                fitter_cache.allocate_column<float>("fdisk",   niter, nband+1);
                fitter_cache.allocate_column<float>("uv",      niter);
                fitter_cache.allocate_column<float>("vj",      niter);
                fitter_cache.allocate_column<float>("ngal",    niter);
                fitter_cache.allocate_column<float>("e1_true", niter);
                fitter_cache.allocate_column<float>("e2_true", niter);
                fitter_cache.allocate_column<float>("r2_true", niter);

                fitter_cache.allocate_column<float>("chi2_obs",    niter, nmc);
                fitter_cache.allocate_column<float>("z_obs",       niter, nmc);
                fitter_cache.allocate_column<float>("z_obsm",      niter, nmc);
                fitter_cache.allocate_column<float>("e1_obs",      niter);
                fitter_cache.allocate_column<float>("e2_obs",      niter);
                fitter_cache.allocate_column<float>("r2_obs",      niter);
                fitter_cache.allocate_column<float>("e1_obsm",     niter);
                fitter_cache.allocate_column<float>("e2_obsm",     niter);
                fitter_cache.allocate_column<float>("r2_obsm",     niter);
                fitter_cache.allocate_column<float>("e1_obs_err",  niter);
                fitter_cache.allocate_column<float>("e2_obs_err",  niter);
                fitter_cache.allocate_column<float>("r2_obs_err",  niter);
                fitter_cache.allocate_column<float>("e1_obsm_err", niter);
                fitter_cache.allocate_column<float>("e2_obsm_err", niter);
                fitter_cache.allocate_column<float>("r2_obsm_err", niter);

                // Save meta data
                fitter_cache.write_columns("m_grid", m, "bt_grid", bt);
                fitter_cache.write_columns("bands", tbands, "lambda", lambda);

                // Let the fitter add its own data to the cache
                fitter.initialize_cache(fitter_cache, cache_mutex);

                fitter_cache.flush();

                note("done.");
            }

            // Initialize individual arrays
            if (keep_individuals_in_memory) {
                indiv_tr.resize(niter);
                indiv_id_mass.resize(niter);
                indiv_id_type.resize(niter);
                indiv_id_disk.resize(niter);
                indiv_id_bulge.resize(niter);
                indiv_id_bt.resize(niter);
                indiv_ngal.resize(niter);
                indiv_uv.resize(niter);
                indiv_vj.resize(niter);
                indiv_fdisk.resize(niter,nband+1);
                indiv_fbulge.resize(niter,nband+1);
                indiv_ml.resize(niter,nmc);
                indiv_ma.resize(niter,nmc);
                indiv_chi2.resize(niter,nmc);
                indiv_zml.resize(niter,nmc);
                indiv_zma.resize(niter,nmc);
            }

            // Compute averages at that redshift
            if (!global_progress_bar) {
                pgi = progress_start(niter);
            }
            generate(zf, dz);

            // Average quantities and store
            if (keep_averages_in_memory) {
                m_tr /= ngal;
                m_ml /= ngal;
                m_ma /= ngal;

                ztr[itz]  = m_tr;
                zml[itz]  = m_ml;
                zma[itz]  = m_ma;
                dndz[itz] = ngal/dz;
            }

            // Write to disk the individual measurements
            if (write_individuals) {
                std::string indiv_filename = cache_dir+"indiv-"+fitter.code_name+"-z"+zid+"-"+cache_id+".fits";
                note("individuals file: ", indiv_filename);

                file::mkdir(cache_dir);
                fits::write_table(indiv_filename,
                    "im",       indiv_id_mass,
                    "it",       indiv_id_type,
                    "ibt",      indiv_id_bt,
                    "idisk",    indiv_id_disk,
                    "ibulge",   indiv_id_bulge,
                    "ngal",     indiv_ngal,
                    "uv",       indiv_uv,
                    "vj",       indiv_vj,
                    "fdisk",    indiv_fdisk,
                    "fbulge",   indiv_fbulge,
                    "e1_true",  get_e1(indiv_tr),
                    "e2_true",  get_e2(indiv_tr),
                    "r2_true",  get_r2(indiv_tr),
                    "e1_obs",   get_e1(indiv_ml),
                    "e2_obs",   get_e2(indiv_ml),
                    "r2_obs",   get_r2(indiv_ml),
                    "e1_obsm",  get_e1(indiv_ma),
                    "e2_obsm",  get_e2(indiv_ma),
                    "r2_obsm",  get_r2(indiv_ma),
                    "chi2_obs", indiv_chi2,
                    "z_obs",    indiv_zml,
                    "z_obsm",   indiv_zma,
                    "bands", tbands, "lambda", lambda, "m_grid", m, "bt_grid", bt, "z_true", zf
                );

                fitter.save_individuals(indiv_filename);
            }

            if (global_progress_bar) {
                progress(pgi);
            }
        }

        // Average over N(z)
        if (keep_averages_in_memory) {
            double ntot = integrate(uzf, dndz);
            m_tr = integrate(uzf, ztr*(dndz/ntot));
            m_ml = integrate(uzf, zml*(dndz/ntot));
            m_ma = integrate(uzf, zma*(dndz/ntot));

            // Save
            if (write_averages) {
                std::string file_base = "psf-mean-z"+to_string(iz)+"-"+fitter.code_name;
                to_fits(file_base+"-tr.fits", m_tr, uzf, dndz, ztr);
                to_fits(file_base+"-ml.fits", m_ml, uzf, dndz, zml);
                to_fits(file_base+"-ma.fits", m_ma, uzf, dndz, zma);
            }
        }

        return true;
    }
};

#endif
