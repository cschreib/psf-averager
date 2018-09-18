#include "egg-analytic.hpp"
#include "metrics.hpp"

struct mock_options {
    double dz = 0.01;
    uint_t nzsplit = 0;
    uint_t nmc = 200;
    uint_t seed = 42;
    vec1f depths;
    double min_mag_err = 0.05;
    bool no_noise = false;
    std::string psf_file;
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
    vec1d mono_lam, mono_q11, mono_q12, mono_q22, mono_w;

    // Observed SEDs
    vec1f save_sed_lambda;     // [nlam]
    vec2f save_sed_egg_fluxes; // [nsed,nlam]

    // EGG PSF library
    vec1d egg_q11, egg_q12, egg_q22;
    vec1d egg_fu, egg_fv, egg_fj, egg_fvis;

    // Internal variables
    vec2d mc;

    bool just_count = false;
    uint_t niter = 0;
    double ngal_tot = 0.0;

    std::mutex pg_mutex;
    progress_t pgi;

    // Config
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    double dz = 0.01;
    uint_t nzsplit = 0;
    seed_t seed = make_seed(42);
    uint_t nband = npos, nmc = npos;
    vec1f phot_err2;
    double rel_err = dnan;
    std::string psf_file;
    filter_t psf_filter;
    double prob_limit = 0.1;

    bool single_pass = false;
    bool global_progress_bar = false;

    const std::string fitter;

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

    psf_averager(const std::string& f) : egg::generator(), fitter(f) {}

    void configure_mock(const mock_options& opts) {
        nband = filters.size();
        dz = opts.dz;
        nzsplit = opts.nzsplit;
        zb = opts.zb;
        no_noise = opts.no_noise;
        psf_file = opts.psf_file;
        keep_individuals_in_memory = opts.keep_individuals_in_memory;
        write_individuals = opts.write_individuals;
        keep_averages_in_memory = opts.keep_averages_in_memory;
        write_averages = opts.write_averages;
        force_cache_id = opts.force_cache_id;
        cache_dir = file::directorize(opts.cache_dir);
        prob_limit = opts.prob_limit;
        if (prob_limit == 0.0) prob_limit = dnan;
        save_seds = opts.save_seds;
        global_progress_bar = opts.global_progress_bar;

        if (save_seds) {
            file::mkdir("seds/");
            save_sed_lambda = rgen_step(0.1,2.0,0.001);
        }

        if (write_individuals) {
            keep_individuals_in_memory = true;
        }
        if (write_averages) {
            keep_averages_in_memory = true;
        }

        // Square of photometric error (Gaussian additive component)
        phot_err2 = sqr(mag2uJy(opts.depths)/10.0);
        vif_check(phot_err2.size() == nband, "mismatch between filters (", nband, ") and depths (",
            phot_err2.size(), ")");

        // Relative error on flux, sets minimum uncertainties
        rel_err = opts.min_mag_err*(log(10.0)/2.5);

        // Cache random noise for re-use (same MC noise will be repeated for all galaxies)
        nmc = opts.nmc;
        seed = make_seed(opts.seed);

        if (no_noise && nmc != 1) {
            note("no_noise set, setting nmc=1");
            nmc = 1;
        }

        // Set PSF filter
        psf_filter = selection_filter;

        // Read monochromatic PSF library
        fits::read_table(psf_file,
            "lambda", mono_lam, "w", mono_w, "q11", mono_q11, "q12", mono_q12, "q22", mono_q22
        );

        // Match it to the PSF filter
        mono_w   = interpolate(mono_w,   mono_lam, selection_filter.lam);
        mono_q11 = interpolate(mono_q11, mono_lam, selection_filter.lam);
        mono_q12 = interpolate(mono_q12, mono_lam, selection_filter.lam);
        mono_q22 = interpolate(mono_q22, mono_lam, selection_filter.lam);
        mono_lam = selection_filter.lam;

        // Ignore data outside of adopted bandpass (450-950)
        mono_w[where(mono_lam < 0.450 || mono_lam > 0.950)] = 0.0;

        // Include PSF weighting flux loss in PSF filter response
        psf_filter.res *= mono_w;
        psf_filter.res /= integrate(psf_filter.lam, psf_filter.res);

        // Compute actual rest-frame colors of EGG SEDs.
        filter_t rest_filter_u, rest_filter_v, rest_filter_j;
        vif_check(read_filter("maiz-U",  rest_filter_u),
            "could not find rest-frame filter U, aborting");
        vif_check(read_filter("maiz-V",  rest_filter_v),
            "could not find rest-frame filter V, aborting");
        vif_check(read_filter("2mass-J", rest_filter_j),
            "could not find rest-frame filter J, aborting");

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

    struct fit_result {
        explicit fit_result(uint_t n) {
            chi2.resize(n);
            z_obs.resize(n);
            z_obsm.resize(n);
            psf_obs.resize(n);
            psf_obsm.resize(n);
        }

        vec1f chi2;
        vec1f z_obs, z_obsm;
        vec<1,metrics> psf_obs, psf_obsm;
    };

    virtual fit_result do_fit(uint_t iter, const vec1d& ftot) = 0;

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

        fit_result fr = do_fit(iter, ftot);

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

    virtual void initialize_fitter(double zf) {}

    virtual void save_individuals(const std::string& filename) {}

    virtual std::string make_cache_hash() = 0;

    virtual void initialize_cache(fits::table& cache) {}

    vec1d resample_sed(const vec1d& tlam, const vec1d& tsed) const {
        const vec1d& blam = save_sed_lambda;
        double dl = blam[1] - blam[0];
        uint_t npt = blam.size();
        vec1d fobs(npt);
        for (uint_t l : range(npt)) {
            if (blam.safe[l]-dl/2.0 > tlam.front() && blam.safe[l]+dl/2.0 < tlam.back()) {
                fobs.safe[l] = integrate(tlam, tsed, blam.safe[l]-dl/2.0, blam.safe[l]+dl/2.0)/dl;
            }
        }

        return fobs;
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

            // Compute PSF moments for each template in the library
            auto set_moments = [&](uint_t ised, vec1d tlam, vec1d tsed) {
                if (!naive_igm) {
                    apply_madau_igm(zf, tlam, tsed);
                }

                tsed = lsun2uJy(zf, df, tlam, tsed);
                tlam *= (1.0 + zf);

                if (save_seds) {
                    save_sed_egg_fluxes(ised,_) = ML_cor*resample_sed(tlam, tsed);
                }

                double fvis    = sed2flux(psf_filter.lam, psf_filter.res, tlam, tsed);
                egg_q11[ised]  = sed2flux(psf_filter.lam, psf_filter.res*mono_q11, tlam, tsed)/fvis;
                egg_q12[ised]  = sed2flux(psf_filter.lam, psf_filter.res*mono_q12, tlam, tsed)/fvis;
                egg_q22[ised]  = sed2flux(psf_filter.lam, psf_filter.res*mono_q22, tlam, tsed)/fvis;
                egg_fvis[ised] = fvis;
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

            // Initialize fitter
            initialize_fitter(zf);

            // Initialize cache
            uint_t precision = ceil(max(-log10(dz), 2.0));
            std::string zid = replace(to_string(format::fixed(format::precision(zf, precision))), ".", "p");
            std::string cache_id;
            if (force_cache_id.empty()) {
                cache_id = hash(make_cache_hash(), bands, phot_err2, nmc, niter);
            } else {
                cache_id = force_cache_id;
            }

            if (write_cache) {
                cache_filename = cache_dir+"cache-"+fitter+"-z"+zid+"-"+cache_id+".fits";
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
                initialize_cache(fitter_cache);

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

            // Initialize random numbers
            mc = randomn(seed, nmc, nband);
            // Subtract mean noise per band to ensure random photometry is unbiased
            for (uint_t i : range(nband)) {
                mc(i,_) -= mean(mc(i,_));
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
                std::string indiv_filename = cache_dir+"indiv-"+fitter+"-z"+zid+"-"+cache_id+".fits";
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

                save_individuals(indiv_filename);
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
                to_fits("psf-mean-z"+to_string(iz)+"-"+fitter+"-tr.fits", m_tr, uzf, dndz, ztr);
                to_fits("psf-mean-z"+to_string(iz)+"-"+fitter+"-ml.fits", m_ml, uzf, dndz, zml);
                to_fits("psf-mean-z"+to_string(iz)+"-"+fitter+"-ma.fits", m_ma, uzf, dndz, zma);
            }
        }

        return true;
    }
};
