#include "egg-analytic.hpp"
#include "metrics.hpp"

struct mock_options {
    double dz = 0.01;
    uint_t nmc = 1000;
    uint_t seed = 42;
    vec1f depths;
    double min_mag_err = 0.05;
    bool no_noise = false;
    std::string psf_file;
};

class psf_averager : public egg::generator {
public :
    // Total number of galaxies
    double ngal = 0, nqu = 0, nsf = 0;

    // Averages
    metrics_set m_tr; // EGG (truth)

    // Monochromatic PSF library
    vec1d mono_lam, mono_q11, mono_q12, mono_q22, mono_w;

    // EGG PSF library
    vec2d egg_q11, egg_q12, egg_q22;
    vec2d egg_fu, egg_fv, egg_fj, egg_fvis;

    // Internal variables
    vec2d mc;
    vec1d uzf;
    vec1d dndz, dndz_qu, dndz_sf;

    bool just_count = false;
    uint_t niter = 0;

    progress_t pgi;

    // Config
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    double dz = 0.01;
    seed_t seed = make_seed(42);
    uint_t nband = npos, nmc = npos;
    vec1f phot_err2;
    double rel_err = dnan;
    std::string psf_file;
    filter_t psf_filter;

    bool single_pass = false;
    bool global_progress_bar = false;

    const std::string fitter;

    bool no_noise = false;

    // Cache
    std::string cache_filename;
    fits::table fitter_cache;
    bool write_cache = true;
    bool use_cache = true;
    bool cache_available = false;

    psf_averager(const std::string& f) : egg::generator(), fitter(f) {}

    void configure_mock(const mock_options& opts) {
        nband = filters.size();
        dz = opts.dz;
        no_noise = opts.no_noise;
        psf_file = opts.psf_file;

        // Square of photometric error (Gaussian additive component)
        phot_err2 = sqr(mag2uJy(opts.depths)/10.0);
        phypp_check(phot_err2.size() == nband, "mismatch between filters (", nband, ") and depths (",
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

        mc = randomn(seed, nmc, nband);

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
        phypp_check(read_filter("maiz-U",  rest_filter_u),
            "could not find rest-frame filter U, aborting");
        phypp_check(read_filter("maiz-V",  rest_filter_v),
            "could not find rest-frame filter V, aborting");
        phypp_check(read_filter("2mass-J", rest_filter_j),
            "could not find rest-frame filter J, aborting");

        egg_fu.resize(use.dims);
        egg_fv.resize(use.dims);
        egg_fj.resize(use.dims);
        for (uint_t iuv : range(use.dims[0]))
        for (uint_t ivj : range(use.dims[1])) {
            if (!use(iuv, ivj)) continue;

            vec1d tlam = lam(iuv,ivj,_);
            vec1d tsed = sed(iuv,ivj,_);

            tsed = lsun2uJy(0.0, 1.0, tlam, tsed);

            egg_fu(iuv, ivj) = sed2flux(rest_filter_u.lam, rest_filter_u.res, tlam, tsed);
            egg_fv(iuv, ivj) = sed2flux(rest_filter_v.lam, rest_filter_v.res, tlam, tsed);
            egg_fj(iuv, ivj) = sed2flux(rest_filter_j.lam, rest_filter_j.res, tlam, tsed);
        }
    }

    virtual void process_cached(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk,
        uint_t id_bulge, uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) = 0;

    virtual void do_fit(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk,
        uint_t id_bulge, uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) = 0;

    virtual void set_priors(const vec1d& fdisk, const vec1d& fbulge) = 0;

    void on_generated(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        const double ftot = fdisk.safe[0] + fbulge.safe[0];

        if (ftot >= flim) {
            // We passed the magnitude cut!

            if (just_count) {
                ++niter;
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

            m_tr.add(id_type, tngal*tr);

            ngal += tngal;
            if (id_type == 0) {
                nqu += tngal;
            } else {
                nsf += tngal;
            }

            // Now generate mocks and fit them
            // -------------------------------

            if (cache_available) {
                // A cache is present and we can reuse it!
                process_cached(id_mass, id_type, id_disk, id_bulge, id_bt, tngal, fdisk, fbulge);
            } else {
                // No cached data, must recompute stuff

                // Setup priors (if any)
                set_priors(fdisk, fbulge);

                // Do the fitting
                do_fit(id_mass, id_type, id_disk, id_bulge, id_bt, tngal, fdisk, fbulge);

                // Save things in the cache
                if (write_cache) {
                    // Compute actual UVJ colors
                    double mbt = bt.safe[id_bt];
                    double mbti = 1.0 - mbt;
                    double rfu = egg_fu.safe[id_disk]*mbti + mbt*egg_fu.safe[id_bulge];
                    double rfv = egg_fv.safe[id_disk]*mbti + mbt*egg_fv.safe[id_bulge];
                    double rfj = egg_fj.safe[id_disk]*mbti + mbt*egg_fj.safe[id_bulge];
                    double rfuv = -2.5*log10(rfu/rfv);
                    double rfvj = -2.5*log10(rfv/rfj);

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
                    fitter_cache.open(cache_filename);
                }
            }

            if (!global_progress_bar) {
                progress(pgi);
            }
        }
    }

    virtual void initialize_redshift_bin(uint_t iz) {}

    virtual void initialize_redshift_slice(uint_t itz) {}

    virtual void finalize_redshift_slice(uint_t itz) {}

    virtual void finalize_redshift_bin(uint_t it, double ntot, double ntot_qu, double ntot_sf) {}

    virtual std::string make_cache_hash() = 0;

    virtual void initialize_cache() {}

    bool average_redshift_bin(uint_t iz) {
        uint_t ntz = max(floor((zb[iz+1]-zb[iz])/dz), 1);
        uzf = zb[iz] + dz*dindgen(ntz);

        // Start averaging
        initialize_redshift_bin(iz);

        vec<1,metrics_set> ztr(ntz);
        dndz.resize(ntz);
        dndz_qu.resize(ntz);
        dndz_sf.resize(ntz);

        if (global_progress_bar) {
            pgi = progress_start(ntz);
        }

        for (uint_t itz : range(ntz)) {
            if (!global_progress_bar) {
                print(itz, "/", ntz);
            }

            double zf = uzf.safe[itz];

            // Compute PSF moments for each template in the library
            egg_q11.resize(use.dims);
            egg_q12.resize(use.dims);
            egg_q22.resize(use.dims);
            egg_fvis.resize(use.dims);

            for (uint_t iuv : range(use.dims[0]))
            for (uint_t ivj : range(use.dims[1])) {
                if (!use(iuv, ivj)) continue;

                vec1d tlam = lam(iuv,ivj,_);
                vec1d tsed = sed(iuv,ivj,_);

                if (!naive_igm) {
                    apply_madau_igm(zf, tlam, tsed);
                }

                tsed = lsun2uJy(zf, 1.0, tlam, tsed);
                tlam *= (1.0 + zf);

                double fvis       = sed2flux(psf_filter.lam, psf_filter.res, tlam, tsed);
                egg_q11(iuv,ivj)  = sed2flux(psf_filter.lam, psf_filter.res*mono_q11, tlam, tsed)/fvis;
                egg_q12(iuv,ivj)  = sed2flux(psf_filter.lam, psf_filter.res*mono_q12, tlam, tsed)/fvis;
                egg_q22(iuv,ivj)  = sed2flux(psf_filter.lam, psf_filter.res*mono_q22, tlam, tsed)/fvis;
                egg_fvis(iuv,ivj) = fvis;
            }

            // Reset averages
            m_tr.reset();
            ngal = 0; nqu = 0; nsf = 0;

            // Pre-compute number of iterations
            if (!single_pass) {
                note("compute number of iterations...");
                niter = 0;
                just_count = true;
                generate(zf, dz);
                note("done: ", niter);
            }

            // Initialize fitter
            initialize_redshift_slice(itz);

            // Initialize cache
            if (use_cache || write_cache) {
                std::string zid = replace(to_string(format::fixed(format::precision(zf, 2))), ".", "p");
                std::string cache_id = hash(make_cache_hash(), bands, phot_err2, nmc, niter);
                cache_filename = "cache-"+fitter+"-z"+zid+"-"+cache_id+".fits";
                note("cache file: ", cache_filename);
            }

            cache_available = false;
            if (use_cache && file::exists(cache_filename)) {
                cache_available = true;
                try {
                    fitter_cache.open(cache_filename);
                } catch (...) {
                    warning("cache file ", cache_filename, " was corrupted");
                    warning("will remove this cache and create a new one");

                    fitter_cache.close();
                    cache_available = false;
                }
            }

            if (!cache_available && write_cache) {
                note("creating cache on disk...");

                file::remove(cache_filename);
                fitter_cache.open(cache_filename);

                // Create cache arrays
                fitter_cache.allocate_column<uint_t>("im",           niter);
                fitter_cache.allocate_column<uint_t>("it",           niter);
                fitter_cache.allocate_column<uint_t>("idisk",        niter);
                fitter_cache.allocate_column<uint_t>("ibulge",       niter);
                fitter_cache.allocate_column<uint_t>("ibt",          niter);
                fitter_cache.allocate_column<float>("fbulge",        niter, nband+1);
                fitter_cache.allocate_column<float>("fdisk",         niter, nband+1);
                fitter_cache.allocate_column<float>("uv",            niter);
                fitter_cache.allocate_column<float>("vj",            niter);
                fitter_cache.allocate_column<float>("ngal",          niter);
                fitter_cache.allocate_column<float>("e1_true",       niter);
                fitter_cache.allocate_column<float>("e2_true",       niter);
                fitter_cache.allocate_column<float>("r2_true",       niter);

                // Save meta data
                fitter_cache.write_columns("m_grid", m, "bt_grid", bt);
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
                fitter_cache.write_columns("bands", tbands, "lambda", lambda);

                // Let the fitter add its own data to the cache
                initialize_cache();

                fitter_cache.flush();

                note("done.");
            }

            // Compute averages at that redshift
            iter = 0;
            just_count = false;
            if (!global_progress_bar) {
                pgi = progress_start(niter);
            }
            generate(zf, dz);

            // Average quantities and store
            m_tr.normalize(ngal, nqu, nsf);
            ztr[itz]     = m_tr;
            dndz[itz]    = ngal/dz;
            dndz_qu[itz] = nqu/dz;
            dndz_sf[itz] = nsf/dz;

            // Store results from fitter
            finalize_redshift_slice(itz);

            if (global_progress_bar) {
                progress(pgi);
            }
        }

        // Average over N(z)
        double ntot    = integrate(uzf, dndz);
        double ntot_qu = integrate(uzf, dndz_qu);
        double ntot_sf = integrate(uzf, dndz_sf);
        m_tr = integrate(uzf, ztr, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);

        // Save
        to_fits("psf-mean-z"+to_string(iz)+"-tr.fits", m_tr,
            uzf, dndz, dndz_qu, dndz_sf, ztr
        );

        // Save results from fitter
        finalize_redshift_bin(iz, ntot, ntot_qu, ntot_sf);

        return true;
    }
};
