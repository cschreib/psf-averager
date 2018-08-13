#include "egg-analytic.hpp"
#include "metrics.hpp"

struct mock_options {
    double dz = 0.01;
    uint_t nmc = 1000;
    uint_t seed = 42;
    vec1f depths;
    double min_mag_err = 0.05;
    bool no_noise = false;
};

class psf_averager : public egg::generator {
public :
    // Total number of galaxies
    double ngal = 0, nqu = 0, nsf = 0;

    // Averages
    metrics_set m_tr; // EGG (truth)

    // EGG PSF library
    vec1s egg_ssed, egg_zz;
    vec1s egg_sz;
    vec1d egg_z;
    vec1d egg_q11, egg_q12, egg_q22;
    vec2d eggz_q11, eggz_q12, eggz_q22;

    // Internal variables
    vec2d mc;
    vec1d uzf;
    vec1d dndz, dndz_qu, dndz_sf;

    bool just_count = false;
    uint_t niter = 0;
    uint_t iter = 0;

    progress_t pgi;

    // Config
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    double dz = 0.01;
    seed_t seed = make_seed(42);
    uint_t nband = npos, nmc = npos;
    vec1f phot_err2;
    double rel_err = dnan;

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

        // Square of photometric error (Gaussian additive component)
        phot_err2 = sqr(mag2uJy(opts.depths)/10.0);
        phypp_check(phot_err2.size() == nband, "mismatch between filters (", nband, ") and depths (",
            phot_err2.size(), ")");

        // Relative error on flux, sets minimum uncertainties
        rel_err = opts.min_mag_err*(log(10.0)/2.5);

        // Cache random noise for re-use (same MC noise will be repeated for all galaxies)
        nmc = opts.nmc;
        seed = make_seed(opts.seed);

        if (no_noise) {
            note("no_noise set, setting nmc=1");
            nmc = 1;
        }

        mc = randomn(seed, nmc, nband);
    }

    virtual void process_cached(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) = 0;

    virtual void do_fit(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) = 0;

    virtual void set_priors(const vec1d& fdisk, const vec1d& fbulge) = 0;

    void on_generated(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
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
            const double fbtn = fbulge.safe[0]/ftot;
            const double fbti = 1.0 - fbtn;

            // Add up PSFs
            metrics tr(
                eggz_q11.safe[id_disk]*fbti + eggz_q11.safe[id_bulge]*fbtn,
                eggz_q12.safe[id_disk]*fbti + eggz_q12.safe[id_bulge]*fbtn,
                eggz_q22.safe[id_disk]*fbti + eggz_q22.safe[id_bulge]*fbtn
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
                    fitter_cache.update_elements("im",        id_mass,      fits::at(iter));
                    fitter_cache.update_elements("it",        id_type,      fits::at(iter));
                    fitter_cache.update_elements("idisk",     id_disk,      fits::at(iter));
                    fitter_cache.update_elements("ibulge",    id_bulge,     fits::at(iter));
                    fitter_cache.update_elements("ibt",       id_bt,        fits::at(iter));
                    fitter_cache.update_elements("ngal",      tngal,        fits::at(iter));
                    fitter_cache.update_elements("fbulge",    fbulge,       fits::at(iter,_));
                    fitter_cache.update_elements("fdisk",     fdisk,        fits::at(iter,_));
                    fitter_cache.update_elements("e1_true",   tr.e1,        fits::at(iter));
                    fitter_cache.update_elements("e2_true",   tr.e2,        fits::at(iter));
                    fitter_cache.update_elements("r2_true",   tr.r2,        fits::at(iter));
                    fitter_cache.open(cache_filename);
                }
            }

            ++iter;
            progress(pgi);
        }
    }

    bool read_egg_psfs(uint_t iz) {
        std::string zdir = "full_z"+to_string(iz)+"/";
        std::string filename = zdir+"psfs-rebin2-cst.txt";
        if (!file::exists(filename)) return false;

        // Read PSF library
        ascii::read_table(filename, egg_ssed, egg_zz, _, _, _, egg_q11, egg_q12, egg_q22);

        phypp_check(!egg_zz.empty(), "empty PSF file '", filename, "'");

        // Get SED ID
        // Make sure SED ID is of the form xx-yy (i.e. convert 9 to 09)
        for (uint_t i : range(egg_ssed)) {
            vec1s spl = split(egg_ssed.safe[i], "-");
            if (spl[0].size() == 1) {
                spl[0] = '0'+spl[0];
            }
            if (spl[1].size() == 1) {
                spl[1] = '0'+spl[1];
            }

            egg_ssed.safe[i] = spl[0]+'-'+spl[1];
        }

        // Sort by z then SED
        {
            vec1u ids = sort(egg_zz+egg_ssed);
            egg_ssed  = egg_ssed[ids];
            egg_zz    = egg_zz[ids];
            egg_q11   = egg_q11[ids];
            egg_q12   = egg_q12[ids];
            egg_q22   = egg_q22[ids];
        }

        // Find redshifts
        egg_sz = unique_values_sorted(egg_zz);
        from_string(replace(egg_sz, "p", "."), egg_z);

        return true;
    }

    virtual void initialize_redshift_bin(uint_t iz) {}

    virtual void initialize_redshift_slice(uint_t itz) {}

    virtual void finalize_redshift_slice(uint_t itz) {}

    virtual void finalize_redshift_bin(uint_t it, double ntot, double ntot_qu, double ntot_sf) {}

    virtual std::string make_cache_hash() = 0;

    virtual void initialize_cache() {}

    bool average_redshift_bin(uint_t iz) {
        if (!read_egg_psfs(iz)) return false;

        uint_t ntz = max(floor((zb[iz+1]-zb[iz])/dz), 1);
        uzf = zb[iz] + dz*dindgen(ntz);

        // Start averaging
        initialize_redshift_bin(iz);

        vec<1,metrics_set> ztr(ntz);
        dndz.resize(ntz);
        dndz_qu.resize(ntz);
        dndz_sf.resize(ntz);

        for (uint_t itz : range(ntz)) {
            print(itz, "/", ntz);

            double zf = uzf.safe[itz];

            // Pre-select redshift slice in EGG PSF library
            vec1u idz = where(egg_zz == egg_sz[min_id(abs(egg_z - zf))]);
            eggz_q11.resize(use.dims);
            eggz_q12.resize(use.dims);
            eggz_q22.resize(use.dims);

            for (uint_t i : idz) {
                vec1s spl = split(egg_ssed.safe[i], "-");

                uint_t iuv, ivj;
                from_string(spl[0], iuv);
                from_string(spl[1], ivj);

                eggz_q11.safe(iuv,ivj) = egg_q11.safe[i];
                eggz_q12.safe(iuv,ivj) = egg_q12.safe[i];
                eggz_q22.safe(iuv,ivj) = egg_q22.safe[i];
            }

            // Reset averages
            m_tr.reset();
            ngal = 0; nqu = 0; nsf = 0;

            // Pre-compute number of iterations
            note("compute number of iterations...");
            niter = 0;
            just_count = true;
            generate(zf, dz);
            note("done: ", niter);

            // Initialize fitter
            initialize_redshift_slice(itz);

            // Initialize cache
            std::string zid = replace(to_string(format::fixed(format::precision(zf, 2))), ".", "p");
            std::string cache_id = hash(make_cache_hash(), bands, phot_err2, nmc, niter);
            cache_filename = "cache-"+fitter+"-z"+zid+"-"+cache_id+".fits";
            note("cache file: ", cache_filename);

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
                tbands[1-_] = bands;
                for (uint_t l : range(filters)) {
                    lambda[l+1] = filters[l].rlam;
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
            pgi = progress_start(niter);
            generate(zf, dz);

            // Average quantities and store
            m_tr.normalize(ngal, nqu, nsf);
            ztr[itz]     = m_tr;
            dndz[itz]    = ngal/dz;
            dndz_qu[itz] = nqu/dz;
            dndz_sf[itz] = nsf/dz;

            // Store results from fitter
            finalize_redshift_slice(iz);
        }

        // Average over N(z)
        double ntot    = integrate(uzf, dndz);
        double ntot_qu = integrate(uzf, dndz_qu);
        double ntot_sf = integrate(uzf, dndz_sf);
        m_tr = integrate(uzf, ztr, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);

        // Save
        to_fits("psf-mean-z"+to_string(iz)+"-bpz-tr.fits", m_tr,
            uzf, dndz, dndz_qu, dndz_sf, ztr
        );

        // Save results from fitter
        finalize_redshift_bin(iz, ntot, ntot_qu, ntot_sf);

        return true;
    }
};
