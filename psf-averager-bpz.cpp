#include "egg-analytic.hpp"
#include "metrics.hpp"

struct fitter_options {
    uint_t nmc = 1000;
    uint_t seed = 42;
    vec1f depths;
    std::string prior_filter;
    double min_mag_err = 0.05;
    double gauss_convolve = 0.021; // 0.03/sqrt(2) to match default BPZ
    uint_t ninterp = 2;
    double dz = 0.01;
    double zfit_max = 7.0;
    double zfit_dz = 0.01;
    bool use_capak_library = true;
    bool apply_igm = true;
    bool force_true_z = false;
    bool no_noise = false;
};

class psf_averager : public egg::generator {
public :
    // Total number of galaxies
    double ngal = 0, nqu = 0, nsf = 0;

    // Averages
    metrics_set m_ml; // maximum likelihood
    metrics_set m_ma; // marginalization
    metrics_set m_tr; // EGG (truth)

    // BPZ PSF library
    vec1s bpz_seds, bpz_seds_sid;
    vec1s bpz_ssed, bpz_zz;
    vec1f bpz_z;
    vec1d bpz_q11, bpz_q12, bpz_q22;
    vec1d zfit_base;
    uint_t nzfit_base = npos;
    uint_t nmodel_base = npos;

    // Subsample of BPZ PSF library
    vec1d bpz_q11_fit, bpz_q12_fit, bpz_q22_fit;
    vec1d zfit;
    uint_t nmodel = npos;
    uint_t nzfit = npos;

    // EGG PSF library
    vec1s egg_ssed, egg_zz;
    vec1s egg_sz;
    vec1d egg_z;
    vec1d egg_q11, egg_q12, egg_q22;
    vec2d eggz_q11, eggz_q12, eggz_q22;

    // Internal variables
    vec1d wflux, rflux, weight;
    vec2d wmodel;
    vec2d mc;
    vec2d tpl_flux;
    uint_t ntemplate = npos;
    uint_t nelliptical = npos;
    uint_t nspiral = npos;
    uint_t nirregular = npos;
    vec1d pmodel;
    vec1d prior;
    vec1d tt_prior;
    vec1d pz, pzc;
    vec1d wmm;

    vec1d zmeas_ml, zmeas_ma, cache_pmodel;
    vec1u cache_bmodel;

    bool just_count = false;
    uint_t niter = 0;
    uint_t iter = 0;

    progress_t pgi;

    // Config
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    double dz = 0.01;
    uint_t nmc = 1000;
    seed_t seed = make_seed(42);
    uint_t nband = npos;
    vec1f phot_err2;
    double rel_err = dnan;
    double gauss_convolve = dnan;
    double dzfit = dnan;
    vec1d gconv;
    uint_t wconv = npos;
    uint_t id_prior = npos;
    uint_t ninterp = 0;
    bool write_cache = true;
    bool use_cache = true;
    bool cache_available = false;
    bool use_capak_library = true;
    bool apply_igm = true;
    bool force_true_z = false;
    bool no_noise = false;

    // Cache
    std::string cache_filename;
    fits::table zmeas_cache;
    uint_t cache_flush_rate = 10;

    psf_averager() : egg::generator() {}

    void configure_fitter(const fitter_options& opts) {
        dz = opts.dz;
        nband = filters.size();
        use_capak_library = opts.use_capak_library;
        apply_igm = opts.apply_igm;
        force_true_z = opts.force_true_z;
        no_noise = opts.no_noise;

        phypp_check(is_any_of(opts.prior_filter, bands),
            "prior filter is not in the filter list ('", opts.prior_filter, "' not found')");
        id_prior = where_first(opts.prior_filter == bands) + 1;

        // Square of photometric error (Gaussian additive component)
        phot_err2 = sqr(mag2uJy(opts.depths)/10.0);
        phypp_check(phot_err2.size() == nband, "mismatch between filters (", nband, ") and depths (",
            phot_err2.size(), ")");

        // Relative error on flux, sets minimum uncertainties
        rel_err = opts.min_mag_err*(log(10.0)/2.5);

        // Cache random noise for re-use (same MC noise will be repeated for all galaxies)
        nmc = opts.nmc;
        seed = make_seed(opts.seed);
        mc = randomn(seed, nmc, nband);

        if (no_noise) {
            note("no_noise set, setting nmc=1");
            nmc = 1;
        }

        // Kernel used to convolve the P(z) to wash out too small fluctuations
        gauss_convolve = opts.gauss_convolve;

        // Interpolation of BPZ SEDs
        ninterp = opts.ninterp;

        // Read BPZ PSF library
        std::string bpz_filename;
        if (use_capak_library) {
            bpz_filename = "BPZ_capak/psfs-rebin2-cst.txt";
        } else {
            bpz_filename = "BPZ/psfs-rebin2-cst.txt";
        }

        ascii::read_table(bpz_filename, bpz_ssed, bpz_zz, _, _, _, bpz_q11, bpz_q12, bpz_q22);

        // Resample BPZ library
        resample_library(opts);

        // Make convolver kernel
        {
            uint_t npt = 18*gauss_convolve/dzfit;
            if (npt % 2 == 0) ++npt;
            wconv = npt/2;

            // This is a gaussian trucated between -9*sigma and +9*sigma
            vec1d kz = dzfit*(findgen(npt)-wconv);
            gconv = integrate_gauss(kz-dzfit/2.0, kz+dzfit/2.0, 0.0, gauss_convolve);
            gconv /= total(gconv);
        }
    }

    void resample_library(const fitter_options& opts) {
        // Sort by z then SED
        {
            vec1u ids = sort(bpz_zz+bpz_ssed);
            bpz_ssed  = bpz_ssed[ids];
            bpz_zz    = bpz_zz[ids];
            bpz_q11   = bpz_q11[ids];
            bpz_q12   = bpz_q12[ids];
            bpz_q22   = bpz_q22[ids];
        }

        // Reduce the size of the BPZ z grid
        vec1s ubz = unique_values_sorted(bpz_zz);
        vec1f tz;
        from_string(replace(ubz, "p", "."), tz);

        if (!force_true_z) {
            // We only do this if not fitting at the true redshift

            uint_t skip_every = max(1, round(opts.zfit_dz/0.001));
            vec1u idr;
            uint_t k = 0;
            vec1s oubz = ubz;
            ubz.clear();
            for (uint_t i : range(oubz)) {
                if (k % skip_every != 0 || tz[i] > opts.zfit_max) {
                    append(idr, where(bpz_zz == oubz[i]));
                } else {
                    ubz.push_back(oubz[i]);
                }
                ++k;
            }

            uint_t old_size = bpz_ssed.size();
            inplace_remove(bpz_ssed, idr);
            inplace_remove(bpz_zz,   idr);
            inplace_remove(bpz_q11,  idr);
            inplace_remove(bpz_q12,  idr);
            inplace_remove(bpz_q22,  idr);

            note("shrunk BPZ PSF library from ", old_size, " to ", bpz_ssed.size(), " (",
                ubz.size(), " redshifts, ", bpz_ssed.size()/ubz.size(), " SEDs)");
        }

        // Find redshifts
        from_string(replace(ubz, "p", "."), zfit_base);
        nzfit_base = zfit_base.size();

        dzfit = zfit_base[1] - zfit_base[0];
        note("redshift step: ", dzfit);

        // Sort BPZ SEDs by color (red to blue)
        std::string sed_dir;
        if (use_capak_library) {
            sed_dir = "/home/cschreib/programming/bpz-1.99.3/SED_capak/";
        } else {
            sed_dir = "/home/cschreib/programming/bpz-1.99.3/SED/";
        }

        bpz_seds = file::list_files(sed_dir, "*.sed");

        {
            vec1d color(bpz_seds.size());
            for (uint_t t : range(bpz_seds)) {
                vec1d rlam, rsed;
                ascii::read_table(sed_dir+bpz_seds[t], rlam, rsed);
                rsed *= 1e-19;
                rsed = cgs2uJy(rlam, rsed);

                vec1u id1 = where(rlam > 2000 && rlam < 3000);
                uint_t id2 = where_first(rlam > 10000);
                color[t] = median(rsed[id1])/rsed[id2];
            }

            vec1u ids = sort(color);
            bpz_seds = bpz_seds[ids];
        }

        ntemplate = bpz_seds.size();
        nelliptical = 1;
        nspiral = 2;
        nirregular = ntemplate - nelliptical - nspiral;

        nmodel_base = ntemplate*nzfit_base;

        for (uint_t i : range(bpz_seds)) {
            note(" - ", i, ": ", bpz_seds[i]);
        }

        phypp_check(nmodel_base == bpz_ssed.size(),
            "mismatch between BPZ SEDs in directory (", ntemplate, ") and PSF library (",
            nzfit_base, "x", bpz_ssed.size()/nzfit_base, ")");

        // Find correspondence with PSF library, which was sorted by file name...
        {
            vec1s psf_seds = bpz_seds;
            inplace_sort(psf_seds);

            vec1u id_translate(psf_seds.size());
            for (uint_t it : range(ntemplate)) {
                id_translate[it] = where_first(psf_seds == bpz_seds[it]);
            }

            vec1u ids(nmodel_base);
            for (uint_t izf : range(nzfit_base))
            for (uint_t it : range(ntemplate)) {
                ids[izf*ntemplate + it] = where_first(
                    bpz_zz == ubz[izf] && bpz_ssed == to_string(id_translate[it])
                );
            }

            // Reshuffle
            // bpz_ssed is left alone, because we want it to remain sorted and map to
            // SED IDs in our new SED ordering, not the one of the PSF library...
            bpz_zz  = bpz_zz[ids];
            bpz_q11 = bpz_q11[ids];
            bpz_q12 = bpz_q12[ids];
            bpz_q22 = bpz_q22[ids];
        }

        bpz_seds = sed_dir+bpz_seds;
        bpz_seds_sid = to_string_vector(indgen(ntemplate));

        // Interpolate the templates
        if (ninterp > 0) {
            std::string tmp_dir = "BPZ_interp/";
            file::mkdir(tmp_dir);

            // Save old PSF library
            vec1s oseds = bpz_seds;
            vec1s ossed = bpz_ssed;
            vec1d oq11  = bpz_q11;
            vec1d oq12  = bpz_q12;
            vec1d oq22  = bpz_q22;

            // Clear library
            bpz_seds.clear();
            bpz_seds_sid.clear();
            bpz_ssed.clear();
            bpz_zz.clear();
            bpz_q11.clear();
            bpz_q12.clear();
            bpz_q22.clear();

            // Add first SED
            bpz_seds.push_back(oseds[0]);
            bpz_seds_sid.push_back("0");
            vec1u id0 = where(ossed == bpz_seds_sid.back());
            append(bpz_ssed, replicate(bpz_seds_sid.back(), nzfit_base));
            append(bpz_zz,   ubz);
            append(bpz_q11,  oq11[id0]);
            append(bpz_q12,  oq12[id0]);
            append(bpz_q22,  oq22[id0]);

            uint_t ntemplate_old = ntemplate;
            ntemplate = 1;

            for (uint_t it : range(ntemplate_old-1)) {
                // Read two adjacent SEDs
                vec1d rlam1, rsed1, rlam2, rsed2;
                ascii::read_table(oseds[it],   rlam1, rsed1);
                ascii::read_table(oseds[it+1], rlam2, rsed2);

                // Convert to uJy (BPZ interpolates fluxes in f_nu)
                rsed1 = cgs2uJy(rlam1, rsed1*1e-19);
                rsed2 = cgs2uJy(rlam2, rsed2*1e-19);

                // Merge wavelength axes in case they don't match and interpolate on merged axis
                vec1d clam = rlam1;
                append(clam, rlam2);
                inplace_sort(clam);
                clam = unique_values_sorted(clam);
                rsed1 = interpolate(rsed1, rlam1, clam);
                rsed2 = interpolate(rsed2, rlam2, clam);

                vec1u id1 = where(ossed == to_string(it+1));

                // Interpolate
                for (uint_t ii : range(ninterp)) {
                    double x = (ii+1.0)/(ninterp+1.0);
                    vec1d tsed1 = (1.0-x)*rsed1;
                    vec1d tsed2 = x*rsed2;
                    vec1d rsed = tsed1 + tsed2;
                    rsed = uJy2cgs(clam, rsed)*1e19;

                    std::string fname = file::remove_extension(file::get_basename(oseds[it]))
                        +"_"+to_string(ii)+"_"+to_string(ninterp)+".sed";
                    ascii::write_table(tmp_dir+fname, clam, rsed);

                    // Insert SED in PSF library
                    bpz_seds.push_back(tmp_dir+fname);
                    std::string istr = align_right(to_string(ii), floor(log10(ninterp)), '0');
                    bpz_seds_sid.push_back(to_string(it)+"."+istr);
                    append(bpz_ssed, replicate(bpz_seds_sid.back(), nzfit_base));
                    append(bpz_zz,   ubz);
                    ++ntemplate;

                    for (uint_t iz : range(nzfit_base)) {
                        // We have to interpolate PSFs with a flux-weighted factor inside the VIS passband
                        vec1d tlam = clam*(1e-4*(1.0 + zfit_base[iz]));
                        double f1 = sed2flux(selection_filter.lam, selection_filter.res, tlam, tsed1);
                        double f2 = sed2flux(selection_filter.lam, selection_filter.res, tlam, tsed2);
                        double xf = f2/(f1 + f2);

                        bpz_q11.push_back(oq11[id0[iz]]*(1.0-xf) + xf*oq11[id1[iz]]);
                        bpz_q12.push_back(oq12[id0[iz]]*(1.0-xf) + xf*oq12[id1[iz]]);
                        bpz_q22.push_back(oq22[id0[iz]]*(1.0-xf) + xf*oq22[id1[iz]]);
                    }
                }

                // Added end SED
                bpz_seds.push_back(oseds[it+1]);
                bpz_seds_sid.push_back(to_string(it+1));
                append(bpz_ssed, replicate(bpz_seds_sid.back(), nzfit_base));
                append(bpz_zz,   ubz);
                append(bpz_q11,  oq11[id1]);
                append(bpz_q12,  oq12[id1]);
                append(bpz_q22,  oq22[id1]);
                ++ntemplate;

                std::swap(id0, id1);
            }

            // Re-sort library, presently sorted by SED and z, into z and SED.
            vec1u ids = sort(bpz_zz+bpz_ssed);
            bpz_ssed  = bpz_ssed[ids];
            bpz_zz    = bpz_zz[ids];
            bpz_q11   = bpz_q11[ids];
            bpz_q12   = bpz_q12[ids];
            bpz_q22   = bpz_q22[ids];

            nmodel_base = nzfit_base*ntemplate;

            note("expanded PSF library from ", ntemplate_old, " to ", ntemplate, " SEDs");
        }

        // Finalize stuff
        from_string(replace(bpz_zz, "p", "."), bpz_z);
    }

    void set_priors(double mag) {
        double momin;
        vec1d alpha, z0, km, kt, ft;

        if (use_capak_library) {
            // Raichoor et al. NGVS prior

            if (mag > 32.0) mag = 32.0;

            if (13.0 <= mag && mag < 17.0) {
                momin = 12.5;
                alpha = {2.69, 2.19, 1.99};
                z0    = {0.005, 0.004, 0.003};
                km    = {0.0256, 0.0200, 0.018};
                kt    = {0.030, -0.048};
                ft    = {0.52, 0.135};
            } else if (17.0 <= mag && mag < 20.0) {
                momin = 17.0;
                alpha = {2.58, 2.00, 2.00};
                z0    = {0.122, 0.094, 0.084};
                km    = {0.103, 0.099, 0.072};
                kt    = {0.138, -0.011};
                ft    = {0.45, 0.17};
            } else {
                momin = 20.0;
                alpha = {2.46, 1.81, 2.00};
                z0    = {0.431, 0.39, 0.3};
                km    = {0.091, 0.10, 0.15};
                kt    = {0.40, 0.3};
                ft    = {0.30, 0.175};
            }
        } else {
            // Benitez 2000

            mag = clamp(mag, 20.0, 32.0);

            momin = 20.0;
            alpha = {2.46, 1.81, 0.91};
            z0    = {0.431, 0.390, 0.063};
            km    = {0.091, 0.0636, 0.123};
            kt    = {0.147, 0.450};
            ft    = {0.35, 0.50};
        }

        if (use_capak_library && mag < 13.0) {
            // Super bright
            for (uint_t iz : range(nzfit)) {
                double val = (zfit.safe[iz] <= 0.1 || (zfit.safe[0] > 0.1 && iz == 0));
                for (uint_t it : range(3)) {
                    uint_t im = iz*3 + it;
                    tt_prior.safe[im] = val;
                }
            }
        } else {
            // Create prior for three SED classes
            vec1d ptype(3);
            ft /= vec1u{nelliptical, nspiral};
            ptype.safe[0] = ft.safe[0]*exp(-kt.safe[0]*(mag - momin));
            ptype.safe[1] = ft.safe[1]*exp(-kt.safe[1]*(mag - momin));
            ptype.safe[2] = (1.0 - nelliptical*ptype.safe[0] - nspiral*ptype.safe[1])/nirregular;

            for (uint_t tt : range(3)) {
                // For each SED class, evaluate function of redshift
                double zm = z0.safe[tt] + km.safe[tt]*(mag - momin);
                double tprior = 0.0;
                for (uint_t iz : range(nzfit)) {
                    uint_t im = iz*3 + tt;
                    tt_prior.safe[im] = pow(zfit.safe[iz], alpha.safe[tt])
                        *exp(-clamp(pow(zfit.safe[iz]/zm, alpha.safe[tt]), 0.0, 700.0));

                    tprior += tt_prior.safe[im];
                }

                double tprior2 = 0.0;
                for (uint_t iz : range(nzfit)) {
                    uint_t im = iz*3 + tt;

                    // Normalize
                    tt_prior.safe[im] /= tprior;

                    if (!use_capak_library) {
                        // Clip low probability wings
                        if (tt_prior.safe[im] < 1e-2/nzfit) {
                            tt_prior.safe[im] = 0.0;
                        } else {
                            tprior2 += tt_prior.safe[im];
                        }
                    }
                }

                if (!use_capak_library) {
                    for (uint_t iz : range(nzfit)) {
                        uint_t im = iz*3 + tt;

                        // Normalize again
                        tt_prior.safe[im] /= tprior2;
                    }
                }
            }
        }

        // Expand to SED library
        for (uint_t iz : range(nzfit))
        for (uint_t it : range(ntemplate)) {
            uint_t im = iz*ntemplate + it;

            uint_t it0 = it/(ninterp+1);
            uint_t tt0 = (it0 < nelliptical ? 0 : it0 < nelliptical+nspiral ? 1 : 2);
            uint_t it1 = it0 + 1;
            uint_t tt1 = (it1 < nelliptical ? 0 : it1 < nelliptical+nspiral ? 1 : 2);
            uint_t itm = it % (ninterp+1);

            if (itm == 0 || tt0 == tt1) {
                // Non-interpolated SEDs get priors from their class,
                // and so do SEDs inteprolated in between two SEDs of the same class
                prior.safe[im] = tt_prior.safe[iz*3 + tt0];
            } else {
                // SEDs interpolated between different classes get interpolated priors
                double x = (itm+1.0)/(ninterp+1.0);
                prior.safe[im] = tt_prior.safe[iz*3 + tt0]*(1.0-x) + x*tt_prior.safe[iz*3 + tt1];
            }
        }
    }

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

            // Now do BPZ
            // ----------

            if (cache_available) {
                // A cache is present and we can reuse it!
                zmeas_cache.read_elements("bmodel", cache_bmodel, fits::at(iter,_));
                zmeas_cache.read_elements("pmodel", cache_pmodel, fits::at(iter,_));

                // Maximum likelihood
                // ------------------

                metrics ml;
                metrics ml2;
                for (uint_t i : range(nmc)) {
                    uint_t iml = cache_bmodel.safe(i);
                    metrics tm(bpz_q11_fit.safe[iml], bpz_q12_fit.safe[iml], bpz_q22_fit.safe[iml]);
                    ml  += tm;
                    ml2 += sqr(tm);
                }
                ml /= nmc;
                ml2 /= nmc;
                ml2 = sqrt(ml2 - sqr(ml));

                m_ml.add(id_type, tngal*ml);

                // Marginalization
                // ---------------

                metrics ma;
                metrics ma2;
                for (uint_t im : range(nmodel)) {
                    metrics tm(bpz_q11_fit.safe[im], bpz_q12_fit.safe[im], bpz_q22_fit.safe[im]);
                    ma  += cache_pmodel.safe(im)*tm;
                    ma2 += cache_pmodel.safe(im)*sqr(tm);
                }
                ma2 = sqrt(ma2 - sqr(ma));

                m_ma.add(id_type, tngal*ma);

                zmeas_cache.update_elements("e1_bpz_ml", ml.e1,      fits::at(iter));
                zmeas_cache.update_elements("e2_bpz_ml", ml.e2,      fits::at(iter));
                zmeas_cache.update_elements("r2_bpz_ml", ml.r2,      fits::at(iter));
                zmeas_cache.update_elements("e1_bpz_ma", ma.e1,      fits::at(iter));
                zmeas_cache.update_elements("e2_bpz_ma", ma.e2,      fits::at(iter));
                zmeas_cache.update_elements("r2_bpz_ma", ma.r2,      fits::at(iter));
                zmeas_cache.update_elements("e1_bpz_ml_err", ml2.e1, fits::at(iter));
                zmeas_cache.update_elements("e2_bpz_ml_err", ml2.e2, fits::at(iter));
                zmeas_cache.update_elements("r2_bpz_ml_err", ml2.r2, fits::at(iter));
                zmeas_cache.update_elements("e1_bpz_ma_err", ma2.e1, fits::at(iter));
                zmeas_cache.update_elements("e2_bpz_ma_err", ma2.e2, fits::at(iter));
                zmeas_cache.update_elements("r2_bpz_ma_err", ma2.r2, fits::at(iter));
            } else {
                // No cached data, must recompute stuff

                // Create noise-free photometry
                for (uint_t l : range(nband)) {
                    double bftot = fdisk.safe[l+1] + fbulge.safe[l+1];
                    weight.safe[l] = 1.0/sqrt(phot_err2.safe[l] + sqr(bftot*rel_err));
                    rflux.safe[l] = wflux.safe[l] = bftot*weight.safe[l];
                }

                // Create models
                for (uint_t im : range(nmodel)) {
                    wmm.safe[im] = 0.0;
                    for (uint_t l : range(nband)) {
                        wmodel.safe(im,l) = tpl_flux.safe(im,l)*weight.safe[l];
                        wmm.safe[im] += sqr(wmodel.safe(im,l));
                    }
                }

                // Create priors
                double pmag = -2.5*log10(fdisk.safe[id_prior] + fbulge.safe[id_prior]) + 23.9;
                set_priors(pmag);

                if (write_cache) {
                    cache_pmodel[_] = 0.0;
                }

                // Simulate redshift measurements
                metrics ma, ml, ma2, ml2;
                for (uint_t i : range(nmc)) {
                    // Create noisy photometry
                    if (!no_noise) {
                        for (uint_t l : range(nband)) {
                            // In weighted units, the random perturbations have a sigma of unity
                            rflux.safe[l] = wflux.safe[l] + mc.safe(i,l);
                        }
                    }

                    // Find best SED and z, and store chi2 to build p(z,SED)
                    uint_t iml = npos;
                    double chi2_best = finf;
                    for (uint_t im : range(nmodel)) {
                        double* wm = &wmodel.safe(im,0);
                        double* wf = &rflux.safe[0];

                        double wfm = 0.0;
                        for (uint_t l : range(nband)) {
                            wfm += wf[l]*wm[l];
                        }

                        double scale = wfm/wmm.safe[im];

                        double tchi2 = 0.0;
                        for (uint_t l : range(nband)) {
                            tchi2 += sqr(wf[l] - scale*wm[l]);
                        }

                        pmodel.safe[im] = tchi2;

                        if (tchi2 < chi2_best) {
                            chi2_best = tchi2;
                            iml = im;
                        }
                    }

                    // Create likelihood, apply prior, and marginalize over SEDs to build p(z)
                    for (uint_t iz : range(nzfit)) {
                        pz.safe[iz] = 0.0;
                        for (uint_t it : range(ntemplate)) {
                            uint_t im = iz*ntemplate + it;
                            pmodel.safe[im] = exp(-0.5*(pmodel.safe[im] - chi2_best))*prior.safe[im];
                            pz.safe[iz] += pmodel.safe[im];
                        }
                    }

                    // Convolve p(z) with Gaussian to avoid too many peaks
                    double pmax = 0.0;
                    for (uint_t iz : range(nzfit)) {
                        // Convolve
                        uint_t i0 = iz > wconv ? iz - wconv : 0;
                        uint_t i1 = iz < nzfit-1-wconv ? iz + wconv : nzfit-1;
                        pzc.safe[iz] = 0.0;
                        for (uint_t tiz : range(i0, i1+1)) {
                            uint_t iconv = (int_t(tiz)-int_t(iz)) + wconv;
                            pzc.safe[iz] += pz.safe[tiz]*gconv.safe[iconv];
                        }

                        if (pzc.safe[iz] > pmax) {
                            pmax = pzc.safe[iz];
                        }
                    }

                    // Clip too low values, and then apply this back to the p(z,SED)
                    // This is done to mimick as closely as possible the behavior of BPZ.
                    double tprob = 0.0;
                    for (uint_t iz : range(nzfit)) {
                        // Clip
                        if (pzc.safe[iz] < 1e-2*pmax) {
                            pzc.safe[iz] = 0;
                        }

                        // Update the p(z,SED)
                        for (uint_t it : range(ntemplate)) {
                            uint_t im = iz*ntemplate + it;

                            if (pz.safe[iz] > 0.0) {
                                pmodel.safe[im] *= pzc.safe[iz]/pz.safe[iz];
                            }

                            tprob += pmodel.safe[im];
                        }
                    }

                    pmodel /= tprob;

                    if (write_cache) {
                        zmeas_ml.safe[i] = zfit.safe[iml/ntemplate];
                        zmeas_ma.safe[i] = total(zfit*pzc)/total(pzc);

                        cache_bmodel.safe[i] = iml;
                        cache_pmodel += pmodel;
                    }

                    // Maximum likelihood
                    // ------------------

                    metrics tml(bpz_q11_fit.safe[iml], bpz_q12_fit.safe[iml], bpz_q22_fit.safe[iml]);
                    ml += tml;
                    ml2 += sqr(tml);

                    // Marginalization
                    // ---------------

                    metrics tma;
                    for (uint_t im : range(nmodel)) {
                        tma += pmodel.safe[im]*metrics(
                            bpz_q11_fit.safe[im], bpz_q12_fit.safe[im], bpz_q22_fit.safe[im]
                        );
                    }

                    ma += tma;
                    ma2 += sqr(tma);
                }

                ml /= nmc;
                ml2 /= nmc;
                ml2 = sqrt(ml2 - sqr(ml));
                ma /= nmc;
                ma2 /= nmc;
                ma2 = sqrt(ma2 - sqr(ma));

                m_ml.add(id_type, tngal*ml);
                m_ma.add(id_type, tngal*ma);

                if (write_cache) {
                    cache_pmodel /= nmc;

                    zmeas_cache.update_elements("im",        id_mass,      fits::at(iter));
                    zmeas_cache.update_elements("it",        id_type,      fits::at(iter));
                    zmeas_cache.update_elements("idisk",     id_disk,      fits::at(iter));
                    zmeas_cache.update_elements("ibulge",    id_bulge,     fits::at(iter));
                    zmeas_cache.update_elements("ibt",       id_bt,        fits::at(iter));
                    zmeas_cache.update_elements("pmag",      pmag,         fits::at(iter));
                    zmeas_cache.update_elements("pmodel",    cache_pmodel, fits::at(iter,_));
                    zmeas_cache.update_elements("bmodel",    cache_bmodel, fits::at(iter,_));
                    zmeas_cache.update_elements("zmeas_ml",  zmeas_ml,     fits::at(iter,_));
                    zmeas_cache.update_elements("zmeas_ma",  zmeas_ma,     fits::at(iter,_));
                    zmeas_cache.update_elements("ngal",      tngal,        fits::at(iter));
                    zmeas_cache.update_elements("fbulge",    fbulge,       fits::at(iter,_));
                    zmeas_cache.update_elements("fdisk",     fdisk,        fits::at(iter,_));
                    zmeas_cache.update_elements("e1_true",   tr.e1,        fits::at(iter));
                    zmeas_cache.update_elements("e2_true",   tr.e2,        fits::at(iter));
                    zmeas_cache.update_elements("r2_true",   tr.r2,        fits::at(iter));
                    zmeas_cache.update_elements("e1_bpz_ml", ml.e1,        fits::at(iter));
                    zmeas_cache.update_elements("e2_bpz_ml", ml.e2,        fits::at(iter));
                    zmeas_cache.update_elements("r2_bpz_ml", ml.r2,        fits::at(iter));
                    zmeas_cache.update_elements("e1_bpz_ma", ma.e1,        fits::at(iter));
                    zmeas_cache.update_elements("e2_bpz_ma", ma.e2,        fits::at(iter));
                    zmeas_cache.update_elements("r2_bpz_ma", ma.r2,        fits::at(iter));
                    zmeas_cache.update_elements("e1_bpz_ml_err", ml2.e1,   fits::at(iter));
                    zmeas_cache.update_elements("e2_bpz_ml_err", ml2.e2,   fits::at(iter));
                    zmeas_cache.update_elements("r2_bpz_ml_err", ml2.r2,   fits::at(iter));
                    zmeas_cache.update_elements("e1_bpz_ma_err", ma2.e1,   fits::at(iter));
                    zmeas_cache.update_elements("e2_bpz_ma_err", ma2.e2,   fits::at(iter));
                    zmeas_cache.update_elements("r2_bpz_ma_err", ma2.r2,   fits::at(iter));
                    zmeas_cache.open(cache_filename);
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

    vec1d get_madau(double z, const vec1d& lam) const {
        // Madau 1995 extinction for a galaxy spectrum at redshift z
        // Taken straight from BPZ

        vec1d tau(lam.size());

        const vec1d& l = {0.1261, 0.1026, 0.0973, 0.0950};
        const vec1d& c = {3.6e-3, 1.7e-3, 1.2e-3, 9.3e-4};

        if (lam[0] > l[0]*(1.0+z)) {
            return replicate(1.0, lam.size());
        }

        // Lyman series absorption
        double ll = 0.0912;
        uint_t i1 = lower_bound(lam, ll);
        for (uint_t i : range(l)) {
            uint_t i2 = lower_bound(lam, l[i]*(1.0+z));
            tau[i1-_-i2] += c[i]*pow(lam[i1-_-i2]/l[i], 3.46);
        }

        if (lam[0] > ll*(1.0+z)) {
            return exp(-tau);
        }

        // Photoelectric absorption
        uint_t i2 = lower_bound(lam, ll*(1.0+z));
        vec1d x = lam[i1-_-i2]/ll;
        vec1d x3 = pow(x, 3.0);
        tau[i1-_-i2] += 0.25*x3*(pow(1.0+z, 0.46) - pow(x, 0.46))
                      + 9.4*pow(x, 1.5)*(pow(1.0+z, 0.18) - pow(x, 0.18))
                      - 0.7*x3*(pow(x, -1.32) - pow(1.0+z, -1.32))
                      - 0.023*(pow(1.0+z, 1.68) - pow(x, 1.68));

        return exp(-clamp(tau, 0, 700));
    }

    bool average_redshift_bin(uint_t iz) {
        if (!read_egg_psfs(iz)) return false;

        uint_t ntz = max(floor((zb[iz+1]-zb[iz])/dz), 1);
        vec1d uzf = zb[iz] + dz*dindgen(ntz);

        // Start averaging
        vec<1,metrics_set> zml(ntz), zma(ntz), ztr(ntz);
        vec1d dndz(ntz), dndz_qu(ntz), dndz_sf(ntz);
        for (uint_t itz : range(ntz)) {
            print(itz, "/", ntz);

            double zf = uzf.safe[itz];
            std::string zid = replace(to_string(format::fixed(format::precision(zf, 2))), ".", "p");

            // Pre-select redshift slice in EGG PSF library
            {
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
            }

            // Enforce redshift if asked
            if (force_true_z) {
                zfit = {zf};
                nzfit = 1;
                nmodel = ntemplate;

                bpz_q11_fit.resize(ntemplate);
                bpz_q12_fit.resize(ntemplate);
                bpz_q22_fit.resize(ntemplate);

                // Interpolate BPZ PSF library at the redshift
                for (uint_t t : range(ntemplate)) {
                    vec1u ids = where(bpz_ssed == bpz_seds_sid[t]);
                    bpz_q11_fit[t] = interpolate_3spline(bpz_q11[ids], zfit_base, zf);
                    bpz_q12_fit[t] = interpolate_3spline(bpz_q12[ids], zfit_base, zf);
                    bpz_q22_fit[t] = interpolate_3spline(bpz_q22[ids], zfit_base, zf);
                }
            } else {
                zfit = zfit_base;
                nzfit = nzfit_base;
                nmodel = nmodel_base;

                bpz_q11_fit = bpz_q11;
                bpz_q12_fit = bpz_q12;
                bpz_q22_fit = bpz_q22;
            }

            // Reset
            m_ml.reset();
            m_ma.reset();
            m_tr.reset();
            ngal = 0; nqu = 0; nsf = 0;

            // Create workspace
            tpl_flux.resize(nmodel, nband);
            wflux.resize(nband);
            rflux.resize(nband);
            weight.resize(nband);
            wmodel.resize(nmodel, nband);
            wmm.resize(nmodel);
            pmodel.resize(nmodel);
            prior.resize(nmodel);
            tt_prior.resize(nzfit*3);
            pz.resize(nzfit);
            pzc.resize(nzfit);

            if (write_cache) {
                zmeas_ml.resize(nmc);
                zmeas_ma.resize(nmc);
                cache_bmodel.resize(nmc);
                cache_pmodel.resize(nmodel);
            }

            // Pre-compute number of iterations
            note("compute number of iterations...");
            niter = 0;
            just_count = true;
            generate(zf, dz);
            note("done: ", niter);

            // Initialize cache
            std::string cache_id = hash(use_capak_library, apply_igm,
                zfit, bpz_seds, ninterp, gauss_convolve, id_prior, phot_err2, nmc, niter);
            cache_filename = "cache-bpz-z"+zid+"-"+cache_id+".fits";
            note("cache file: ", cache_filename);

            cache_available = false;
            if (use_cache && file::exists(cache_filename)) {
                cache_available = true;
                try {
                    zmeas_cache.open(cache_filename);
                } catch (...) {
                    warning("cache file ", cache_filename, " was corrupted");
                    warning("will remove this cache and create a new one");

                    zmeas_cache.close();
                    cache_available = false;
                }
            }

            if (!cache_available && write_cache) {
                note("creating cache on disk...");

                file::remove(cache_filename);
                zmeas_cache.open(cache_filename);

                // Create cache arrays
                zmeas_cache.allocate_column<uint_t>("im",           niter);
                zmeas_cache.allocate_column<uint_t>("it",           niter);
                zmeas_cache.allocate_column<uint_t>("idisk",        niter);
                zmeas_cache.allocate_column<uint_t>("ibulge",       niter);
                zmeas_cache.allocate_column<uint_t>("ibt",          niter);
                zmeas_cache.allocate_column<float>("pmag",          niter);
                zmeas_cache.allocate_column<uint_t>("bmodel",       niter, nmc);
                zmeas_cache.allocate_column<float>("pmodel",        niter, nmodel);
                zmeas_cache.allocate_column<float>("zmeas_ml",      niter, nmc);
                zmeas_cache.allocate_column<float>("zmeas_ma",      niter, nmc);
                zmeas_cache.allocate_column<float>("fbulge",        niter, nband+1);
                zmeas_cache.allocate_column<float>("fdisk",         niter, nband+1);
                zmeas_cache.allocate_column<float>("ngal",          niter);
                zmeas_cache.allocate_column<float>("e1_true",       niter);
                zmeas_cache.allocate_column<float>("e2_true",       niter);
                zmeas_cache.allocate_column<float>("r2_true",       niter);
                zmeas_cache.allocate_column<float>("e1_bpz_ml",     niter);
                zmeas_cache.allocate_column<float>("e2_bpz_ml",     niter);
                zmeas_cache.allocate_column<float>("r2_bpz_ml",     niter);
                zmeas_cache.allocate_column<float>("e1_bpz_ma",     niter);
                zmeas_cache.allocate_column<float>("e2_bpz_ma",     niter);
                zmeas_cache.allocate_column<float>("r2_bpz_ma",     niter);
                zmeas_cache.allocate_column<float>("e1_bpz_ml_err", niter);
                zmeas_cache.allocate_column<float>("e2_bpz_ml_err", niter);
                zmeas_cache.allocate_column<float>("r2_bpz_ml_err", niter);
                zmeas_cache.allocate_column<float>("e1_bpz_ma_err", niter);
                zmeas_cache.allocate_column<float>("e2_bpz_ma_err", niter);
                zmeas_cache.allocate_column<float>("r2_bpz_ma_err", niter);

                // Save meta data
                zmeas_cache.write_columns("z_grid", zfit, "sed_grid", bpz_seds);
                zmeas_cache.write_columns("m_grid", m, "bt_grid", bt);
                vec1s tbands(nband+1);
                vec1f lambda(nband+1);
                tbands[0] = selection_band;
                lambda[0] = selection_filter.rlam;
                tbands[1-_] = bands;
                for (uint_t l : range(filters)) {
                    lambda[l+1] = filters[l].rlam;
                }
                zmeas_cache.write_columns("bands", tbands, "lambda", lambda);
                zmeas_cache.flush();

                note("done.");
            }

            if (!cache_available) {
                // Pre-compute BPZ template fluxes
                note("pre-compute BPZ template fluxes...");
                for (uint_t t : range(ntemplate)) {
                    vec1d rlam, rsed;
                    ascii::read_table(bpz_seds[t], rlam, rsed);
                    rsed *= 1e-19;

                    for (uint_t izf : range(nzfit)) {
                        // Apply IGM absorption
                        vec1d olam = rlam*(1.0 + zfit[izf]);
                        vec1d osed;
                        if (apply_igm) {
                            osed = get_madau(zfit[izf], olam)*rsed;
                        } else {
                            osed = rsed;
                        }

                        osed = cgs2uJy(olam, osed);
                        olam *= 1e-4;

                        for (uint_t l : range(nband)) {
                            double flx = sed2flux(filters[l].lam, filters[l].res, olam, osed);
                            if (!is_finite(flx)) {
                                // Falling out of the filter, assuming zero flux
                                flx = 0.0;
                            }
                            tpl_flux.safe(izf*ntemplate+t,l) = flx;
                        }
                    }
                }
                note("done.");
            }

            // Compute average at that redshift
            iter = 0;
            just_count = false;
            pgi = progress_start(niter);
            generate(zf, dz);

            // Average quantities
            m_ml.normalize(ngal, nqu, nsf);
            m_ma.normalize(ngal, nqu, nsf);
            m_tr.normalize(ngal, nqu, nsf);

            // Store
            zml[itz] = m_ml;
            zma[itz] = m_ma;
            ztr[itz] = m_tr;

            dndz[itz]    = ngal/dz;
            dndz_qu[itz] = nqu/dz;
            dndz_sf[itz] = nsf/dz;
        }

        // Average over N(z)
        double ntot    = integrate(uzf, dndz);
        double ntot_qu = integrate(uzf, dndz_qu);
        double ntot_sf = integrate(uzf, dndz_sf);
        m_ml = integrate(uzf, zml, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);
        m_ma = integrate(uzf, zma, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);
        m_tr = integrate(uzf, ztr, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);

        // Save
        to_fits("psf-mean-z"+to_string(iz)+"-bpz-ml.fits", m_ml,
            uzf, dndz, dndz_qu, dndz_sf, zml
        );
        to_fits("psf-mean-z"+to_string(iz)+"-bpz-ma.fits", m_ma,
            uzf, dndz, dndz_qu, dndz_sf, zma
        );
        to_fits("psf-mean-z"+to_string(iz)+"-bpz-tr.fits", m_tr,
            uzf, dndz, dndz_qu, dndz_sf, ztr
        );

        return true;
    }
};

int phypp_main(int argc, char* argv[]) {
    // Survey definition
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";

    // Mock photometry
    vec1s filters = {"sdss-u", "sdss-g", "sdss-r", "sdss-i", "sdss-z", "vista-Y", "vista-J", "vista-H"};
    vec1f depths =  {24.5,     24.4,     24.1,     24.1,     23.7,     23.2,     23.2,      23.2}; // 10 sigma
    double min_mag_err = 0.03;
    uint_t nmc = 200;
    uint_t ninterp = 2;
    double dz = 0.01;
    uint_t seds_step = 5;
    double zfit_max = 7.0;
    double zfit_dz = 0.01;
    std::string prior_filter = "sdss-i";
    bool use_capak_library = true;
    bool apply_igm = true;
    bool force_true_z = false;
    bool no_noise = false;
    bool write_cache = true;
    bool use_cache = true;
    uint_t iz = 5;

    read_args(argc, argv, arg_list(
        maglim, selection_band, filters, depths, nmc, min_mag_err, prior_filter, ninterp, dz,
        seds_step, use_capak_library, apply_igm, zfit_max, zfit_dz, write_cache, use_cache, iz,
        force_true_z, no_noise
    ));

    psf_averager pavg;
    pavg.write_cache = write_cache;
    pavg.use_cache = use_cache;

    // Setup survey
    egg::generator_options opts;
    opts.share_dir = "/home/cschreib/code/egg/share/";
    opts.filter_db = "/home/cschreib/work_psf/scripts/psf-averager/filters.dat";
    opts.filter_flambda = true;
    opts.filter_photons = true;
    opts.trim_filters = true;
    opts.selection_band = selection_band;
    opts.filters = filters;
    opts.maglim = maglim;
    opts.logmass_steps = 50;
    opts.bt_steps = 5;
    // opts.mass_max = 12.0;
    opts.seds_step = seds_step;
    pavg.initialize(opts);

    // Setup redshift fitting
    fitter_options fopts;
    fopts.depths = depths;
    fopts.nmc = nmc;
    fopts.min_mag_err = min_mag_err;
    fopts.prior_filter = prior_filter;
    fopts.ninterp = ninterp;
    fopts.dz = dz;
    fopts.zfit_max = zfit_max;
    fopts.zfit_dz = zfit_dz;
    fopts.use_capak_library = use_capak_library;
    fopts.apply_igm = apply_igm;
    fopts.force_true_z = force_true_z;
    fopts.no_noise = no_noise;
    pavg.configure_fitter(fopts);

    // Average PSF metrics
    // for (uint_t iz : range(pavg.zb)) {
    //     if (!pavg.average_redshift_bin(iz)) continue;
    // }

    if (!pavg.average_redshift_bin(iz)) return 1;

    return 0;
}
