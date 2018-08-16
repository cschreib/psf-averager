#include "psf-averager-photoz.hpp"

struct fitter_options {
    std::string prior_filter;
    double gauss_convolve = 0.021; // 0.03/sqrt(2) to match default BPZ
    uint_t ninterp = 2;
    double zfit_max = 7.0;
    double zfit_dz = 0.01;
    bool use_capak_library = true;
    bool use_noline_library = false;
    bool apply_igm = true;
    bool force_true_z = false;
    bool cache_save_pmodel = true;
};

class bpz_averager : public psf_averager {
public :
    // Averages
    metrics_set m_ml; // maximum likelihood
    metrics_set m_ma; // marginalization

    vec<1,metrics_set> zml, zma;

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
    uint_t ntemplate = npos;
    uint_t nelliptical = npos;
    uint_t nspiral = npos;
    uint_t nirregular = npos;

    // Internal variables
    vec1d pmodel;
    vec1d prior;
    vec1d tt_prior;
    vec1d pz, pzc;
    vec1d chi2_best;

    vec1d wflux, rflux, weight;
    vec2d wmodel;
    vec1d wmm;
    vec2d tpl_flux;

    vec1d zmeas_ml, zmeas_ma, cache_pmodel;
    vec1u cache_bmodel;

    // Config
    double gauss_convolve = dnan;
    double dzfit = dnan;
    vec1d gconv;
    uint_t wconv = npos;
    uint_t id_prior = npos;
    uint_t ninterp = 0;
    bool use_capak_library = true;
    bool use_noline_library = false;
    bool apply_igm = true;
    bool force_true_z = false;
    bool cache_save_pmodel = true;


    bpz_averager() : psf_averager("bpz") {}

    void configure_fitter(const fitter_options& opts) {
        nband = filters.size();
        use_capak_library = opts.use_capak_library;
        use_noline_library = opts.use_noline_library;
        apply_igm = opts.apply_igm;
        force_true_z = opts.force_true_z;

        phypp_check(is_any_of(opts.prior_filter, bands),
            "prior filter is not in the filter list ('", opts.prior_filter, "' not found')");
        id_prior = where_first(opts.prior_filter == bands) + 1;
        cache_save_pmodel = opts.cache_save_pmodel;

        // Kernel used to convolve the P(z) to wash out too small fluctuations
        gauss_convolve = opts.gauss_convolve;

        // Interpolation of BPZ SEDs
        ninterp = opts.ninterp;

        // Read BPZ PSF library
        std::string bpz_filename;
        if (use_capak_library) {
            if (use_noline_library) {
                bpz_filename = psf_dir+"BPZ_capak_noline/psfs-rebin2-cst.txt";
            } else {
                bpz_filename = psf_dir+"BPZ_capak/psfs-rebin2-cst.txt";
            }
        } else {
            bpz_filename = psf_dir+"BPZ/psfs-rebin2-cst.txt";
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
        } else {
            note("read BPZ PSF library, ", bpz_ssed.size(), " elements (",
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
            if (use_noline_library) {
                sed_dir = "/home/cschreib/programming/bpz-1.99.3/SED_capak_noline/";
            } else {
                sed_dir = "/home/cschreib/programming/bpz-1.99.3/SED_capak/";
            }
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
            vec1u sid = sort(bpz_seds);
            vec1u ids(nmodel_base);
            for (uint_t izf : range(nzfit_base)) {
                uint_t i0 = izf*ntemplate;
                for (uint_t it : range(ntemplate)) {
                    ids.safe[i0+sid.safe[it]] = i0+it;
                }
            }

            // Reshuffle
            // bpz_ssed is left alone, because we want it to remain sorted and map to
            // SED IDs in our new SED ordering, not the one of the PSF library...
            // bpz_zz is also left alone because the reshufling does not affect it.
            bpz_q11 = bpz_q11[ids];
            bpz_q12 = bpz_q12[ids];
            bpz_q22 = bpz_q22[ids];
        }

        bpz_seds = sed_dir+bpz_seds;
        bpz_seds_sid = to_string_vector(indgen(ntemplate));

        // Interpolate the templates
        if (ninterp > 0) {
            note("interpolating library...");

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
            append(bpz_q11,  oq11.safe[id0]);
            append(bpz_q12,  oq12.safe[id0]);
            append(bpz_q22,  oq22.safe[id0]);

            uint_t ntemplate_old = ntemplate;
            ntemplate = 1;

            // Cache VIS fluxes
            vec2d vis(ntemplate_old, nzfit_base);
            for (uint_t it : range(ntemplate_old)) {
                vec1d rlam, rsed;
                ascii::read_table(oseds[it], rlam, rsed);
                rsed = cgs2uJy(rlam, rsed*1e-19);

                for (uint_t iz : range(nzfit_base)) {
                    vec1d tlam = rlam*(1e-4*(1.0 + zfit_base[iz]));
                    double flx = sed2flux(psf_filter.lam, psf_filter.res, tlam, rsed);
                    if (!is_finite(flx)) {
                        // Falling out of the filter, assuming zero flux
                        flx = 0.0;
                    }

                    vis.safe(it,iz) = flx;
                }
            }

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
                        double f1 = (1.0-x)*vis.safe(it,iz);
                        double f2 = x*vis.safe(it+1,iz);
                        double xf = f2/(f1 + f2);

                        bpz_q11.push_back(oq11.safe[id0.safe[iz]]*(1.0-xf) + xf*oq11.safe[id1.safe[iz]]);
                        bpz_q12.push_back(oq12.safe[id0.safe[iz]]*(1.0-xf) + xf*oq12.safe[id1.safe[iz]]);
                        bpz_q22.push_back(oq22.safe[id0.safe[iz]]*(1.0-xf) + xf*oq22.safe[id1.safe[iz]]);
                    }
                }

                // Add last SED
                bpz_seds.push_back(oseds[it+1]);
                bpz_seds_sid.push_back(to_string(it+1));
                append(bpz_ssed, replicate(bpz_seds_sid.back(), nzfit_base));
                append(bpz_zz,   ubz);
                append(bpz_q11,  oq11.safe[id1]);
                append(bpz_q12,  oq12.safe[id1]);
                append(bpz_q22,  oq22.safe[id1]);
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

    void set_priors(const vec1d& fdisk, const vec1d& fbulge) override {
        double mag = -2.5*log10(fdisk.safe[id_prior] + fbulge.safe[id_prior]) + 23.9;

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

    void compute_averages(uint_t id_type, double tngal) {
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
            ma  += cache_pmodel.safe[im]*tm;
            ma2 += cache_pmodel.safe[im]*sqr(tm);
        }
        ma2 = sqrt(ma2 - sqr(ma));

        m_ma.add(id_type, tngal*ma);

        if (write_cache) {
            fitter_cache.update_elements("e1_bpz_ml", ml.e1,      fits::at(iter));
            fitter_cache.update_elements("e2_bpz_ml", ml.e2,      fits::at(iter));
            fitter_cache.update_elements("r2_bpz_ml", ml.r2,      fits::at(iter));
            fitter_cache.update_elements("e1_bpz_ma", ma.e1,      fits::at(iter));
            fitter_cache.update_elements("e2_bpz_ma", ma.e2,      fits::at(iter));
            fitter_cache.update_elements("r2_bpz_ma", ma.r2,      fits::at(iter));
            fitter_cache.update_elements("e1_bpz_ml_err", ml2.e1, fits::at(iter));
            fitter_cache.update_elements("e2_bpz_ml_err", ml2.e2, fits::at(iter));
            fitter_cache.update_elements("r2_bpz_ml_err", ml2.r2, fits::at(iter));
            fitter_cache.update_elements("e1_bpz_ma_err", ma2.e1, fits::at(iter));
            fitter_cache.update_elements("e2_bpz_ma_err", ma2.e2, fits::at(iter));
            fitter_cache.update_elements("r2_bpz_ma_err", ma2.r2, fits::at(iter));
        }
    }

    void process_cached(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        fitter_cache.read_elements("bmodel", cache_bmodel, fits::at(iter,_));
        fitter_cache.read_elements("pmodel", cache_pmodel, fits::at(iter,_));

        compute_averages(id_type, tngal);
    }

    void do_fit(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

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

        for (uint_t im : range(nmodel)) {
            cache_pmodel.safe[im] = 0.0;
        }

        // Simulate redshift measurements
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
            chi2_best.safe[i] = finf;
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

                if (tchi2 < chi2_best.safe[i]) {
                    chi2_best.safe[i] = tchi2;
                    iml = im;
                }
            }

            // Create likelihood, apply prior, and marginalize over SEDs to build p(z)
            for (uint_t iz : range(nzfit)) {
                pz.safe[iz] = 0.0;
                for (uint_t it : range(ntemplate)) {
                    uint_t im = iz*ntemplate + it;
                    pmodel.safe[im] = exp(-0.5*(pmodel.safe[im] - chi2_best.safe[i]))*prior.safe[im];
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
            }

            cache_bmodel.safe[i] = iml;
            cache_pmodel += pmodel;
        }

        cache_pmodel /= nmc;

        if (write_cache) {
            if (cache_save_pmodel) {
                fitter_cache.update_elements("pmodel", cache_pmodel, fits::at(iter,_));
            }
            fitter_cache.update_elements("bmodel",    cache_bmodel, fits::at(iter,_));
            fitter_cache.update_elements("best_chi2", chi2_best,    fits::at(iter,_));
            fitter_cache.update_elements("zmeas_ml",  zmeas_ml,     fits::at(iter,_));
            fitter_cache.update_elements("zmeas_ma",  zmeas_ma,     fits::at(iter,_));
        }

        compute_averages(id_type, tngal);
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

    void initialize_redshift_bin(uint_t iz) override {
        zml.resize(uzf.size());
        zma.resize(uzf.size());
    }

    std::string make_cache_hash() override {
        return hash(use_capak_library, use_egg_library, use_noline_library,
            apply_igm, zfit, bpz_seds, ninterp, gauss_convolve, cache_save_pmodel);
    }

    void initialize_redshift_slice(uint_t itz) override {
        double zf = uzf[itz];

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

        // Reset averages
        m_ml.reset();
        m_ma.reset();

        // Create workspace
        pmodel.resize(nmodel);
        prior.resize(nmodel);
        tt_prior.resize(nzfit*3);
        pz.resize(nzfit);
        pzc.resize(nzfit);
        tpl_flux.resize(nmodel, nband);
        chi2_best.resize(nmc);

        wflux.resize(nband);
        rflux.resize(nband);
        weight.resize(nband);
        wmodel.resize(nmodel, nband);
        wmm.resize(nmodel);

        cache_bmodel.resize(nmc);
        cache_pmodel.resize(nmodel);

        if (write_cache) {
            zmeas_ml.resize(nmc);
            zmeas_ma.resize(nmc);
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
    }

    void initialize_cache() override {
        // Initialize arrays
        fitter_cache.allocate_column<uint_t>("bmodel",       niter, nmc);
        if (cache_save_pmodel) {
            fitter_cache.allocate_column<float>("pmodel",    niter, nmodel);
        }
        fitter_cache.allocate_column<float>("best_chi2",     niter, nmc);
        fitter_cache.allocate_column<float>("zmeas_ml",      niter, nmc);
        fitter_cache.allocate_column<float>("zmeas_ma",      niter, nmc);
        fitter_cache.allocate_column<float>("e1_bpz_ml",     niter);
        fitter_cache.allocate_column<float>("e2_bpz_ml",     niter);
        fitter_cache.allocate_column<float>("r2_bpz_ml",     niter);
        fitter_cache.allocate_column<float>("e1_bpz_ma",     niter);
        fitter_cache.allocate_column<float>("e2_bpz_ma",     niter);
        fitter_cache.allocate_column<float>("r2_bpz_ma",     niter);
        fitter_cache.allocate_column<float>("e1_bpz_ml_err", niter);
        fitter_cache.allocate_column<float>("e2_bpz_ml_err", niter);
        fitter_cache.allocate_column<float>("r2_bpz_ml_err", niter);
        fitter_cache.allocate_column<float>("e1_bpz_ma_err", niter);
        fitter_cache.allocate_column<float>("e2_bpz_ma_err", niter);
        fitter_cache.allocate_column<float>("r2_bpz_ma_err", niter);

        // Save meta data
        fitter_cache.write_columns("z_grid", zfit, "sed_grid", bpz_seds);
    }

    void finalize_redshift_slice(uint_t itz) override {
        // Average quantities
        m_ml.normalize(ngal, nqu, nsf);
        m_ma.normalize(ngal, nqu, nsf);

        // Store
        zml[itz] = m_ml;
        zma[itz] = m_ma;
    }

    void finalize_redshift_bin(uint_t iz, double ntot, double ntot_qu, double ntot_sf) override {
        // Average over N(z)
        m_ml = integrate(uzf, zml, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);
        m_ma = integrate(uzf, zma, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);

        // Save
        to_fits("psf-mean-z"+to_string(iz)+"-bpz-ml.fits", m_ml,
            uzf, dndz, dndz_qu, dndz_sf, zml
        );
        to_fits("psf-mean-z"+to_string(iz)+"-bpz-ma.fits", m_ma,
            uzf, dndz, dndz_qu, dndz_sf, zma
        );
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
    bool use_noline_library = false;
    bool apply_igm = true;
    bool force_true_z = false;
    bool no_noise = false;
    bool write_cache = true;
    bool use_cache = true;
    bool cache_save_pmodel = true;
    uint_t iz = 5;

    read_args(argc, argv, arg_list(
        maglim, selection_band, filters, depths, nmc, min_mag_err, prior_filter, ninterp, dz,
        seds_step, use_capak_library, use_noline_library, apply_igm, zfit_max, zfit_dz, write_cache,
        use_cache, iz, force_true_z, no_noise, cache_save_pmodel
    ));

    bpz_averager pavg;
    pavg.write_cache = write_cache;
    pavg.use_cache = use_cache;

    // Setup survey
    egg::generator_options opts;
    opts.share_dir = "/home/cschreib/code/egg/share/";
    opts.filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
    opts.filter_flambda = true;
    opts.filter_photons = true;
    opts.trim_filters = true;
    opts.selection_band = selection_band;
    opts.filters = filters;
    opts.maglim = maglim;
    opts.logmass_steps = 50;
    opts.bt_steps = 5;
    opts.logmass_max = 12.0;
    opts.seds_step = seds_step;
    pavg.initialize(opts);

    // Setup mock
    mock_options mopts;
    mopts.depths = depths;
    mopts.nmc = nmc;
    mopts.min_mag_err = min_mag_err;
    mopts.dz = dz;
    mopts.no_noise = no_noise;
    mopts.psf_dir = "../../psf-library/";
    pavg.configure_mock(mopts);

    // Setup redshift fitting
    fitter_options fopts;
    fopts.prior_filter = prior_filter;
    fopts.ninterp = ninterp;
    fopts.zfit_max = zfit_max;
    fopts.zfit_dz = zfit_dz;
    fopts.use_capak_library = use_capak_library;
    fopts.use_noline_library = use_noline_library;
    fopts.apply_igm = apply_igm;
    fopts.force_true_z = force_true_z;
    fopts.cache_save_pmodel = cache_save_pmodel;
    pavg.configure_fitter(fopts);

    // Average PSF metrics
    // for (uint_t iz : range(pavg.zb)) {
    //     if (!pavg.average_redshift_bin(iz)) continue;
    // }

    if (!pavg.average_redshift_bin(iz)) return 1;

    return 0;
}
