#include "psf-averager-photoz.hpp"

struct fitter_options {
    std::string prior_filter;
    double gauss_convolve = 0.021; // 0.03/sqrt(2) to match default BPZ
    uint_t ninterp = 2;
    double zfit_max = 7.0;
    double zfit_dz = 0.01;
    bool use_capak_library = true;
    bool use_egg_library = false;
    bool use_noline_library = false;
    bool apply_igm = true;
    bool force_true_z = false;
    bool cache_save_pmodel = true;
    std::string sed_dir;
};

class bpz_averager : public psf_averager {
public :
    // Averages
    metrics_set m_ml; // maximum likelihood
    metrics_set m_ma; // marginalization

    vec<1,metrics_set> zml, zma;

    // Individuals
    vec<2,metrics> indiv_ml, indiv_ma;
    vec2f indiv_chi2, indiv_zml, indiv_zma;

    // BPZ PSF library
    vec1d bpz_q11, bpz_q12, bpz_q22;

    // Fit models
    vec1s bpz_seds;
    vec1d zfit_base;
    uint_t nzfit_base = npos;
    uint_t nmodel_base = npos;
    vec1d zfit;
    uint_t nmodel = npos;
    uint_t nzfit = npos;
    uint_t ntemplate = npos;
    uint_t nelliptical = npos;
    uint_t nspiral = npos;
    uint_t nirregular = npos;

    // Internal variables (constant per slice)
    vec2d tpl_flux;

    // Internal variables (workspace)
    struct workspace {
        vec1d pmodel;
        vec1d prior;
        vec1d tt_prior;
        vec1d pz, pzc;
        vec1d chi2_best;

        vec1d wflux, rflux, weight;
        vec2d wmodel;
        vec1d wmm;

        vec1d zmeas_ml, zmeas_ma;
        vec2d cache_pmodel;
        vec1u cache_bmodel;
    };

    workspace global;

    // Config
    double gauss_convolve = dnan;
    vec1d gconv;
    uint_t wconv = npos;
    uint_t id_prior = npos;
    uint_t ninterp = 0;
    bool use_capak_library = true;
    bool use_egg_library = false;
    bool use_noline_library = false;
    bool apply_igm = true;
    bool force_true_z = false;
    bool cache_save_pmodel = true;
    std::string sed_dir;

    bpz_averager() : psf_averager("bpz") {}

    void configure_fitter(fitter_options opts) {
        nband = filters.size();

        // Choice of template library
        use_capak_library = opts.use_capak_library;
        use_noline_library = opts.use_noline_library;
        use_egg_library = opts.use_egg_library;

        // Options on the fit
        apply_igm = opts.apply_igm;
        force_true_z = opts.force_true_z;
        cache_save_pmodel = opts.cache_save_pmodel;
        sed_dir = opts.sed_dir;

        // Kernel used to convolve the P(z) to wash out too small fluctuations
        gauss_convolve = opts.gauss_convolve;

        // Interpolation of BPZ SEDs
        ninterp = opts.ninterp;

        if (use_egg_library) {
            use_capak_library = true;
            use_noline_library = true;
            ninterp = 0;
        }

        phypp_check(is_any_of(opts.prior_filter, bands),
            "prior filter is not in the filter list ('", opts.prior_filter, "' not found')");
        id_prior = where_first(opts.prior_filter == bands) + 1;

        // Setup redshift grid
        nzfit_base = ceil(opts.zfit_max/opts.zfit_dz);
        zfit_base = opts.zfit_dz*dindgen(nzfit_base);

        // Build (and resample) BPZ library
        make_sed_library();

        // Make convolver kernel
        {
            uint_t npt = 18*gauss_convolve/opts.zfit_dz;
            if (npt % 2 == 0) ++npt;
            wconv = npt/2;

            // This is a gaussian trucated between -9*sigma and +9*sigma
            vec1d kz = opts.zfit_dz*(findgen(npt)-wconv);
            gconv = integrate_gauss(kz-opts.zfit_dz/2.0, kz+opts.zfit_dz/2.0, 0.0, gauss_convolve);
            gconv /= total(gconv);
        }
    }

    void make_sed_library() {
        // Locate SEDs
        std::string tsed_dir;
        if (use_egg_library) {
            tsed_dir = "SED_egg/";
        } else {
            if (use_capak_library) {
                if (use_noline_library) {
                    tsed_dir = "SED_capak_noline/";
                } else {
                    tsed_dir = "SED_capak/";
                }
            } else {
                tsed_dir = "SED/";
            }
        }

        tsed_dir = sed_dir+tsed_dir;

        bpz_seds = file::list_files(tsed_dir, "*.sed");
        if (bpz_seds.empty()) {
            bpz_seds = file::list_files(tsed_dir, "*.dat");
        }

        // Sort BPZ SEDs by color (red to blue)
        vec1d color(bpz_seds.size()); {
            for (uint_t t : range(bpz_seds)) {
                vec1d rlam, rsed;
                ascii::read_table(tsed_dir+bpz_seds[t], rlam, rsed);
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

        // Split SEDs in classes
        if (use_egg_library) {
            // EGG SEDs, we do our best to find something that matches...
            nelliptical = count(color <= 0.02);
            nspiral = count(color > 0.02 && color <= 0.15);
        } else {
            // Standard set
            // Capak: color = {0.0023, 0.036, 0.097, 0.21, 0.22, 0.41}
            //                 Ell        Spirals       Irregulars
            nelliptical = 1;
            nspiral = 2;
        }

        nirregular = ntemplate - nelliptical - nspiral;

        note("SEDs: ", nelliptical, " ellipticals, ", nspiral, " spirals, ", nirregular, " irregulars");

        nmodel_base = ntemplate*nzfit_base;

        for (uint_t i : range(bpz_seds)) {
            note(" - ", i, ": ", bpz_seds[i]);
        }

        bpz_seds = tsed_dir+bpz_seds;

        // Interpolate the templates
        if (ninterp > 0) {
            note("interpolating library...");

            std::string tmp_dir = "BPZ_interp/";
            file::mkdir(tmp_dir);

            // Clear library
            vec1s oseds = bpz_seds;
            bpz_seds.clear();

            // Add first SED
            bpz_seds.push_back(oseds[0]);

            for (uint_t it : range(ntemplate-1)) {
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

                for (uint_t ii : range(ninterp)) {
                    // Interpolate
                    double x = (ii+1.0)/(ninterp+1.0);
                    vec1d tsed1 = (1.0-x)*rsed1;
                    vec1d tsed2 = x*rsed2;
                    vec1d rsed = tsed1 + tsed2;
                    rsed = uJy2cgs(clam, rsed)*1e19;

                    // Save to disk
                    std::string fname = file::remove_extension(file::get_basename(oseds[it]))
                        +"_"+to_string(ii)+"_"+to_string(ninterp)+".sed";
                    ascii::write_table(tmp_dir+fname, clam, rsed);

                    // Insert SED in library
                    bpz_seds.push_back(tmp_dir+fname);
                }

                // Add next SED
                bpz_seds.push_back(oseds[it+1]);
            }

            uint_t ntemplate_old = ntemplate;
            ntemplate = bpz_seds.size();

            nmodel_base = nzfit_base*ntemplate;

            note("expanded BPZ template library from ", ntemplate_old, " to ", ntemplate, " SEDs");
        }
    }

    void set_priors(workspace& w, const vec1d& fdisk, const vec1d& fbulge) {
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
                    w.tt_prior.safe[im] = val;
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
                    w.tt_prior.safe[im] = pow(zfit.safe[iz], alpha.safe[tt])
                        *exp(-clamp(pow(zfit.safe[iz]/zm, alpha.safe[tt]), 0.0, 700.0));

                    tprior += w.tt_prior.safe[im];
                }

                double tprior2 = 0.0;
                for (uint_t iz : range(nzfit)) {
                    uint_t im = iz*3 + tt;

                    // Normalize
                    w.tt_prior.safe[im] /= tprior;

                    if (!use_capak_library) {
                        // Clip low probability wings
                        if (w.tt_prior.safe[im] < 1e-2/nzfit) {
                            w.tt_prior.safe[im] = 0.0;
                        } else {
                            tprior2 += w.tt_prior.safe[im];
                        }
                    }
                }

                if (!use_capak_library) {
                    for (uint_t iz : range(nzfit)) {
                        uint_t im = iz*3 + tt;

                        // Normalize again
                        w.tt_prior.safe[im] /= tprior2;
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
                w.prior.safe[im] = w.tt_prior.safe[iz*3 + tt0];
            } else {
                // SEDs interpolated between different classes get interpolated priors
                double x = (itm+1.0)/(ninterp+1.0);
                w.prior.safe[im] = w.tt_prior.safe[iz*3 + tt0]*(1.0-x) + x*w.tt_prior.safe[iz*3 + tt1];
            }
        }
    }

    void compute_averages(uint_t iter, const workspace& w, uint_t id_type, double tngal) {
        // Maximum likelihood
        // ------------------

        metrics ml;
        metrics ml2;
        for (uint_t i : range(nmc)) {
            uint_t iml = w.cache_bmodel.safe(i);
            metrics tm(bpz_q11.safe[iml], bpz_q12.safe[iml], bpz_q22.safe[iml]);
            ml  += tm;
            ml2 += sqr(tm);

            if (keep_individuals_in_memory) {
                indiv_ml.safe(iter,i) = tm;
            }
        }

        ml /= nmc;
        ml2 /= nmc;
        ml2 = sqrt(ml2 - sqr(ml));


        // Marginalization
        // ---------------

        metrics ma;
        metrics ma2;
        for (uint_t i : range(nmc)) {
            metrics tm;
            for (uint_t im : range(nmodel)) {
                metrics ttm(bpz_q11.safe[im], bpz_q12.safe[im], bpz_q22.safe[im]);
                tm += w.cache_pmodel.safe(i,im)*ttm;
            }

            ma  += tm;
            ma2 += sqr(tm);

            if (keep_individuals_in_memory) {
                indiv_ma.safe(iter,i) = tm;
            }
        }

        ma /= nmc;
        ma2 /= nmc;
        ma2 = sqrt(ma2 - sqr(ma));

        // Store average values if asked
        if (keep_averages_in_memory) {
            // Could be executed concurrently, use mutex when in multithreading context
            auto lock = (nthread > 0 ?
                std::unique_lock<std::mutex>(avg_mutex) : std::unique_lock<std::mutex>());

            m_ml.add(id_type, tngal*ml);
            m_ma.add(id_type, tngal*ma);
        }

        if (write_cache) {
            fitter_cache.update_elements("e1_bpz_ml",     ml.e1,  fits::at(iter));
            fitter_cache.update_elements("e2_bpz_ml",     ml.e2,  fits::at(iter));
            fitter_cache.update_elements("r2_bpz_ml",     ml.r2,  fits::at(iter));
            fitter_cache.update_elements("e1_bpz_ma",     ma.e1,  fits::at(iter));
            fitter_cache.update_elements("e2_bpz_ma",     ma.e2,  fits::at(iter));
            fitter_cache.update_elements("r2_bpz_ma",     ma.r2,  fits::at(iter));
            fitter_cache.update_elements("e1_bpz_ml_err", ml2.e1, fits::at(iter));
            fitter_cache.update_elements("e2_bpz_ml_err", ml2.e2, fits::at(iter));
            fitter_cache.update_elements("r2_bpz_ml_err", ml2.r2, fits::at(iter));
            fitter_cache.update_elements("e1_bpz_ma_err", ma2.e1, fits::at(iter));
            fitter_cache.update_elements("e2_bpz_ma_err", ma2.e2, fits::at(iter));
            fitter_cache.update_elements("r2_bpz_ma_err", ma2.r2, fits::at(iter));
        }
    }

    void process_cached(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        workspace local;
        workspace& w = (nthread == 0 ? global : local);

        fitter_cache.read_elements("bmodel", w.cache_bmodel, fits::at(iter,_));
        fitter_cache.read_elements("pmodel", w.cache_pmodel, fits::at(iter,_,_));

        compute_averages(iter, w, id_type, tngal);
    }

    void do_fit(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        workspace local;
        if (nthread > 0) {
            reset_workspace(local);
        }

        workspace& w = (nthread == 0 ? global : local);

        // Setup priors
        set_priors(w, fdisk, fbulge);

        // Create noise-free photometry
        for (uint_t l : range(nband)) {
            double bftot = fdisk.safe[l+1] + fbulge.safe[l+1];
            w.weight.safe[l] = 1.0/sqrt(phot_err2.safe[l] + sqr(bftot*rel_err));
            w.rflux.safe[l] = w.wflux.safe[l] = bftot*w.weight.safe[l];
        }

        // Create models
        for (uint_t im : range(nmodel)) {
            w.wmm.safe[im] = 0.0;
            for (uint_t l : range(nband)) {
                w.wmodel.safe(im,l) = tpl_flux.safe(im,l)*w.weight.safe[l];
                w.wmm.safe[im] += sqr(w.wmodel.safe(im,l));
            }
        }

        // Simulate redshift measurements
        for (uint_t i : range(nmc)) {
            for (uint_t im : range(nmodel)) {
                w.cache_pmodel.safe(i,im) = 0.0;
            }

            // Create noisy photometry
            if (!no_noise) {
                for (uint_t l : range(nband)) {
                    // In weighted units, the random perturbations have a sigma of unity
                    w.rflux.safe[l] = w.wflux.safe[l] + mc.safe(i,l);
                }
            }

            // Find best SED and z, and store chi2 to build p(z,SED)
            uint_t iml = npos;
            w.chi2_best.safe[i] = finf;
            for (uint_t im : range(nmodel)) {
                double* wm = &w.wmodel.safe(im,0);
                double* wf = &w.rflux.safe[0];

                double wfm = 0.0;
                for (uint_t l : range(nband)) {
                    wfm += wf[l]*wm[l];
                }

                double scale = wfm/w.wmm.safe[im];

                double tchi2 = 0.0;
                for (uint_t l : range(nband)) {
                    tchi2 += sqr(wf[l] - scale*wm[l]);
                }

                w.pmodel.safe[im] = tchi2;

                if (tchi2 < w.chi2_best.safe[i]) {
                    w.chi2_best.safe[i] = tchi2;
                    iml = im;
                }
            }

            // Create likelihood, apply prior, and marginalize over SEDs to build p(z)
            for (uint_t iz : range(nzfit)) {
                w.pz.safe[iz] = 0.0;
                for (uint_t it : range(ntemplate)) {
                    uint_t im = iz*ntemplate + it;
                    w.pmodel.safe[im] = exp(-0.5*(w.pmodel.safe[im] - w.chi2_best.safe[i]))*w.prior.safe[im];
                    w.pz.safe[iz] += w.pmodel.safe[im];
                }
            }

            // Convolve p(z) with Gaussian to avoid too many peaks
            double pmax = 0.0;
            for (uint_t iz : range(nzfit)) {
                // Convolve
                uint_t i0 = iz > wconv ? iz - wconv : 0;
                uint_t i1 = iz < nzfit-1-wconv ? iz + wconv : nzfit-1;
                w.pzc.safe[iz] = 0.0;
                for (uint_t tiz : range(i0, i1+1)) {
                    uint_t iconv = (int_t(tiz)-int_t(iz)) + wconv;
                    w.pzc.safe[iz] += w.pz.safe[tiz]*gconv.safe[iconv];
                }

                if (w.pzc.safe[iz] > pmax) {
                    pmax = w.pzc.safe[iz];
                }
            }

            // Clip too low values, and then apply this back to the p(z,SED)
            // This is done to mimick as closely as possible the behavior of BPZ.
            double tprob = 0.0;
            for (uint_t iz : range(nzfit)) {
                // Clip
                if (w.pzc.safe[iz] < 1e-2*pmax) {
                    w.pzc.safe[iz] = 0;
                }

                // Update the p(z,SED)
                for (uint_t it : range(ntemplate)) {
                    uint_t im = iz*ntemplate + it;

                    if (w.pz.safe[iz] > 0.0) {
                        w.pmodel.safe[im] *= w.pzc.safe[iz]/w.pz.safe[iz];
                    }

                    tprob += w.pmodel.safe[im];
                }
            }

            w.pmodel /= tprob;

            double zma = total(zfit*w.pzc)/total(w.pzc);

            if (keep_individuals_in_memory) {
                indiv_chi2.safe(iter,i) = w.chi2_best.safe[i];
                indiv_zml.safe(iter,i) = zfit.safe[iml/ntemplate];
                indiv_zma.safe(iter,i) = zma;
            }

            if (write_cache) {
                w.zmeas_ml.safe[i] = zfit.safe[iml/ntemplate];
                w.zmeas_ma.safe[i] = zma;
            }

            w.cache_bmodel.safe[i] = iml;

            for (uint_t im : range(nmodel)) {
                w.cache_pmodel.safe(i,im) = w.pmodel.safe[im];
            }
        }

        w.cache_pmodel /= nmc;

        if (write_cache) {
            if (cache_save_pmodel) {
                fitter_cache.update_elements("pmodel", w.cache_pmodel, fits::at(iter,_,_));
            }
            fitter_cache.update_elements("bmodel",    w.cache_bmodel, fits::at(iter,_));
            fitter_cache.update_elements("best_chi2", w.chi2_best,    fits::at(iter,_));
            fitter_cache.update_elements("zmeas_ml",  w.zmeas_ml,     fits::at(iter,_));
            fitter_cache.update_elements("zmeas_ma",  w.zmeas_ma,     fits::at(iter,_));
        }

        compute_averages(iter, w, id_type, tngal);
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
        if (keep_averages_in_memory) {
            zml.resize(uzf.size());
            zma.resize(uzf.size());
        }
    }

    std::string make_cache_hash() override {
        return hash(use_capak_library, use_egg_library, use_noline_library,
            apply_igm, zfit, bpz_seds, ninterp, gauss_convolve, cache_save_pmodel);
    }

    void reset_workspace(workspace& w) {
        w.pmodel.resize(nmodel);
        w.prior.resize(nmodel);
        w.tt_prior.resize(nzfit*3);
        w.pz.resize(nzfit);
        w.pzc.resize(nzfit);
        w.chi2_best.resize(nmc);

        w.wflux.resize(nband);
        w.rflux.resize(nband);
        w.weight.resize(nband);
        w.wmodel.resize(nmodel, nband);
        w.wmm.resize(nmodel);

        w.cache_bmodel.resize(nmc);
        w.cache_pmodel.resize(nmc, nmodel);

        if (write_cache) {
            w.zmeas_ml.resize(nmc);
            w.zmeas_ma.resize(nmc);
        }
    }

    void initialize_redshift_slice(uint_t itz) override {
        double zf = uzf[itz];

        // Enforce redshift if asked
        if (force_true_z) {
            zfit = {zf};
            nzfit = 1;
            nmodel = ntemplate;
        } else {
            zfit = zfit_base;
            nzfit = nzfit_base;
            nmodel = nmodel_base;
        }

        // Pre-compute BPZ template fluxes and PSF moments
        bool compute_moments = false;
        if (force_true_z || bpz_q11.empty()) {
            bpz_q11.resize(nmodel);
            bpz_q12.resize(nmodel);
            bpz_q22.resize(nmodel);
            compute_moments = true;
        }

        bool compute_fluxes = false;
        if (!cache_available || tpl_flux.empty()) {
            tpl_flux.resize(nmodel, nband);
            compute_fluxes = true;
        }

        if (compute_moments || compute_fluxes) {
            for (uint_t it : range(ntemplate)) {
                vec1d rlam, rsed;
                ascii::read_table(bpz_seds[it], rlam, rsed);
                rsed *= 1e-19;

                for (uint_t iz : range(nzfit)) {
                    // Apply IGM absorption
                    vec1d olam = rlam*(1.0 + zfit[iz]);
                    vec1d osed = rsed;
                    if (apply_igm) {
                        osed *= get_madau(zfit[iz], olam);
                    }

                    osed = cgs2uJy(olam, osed);
                    olam *= 1e-4;

                    if (compute_fluxes) {
                        // Compute model fluxes
                        for (uint_t l : range(nband)) {
                            double flx = sed2flux(filters[l].lam, filters[l].res, olam, osed);
                            if (!is_finite(flx)) {
                                // Falling out of the filter, assuming zero flux
                                flx = 0.0;
                            }
                            tpl_flux.safe(iz*ntemplate+it,l) = flx;
                        }
                    }

                    if (compute_moments) {
                        // Compute PSF moments
                        double fvis = sed2flux(psf_filter.lam, psf_filter.res, olam, osed);
                        bpz_q11.safe[iz*ntemplate+it] = sed2flux(
                            psf_filter.lam, psf_filter.res*mono_q11, olam, osed
                        )/fvis;
                        bpz_q12.safe[iz*ntemplate+it] = sed2flux(
                            psf_filter.lam, psf_filter.res*mono_q12, olam, osed
                        )/fvis;
                        bpz_q22.safe[iz*ntemplate+it] = sed2flux(
                            psf_filter.lam, psf_filter.res*mono_q22, olam, osed
                        )/fvis;
                    }
                }
            }
        }

        // Reset averages
        if (keep_averages_in_memory) {
            m_ml.reset();
            m_ma.reset();
        }

        // Create workspace
        if (nthread == 0) {
            reset_workspace(global);
        }

        // Initialize individual arrays
        if (keep_individuals_in_memory) {
            indiv_ml.resize(niter,nmc);
            indiv_ma.resize(niter,nmc);
            indiv_chi2.resize(niter,nmc);
            indiv_zml.resize(niter,nmc);
            indiv_zma.resize(niter,nmc);
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
        if (keep_averages_in_memory) {
            // Average quantities
            m_ml.normalize(ngal, nqu, nsf);
            m_ma.normalize(ngal, nqu, nsf);

            // Store
            zml[itz] = m_ml;
            zma[itz] = m_ma;
        }

        if (write_individuals) {
            // Write to disk the individual measurements
            fits::update_table(indiv_filename,
                "e1_obs",   get_e1(indiv_ml),
                "e2_obs",   get_e2(indiv_ml),
                "r2_obs",   get_r2(indiv_ml),
                "e1_obsm",  get_e1(indiv_ma),
                "e2_obsm",  get_e2(indiv_ma),
                "r2_obsm",  get_r2(indiv_ma),
                "chi2_obs", indiv_chi2,
                "z_obs",    indiv_zml,
                "z_obsm",   indiv_zma,
                "z_grid",   zfit,
                "sed_grid", bpz_seds
            );
        }
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
    // External data
    std::string share_dir = "/home/cschreib/code/egg-analytic/share/";
    std::string filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
    std::string sed_lib = "/home/cschreib/code/egg-analytic/share/opt_lib_fastpp_hd_noigm.fits";
    std::string sed_imf = "chabrier";
    std::string psf_file  = "/home/cschreib/code/euclid_psf/psf-averager/psf-mono.fits";
    std::string sed_dir   = "/home/cschreib/code/euclid_psf/psf-averager/seds/bpz/";

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
    bool use_egg_library = false;
    bool apply_igm = true;
    bool force_true_z = false;
    bool no_noise = false;
    bool write_cache = false;
    bool use_cache = false;
    bool cache_save_pmodel = true;
    bool write_individuals = false;
    bool write_averages = true;
    uint_t nthread = 0;
    uint_t iz = 5;
    std::string cache_id;

    read_args(argc, argv, arg_list(
        maglim, selection_band, filters, depths, nmc, min_mag_err, prior_filter, ninterp, dz,
        seds_step, use_capak_library, use_noline_library, apply_igm, zfit_max, zfit_dz, write_cache,
        use_cache, iz, force_true_z, no_noise, use_egg_library, cache_save_pmodel, share_dir,
        filter_db, psf_file, sed_dir, nthread, write_individuals, write_averages, cache_id, sed_lib,
        sed_imf
    ));

    bpz_averager pavg;
    pavg.write_cache = write_cache;
    pavg.use_cache = use_cache;

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
    opts.filters = filters;
    opts.maglim = maglim;
    opts.logmass_steps = 50;
    opts.bt_steps = 5;
    opts.logmass_max = 12.0;
    opts.seds_step = seds_step;
    opts.nthread = nthread;
    pavg.initialize(opts);

    // Setup mock
    mock_options mopts;
    mopts.depths = depths;
    mopts.nmc = nmc;
    mopts.min_mag_err = min_mag_err;
    mopts.dz = dz;
    mopts.no_noise = no_noise;
    mopts.psf_file = psf_file;
    mopts.write_individuals = write_individuals;
    mopts.write_averages = write_averages;
    mopts.force_cache_id = cache_id;
    pavg.configure_mock(mopts);

    // Setup redshift fitting
    fitter_options fopts;
    fopts.prior_filter = prior_filter;
    fopts.ninterp = ninterp;
    fopts.zfit_max = zfit_max;
    fopts.zfit_dz = zfit_dz;
    fopts.use_capak_library = use_capak_library;
    fopts.use_noline_library = use_noline_library;
    fopts.use_egg_library = use_egg_library;
    fopts.apply_igm = apply_igm;
    fopts.force_true_z = force_true_z;
    fopts.cache_save_pmodel = cache_save_pmodel;
    fopts.sed_dir = sed_dir;
    pavg.configure_fitter(fopts);

    // Average PSF metrics
    // for (uint_t iz : range(pavg.zb)) {
    //     if (!pavg.average_redshift_bin(iz)) continue;
    // }

    if (!pavg.average_redshift_bin(iz)) return 1;

    return 0;
}
