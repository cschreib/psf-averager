#ifndef FITTER_BPZ_INCLUDED
#define FITTER_BPZ_INCLUDED

#include "fitter-base.hpp"
#include "common.hpp"

struct fitter_options {
    std::string prior_filter;
    std::string sed_dir;
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
    double min_mag_err = 0.0;
};

class bpz_fitter : public fitter_base {
public :
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

        vec1d wflux, rflux, weight;
        vec2d wmodel;
        vec1d wmm;

        vec2d cache_pmodel;
        vec1u cache_bmodel;
    };

    workspace global;

    // Config
    double gauss_convolve = dnan;
    double rel_err = 0.0;
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
    bool no_prior = false;
    std::string sed_dir;


    explicit bpz_fitter(filter_database& db, psf_moments& pm) : fitter_base(db, pm) {
        code_name = "bpz";
    }

    void do_read_options(program_arguments& opts) override {
        fitter_options fopts;

        // Behavior
        fopts.zfit_max = 7.0;
        fopts.zfit_dz = 0.01;
        fopts.prior_filter = "sdss-i";
        fopts.apply_igm = true;
        fopts.force_true_z = false;
        fopts.min_mag_err = 0.03;
        fopts.gauss_convolve = 0.021; // 0.03/sqrt(2) to match default BPZ

        // Library
        fopts.sed_dir = "/home/cschreib/code/euclid_psf/psf-averager/seds/bpz/";
        fopts.ninterp = 2;
        fopts.use_capak_library = true;
        fopts.use_noline_library = false;
        fopts.use_egg_library = false;

        // Output
        fopts.cache_save_pmodel = false;

        opts.read(arg_list(
            fopts.prior_filter, fopts.zfit_max, fopts.zfit_dz,
            fopts.apply_igm, fopts.force_true_z, fopts.min_mag_err, fopts.sed_dir, fopts.ninterp,
            fopts.use_capak_library, fopts.use_noline_library, fopts.use_egg_library,
            fopts.cache_save_pmodel, fopts.gauss_convolve, no_prior
        ));

        initialize(fopts);
    }

    void initialize(const fitter_options& fopts) {
        nband = filters.size();

        // Choice of template library
        use_capak_library = fopts.use_capak_library;
        use_noline_library = fopts.use_noline_library;
        use_egg_library = fopts.use_egg_library;

        // Options on the fit
        apply_igm = fopts.apply_igm;
        force_true_z = fopts.force_true_z;
        cache_save_pmodel = fopts.cache_save_pmodel;
        sed_dir = file::directorize(fopts.sed_dir);

        // Kernel used to convolve the P(z) to wash out too small fluctuations
        gauss_convolve = fopts.gauss_convolve;

        // Interpolation of BPZ SEDs
        ninterp = fopts.ninterp;

        if (use_egg_library) {
            use_capak_library = true;
            use_noline_library = true;
            ninterp = 0;
        }

        if (!no_prior) {
            vif_check(is_any_of(fopts.prior_filter, bands),
                "prior filter is not in the filter list ('", fopts.prior_filter, "' not found')");
            id_prior = where_first(fopts.prior_filter == bands);
        }

        // Setup redshift grid
        nzfit_base = ceil(fopts.zfit_max/fopts.zfit_dz);
        zfit_base = fopts.zfit_dz*indgen<double>(nzfit_base);

        // Relative error on flux, sets minimum uncertainties
        rel_err = fopts.min_mag_err*(log(10.0)/2.5);

        // Build (and resample) BPZ library
        make_sed_library();

        // Make convolver kernel
        {
            uint_t npt = 18*gauss_convolve/fopts.zfit_dz;
            if (npt % 2 == 0) ++npt;
            wconv = npt/2;

            // This is a gaussian trucated between -9*sigma and +9*sigma
            vec1d kz = fopts.zfit_dz*(indgen<double>(npt)-wconv);
            gconv = integrate_gauss(kz-fopts.zfit_dz/2.0, kz+fopts.zfit_dz/2.0, 0.0, gauss_convolve);
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

        ntemplate = bpz_seds.size();

        vif_check(ntemplate != 0, "could not find any template in '", tsed_dir, "'");

        // Sort SEDs by color (red to blue)
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

        if (!no_prior) {
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

            note("SEDs: ", nelliptical, " ellipticals, ",
                           nspiral, " spirals, ",
                           nirregular, " irregulars");
        } else {
            nelliptical = ntemplate;
            nspiral = 0;
            nirregular = 0;
        }

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

    void set_priors(workspace& w, const vec1d& ftot) {
        double mag = -2.5*log10(ftot.safe[id_prior]) + 23.9;

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

    fit_result do_fit(uint_t iter, const vec1d& ftot) override {
        fit_result fr(nmc);

        workspace local;
        if (multi_threaded) {
            reset_workspace(local);
        }

        workspace& w = (multi_threaded ? local : global);

        // Setup priors
        if (!no_prior) {
            set_priors(w, ftot);
        }

        // Create noise-free photometry
        for (uint_t l : range(nband)) {
            double bftot = ftot.safe[l];
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
            double bchi2 = finf;
            for (uint_t im : range(nmodel)) {
                const double* wm = &w.wmodel.safe(im,0);
                const double* wf = &w.rflux.safe[0];

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

                if (tchi2 < bchi2) {
                    bchi2 = tchi2;
                    iml = im;
                }
            }

            // Create likelihood, apply prior, and marginalize over SEDs to build p(z)
            for (uint_t iz : range(nzfit)) {
                w.pz.safe[iz] = 0.0;
                for (uint_t it : range(ntemplate)) {
                    uint_t im = iz*ntemplate + it;
                    w.pmodel.safe[im] = exp(-0.5*(w.pmodel.safe[im] - bchi2))*w.prior.safe[im];
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

            fr.chi2.safe[i] = bchi2;
            fr.z_obs.safe[i] = zfit.safe[iml/ntemplate];
            fr.z_obsm.safe[i] = total(zfit*w.pzc)/total(w.pzc);

            w.cache_bmodel.safe[i] = iml;

            for (uint_t im : range(nmodel)) {
                w.cache_pmodel.safe(i,im) = w.pmodel.safe[im];
            }

            if (!no_psf) {
                // Maximum likelihood
                fr.psf_obs.safe[i] = metrics(bpz_q11.safe[iml], bpz_q12.safe[iml], bpz_q22.safe[iml]);

                // Marginalization
                metrics& tm = fr.psf_obsm.safe[i];
                for (uint_t im : range(nmodel)) {
                    metrics ttm(bpz_q11.safe[im], bpz_q12.safe[im], bpz_q22.safe[im]);
                    tm += w.cache_pmodel.safe(i,im)*ttm;
                }
            }
        }

        w.cache_pmodel /= nmc;

        // Save things in the cache
        if (write_cache) {
            // Could be executed concurrently, use mutex when in multithreading context
            auto lock = (multi_threaded ?
                std::unique_lock<std::mutex>(*cache_mutex_ptr) : std::unique_lock<std::mutex>());

            if (cache_save_pmodel) {
                cache_ptr->update_elements("pmodel", w.cache_pmodel, fits::at(iter,_,_));
            }
            cache_ptr->update_elements("bmodel", w.cache_bmodel, fits::at(iter,_));
        }

        return fr;
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

    std::string make_cache_hash() override {
        return hash(use_capak_library, use_egg_library, use_noline_library,
            apply_igm, zfit, bpz_seds, ninterp, gauss_convolve, cache_save_pmodel,
            bands, phot_err2, nmc, rel_err);
    }

    void reset_workspace(workspace& w) {
        w.pmodel.resize(nmodel);
        w.prior.resize(nmodel);
        w.tt_prior.resize(nzfit*3);
        w.pz.resize(nzfit);
        w.pzc.resize(nzfit);

        w.wflux.resize(nband);
        w.rflux.resize(nband);
        w.weight.resize(nband);
        w.wmodel.resize(nmodel, nband);
        w.wmm.resize(nmodel);

        w.cache_bmodel.resize(nmc);
        w.cache_pmodel.resize(nmc, nmodel);

        if (no_prior) {
            w.prior[_] = 1.0;
        }
    }

    void do_prepare_fit(double zf) override {
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
        if (!no_psf && (force_true_z || bpz_q11.empty())) {
            bpz_q11.resize(nmodel);
            bpz_q12.resize(nmodel);
            bpz_q22.resize(nmodel);
            compute_moments = true;
        }

        bool compute_fluxes = false;
        if (tpl_flux.empty()) {
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
                        psf.get_moments(
                            olam, osed,
                            bpz_q11.safe[iz*ntemplate+it], bpz_q12.safe[iz*ntemplate+it],
                            bpz_q22.safe[iz*ntemplate+it]
                        );
                    }
                }
            }
        }

        // Create workspace
        if (!multi_threaded) {
            reset_workspace(global);
        }
    }

    void do_initialize_cache() override {
        // Initialize arrays
        cache_ptr->allocate_column<uint_t>("bmodel", niter, nmc);
        if (cache_save_pmodel) {
            cache_ptr->allocate_column<float>("pmodel",  niter, nmodel);
        }

        // Save meta data
        cache_ptr->write_columns("z_grid", zfit, "sed_grid", bpz_seds);
    }

    void save_individuals(const std::string& filename) override {
        // Write to disk the individual measurements
        fits::table otbl(filename);
        otbl.update_column("z_grid",   zfit);
        otbl.update_column("sed_grid", bpz_seds);
    }
};

#endif
