#include "psf-averager-photoz.hpp"

struct fitter_options {
    std::string prior_file;
    std::string prior_filter;
    double zfit_max = 7.0;
    double zfit_dz = 0.01;
    std::string template_error;
    double template_error_amp = 0.5;
    double fit_tftol = 1e-4;
    bool apply_igm = true;
    bool use_noline_library = false;
    bool use_egg_library = false;
    uint_t egg_sed_step = 1;
    bool force_true_z = false;
    uint_t limited_set = 0;
    bool cache_save_pmodel = true;
};

class eazy_averager : public psf_averager {
public :
    // Averages
    metrics_set m_ml; // maximum likelihood
    metrics_set m_ma; // marginalization

    vec<1,metrics_set> zml, zma;

    // EAzY PSF library
    vec1d eazy_q11, eazy_q12, eazy_q22;
    vec1d eazy_fvis;

    // Fit models
    vec1s eazy_seds;
    vec1d zfit_base;
    uint_t nzfit_base = npos;
    uint_t nmodel_base = npos;
    vec1d zfit;
    uint_t nmodel = npos;
    uint_t nzfit = npos;
    uint_t ntemplate = npos;

    // Internal variables
    vec2d prior_table_base;
    vec2d prior_table;
    vec1d prior_mag;
    vec1d prior;

    vec1d tpl_error_x, tpl_error_y;

    vec2d pz;
    vec1d tpz;
    vec1d chi2_best;

    vec1d wflux, rflux, weight;
    vec2d wmodel;
    vec2d tpl_flux;
    vec2d tpl_err;

    vec1u cache_bestz;
    vec2d cache_bestc;
    vec1d cache_pz;
    vec1d cache_pc;

    // Config
    double dzfit = dnan;
    std::string prior_file, prior_filter;
    uint_t id_prior = npos;
    bool apply_igm = true;
    bool force_true_z = false;
    bool use_noline_library = false;
    bool use_egg_library = false;
    bool cache_save_pmodel = true;
    uint_t egg_sed_step = 1;
    uint_t limited_set = 0;
    double fit_tftol = 1e-4;


    eazy_averager() : psf_averager("eazy") {}

    void configure_fitter(const fitter_options& opts) {
        nband = filters.size();
        apply_igm = opts.apply_igm;
        force_true_z = opts.force_true_z;
        use_noline_library = opts.use_noline_library;
        use_egg_library = opts.use_egg_library;
        egg_sed_step = opts.egg_sed_step;
        limited_set = opts.limited_set;
        fit_tftol = opts.fit_tftol;
        cache_save_pmodel = opts.cache_save_pmodel;

        if (limited_set > 3) {
            warning("'limited_set' can be at most equal to 3, set to zero to use all set");
            limited_set = 3;
        }

        prior_filter = opts.prior_filter;
        prior_file = opts.prior_file;
        phypp_check(is_any_of(prior_filter, bands),
            "prior filter is not in the filter list ('", prior_filter, "' not found')");
        id_prior = where_first(prior_filter == bands) + 1;

        // Read template error file
        ascii::read_table(opts.template_error, tpl_error_x, tpl_error_y);
        tpl_error_x /= 1e4;
        tpl_error_y *= opts.template_error_amp;

        // Setup redshift grid
        nzfit_base = ceil(opts.zfit_max/opts.zfit_dz);
        zfit_base = opts.zfit_dz*dindgen(nzfit_base);

        // Setup EAzY library
        make_sed_library();

        // Read prior
        read_prior();
    }

    void make_sed_library() {
        // List SEDs
        std::string sed_dir = "/home/cschreib/programming/eazy-photoz/templates/";
        if (use_egg_library) {
            eazy_seds = sed_dir+"EGG/"+file::list_files(sed_dir+"EGG/", "*.dat");
            inplace_sort(eazy_seds);

            if (egg_sed_step > 1) {
                // Remove some SEDs to save time
                vec1b keep(eazy_seds.size());
                vec1s lib_sid(eazy_seds.size());
                for (uint_t i : range(eazy_seds)) {
                    lib_sid[i] = eazy_seds[i].substr((sed_dir+"EGG/egg-").size(), 5);
                }

                for (uint_t iuv : range(use.dims[0]))
                for (uint_t ivj : range(use.dims[1])) {
                    std::string suv = align_right(to_string(iuv), 2, '0');
                    std::string svj = align_right(to_string(ivj), 2, '0');
                    std::string sid = suv+"-"+svj;
                    uint_t i = where_first(lib_sid == sid);
                    if (i != npos) {
                        keep[i] = (iuv+ivj) % egg_sed_step == 0;
                    }
                }

                eazy_seds = eazy_seds[where(keep)];
            }
        } else {
            if (use_noline_library) {
                eazy_seds = sed_dir + vec1s{
                    "EAZY_v1.1_noline/eazy_v1.1_sed1.dat",
                    "EAZY_v1.1_noline/eazy_v1.1_sed2.dat",
                    "EAZY_v1.1_noline/eazy_v1.1_sed3.dat",
                    "EAZY_v1.1_noline/eazy_v1.1_sed4.dat",
                    "EAZY_v1.1_noline/eazy_v1.1_sed5.dat",
                    "EAZY_v1.1_noline/eazy_v1.1_sed6.dat",
                    "EAZY_v1.1_noline/eazy_v1.1_sed7.dat",
                    "Dusty/c09_del_8.6_z_0.019_chab_age09.40_av2.0.dat",
                    "erb2010_highEW.dat"
                };
            } else {
                eazy_seds = sed_dir + vec1s{
                    "EAZY_v1.1_lines/eazy_v1.1_sed1.dat",
                    "EAZY_v1.1_lines/eazy_v1.1_sed2.dat",
                    "EAZY_v1.1_lines/eazy_v1.1_sed3.dat",
                    "EAZY_v1.1_lines/eazy_v1.1_sed4.dat",
                    "EAZY_v1.1_lines/eazy_v1.1_sed5.dat",
                    "EAZY_v1.1_lines/eazy_v1.1_sed6.dat",
                    "EAZY_v1.1_lines/eazy_v1.1_sed7.dat",
                    "Dusty/c09_del_8.6_z_0.019_chab_age09.40_av2.0.dat",
                    "erb2010_highEW.dat"
                };
            }
        }

        ntemplate = eazy_seds.size();

        nmodel_base = ntemplate*nzfit_base;

        for (uint_t i : range(eazy_seds)) {
            note(" - ", i, ": ", eazy_seds[i]);
        }
    }

    void read_prior() {
        // Read and re-sample prior
        {
            std::ifstream in(prior_file);
            std::string line;
            std::getline(in, line);
            vec1s spl = split_any_of(line, " ");
            phypp_check(spl.size() > 2 && spl[0] == "#" && spl[1] == "z",
                "ill-formed header line for prior file");

            spl = spl[2-_];
            phypp_check(count(!from_string(spl, prior_mag)) == 0, "could not read prior magnitudes "
                "(", spl, ")");
        }

        vec1d tmp_z;
        vec2d tmp_prior;
        ascii::read_table(prior_file, tmp_z, ascii::columns(prior_mag.size(), tmp_prior));

        prior_table_base.resize(prior_mag.size(), zfit_base.size());
        for (uint_t im : range(prior_mag)) {
            prior_table_base.safe(im,_) = interpolate(tmp_prior.safe(_,im), tmp_z, zfit_base);
        }

        // Erase extrapolated prior values
        // Note: as in EAzY. Not sure this is the right thing to do though.
        prior_table_base(_,where(zfit_base < min(tmp_z) || zfit_base > max(tmp_z))) = 1.0;
    }

    void set_priors(const vec1d& fdisk, const vec1d& fbulge) override {
        double mag = -2.5*log10(fdisk.safe[id_prior] + fbulge.safe[id_prior]) + 23.9;

        uint_t imag = lower_bound(prior_mag, mag);
        if (imag == npos) {
            for (uint_t iz : range(nzfit)) {
                prior.safe[iz] = prior_table.safe(0,iz);
            }
        } else if (imag == prior_mag.size()-1) {
            for (uint_t iz : range(nzfit)) {
                prior.safe[iz] = prior_table.safe(imag,iz);
            }
        } else {
            double x = (mag - prior_mag.safe[imag])/(prior_mag.safe[imag+1] - prior_mag.safe[imag]);
            for (uint_t iz : range(nzfit)) {
                prior.safe[iz] = prior_table.safe(imag,iz)*(1-x) + x*prior_table.safe(imag+1,iz);
            }
        }
    }

    void compute_averages(uint_t id_type, double tngal) {
        // Maximum likelihood
        // ------------------

        metrics ml;
        metrics ml2;
        for (uint_t i : range(nmc)) {
            uint_t iz = cache_bestz.safe[i];

            // Compute flux-weighted average moments for this model
            double q11 = 0.0, q12 = 0.0, q22 = 0.0;
            double wtot = 0.0;
            for (uint_t it : range(ntemplate)) {
                uint_t im = iz*ntemplate + it;

                double w = cache_bestc.safe(i,it)*eazy_fvis.safe[im];
                wtot += w;

                q11 += w*eazy_q11.safe[im];
                q12 += w*eazy_q12.safe[im];
                q22 += w*eazy_q22.safe[im];
            }

            metrics tm(q11/wtot, q12/wtot, q22/wtot);
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
        for (uint_t iz : range(nzfit)) {
            // Compute flux-weighted average moments for this model
            double q11 = 0.0, q12 = 0.0, q22 = 0.0;
            double wtot = 0.0;
            for (uint_t it : range(ntemplate)) {
                uint_t im = iz*ntemplate + it;

                double w = cache_pc.safe(im)*eazy_fvis.safe[im];
                wtot += w;

                q11 += w*eazy_q11.safe[im];
                q12 += w*eazy_q12.safe[im];
                q22 += w*eazy_q22.safe[im];
            }

            metrics tm(q11, q12, q22);
            ma  += cache_pz.safe[iz]*tm;
            ma2 += cache_pz.safe[iz]*sqr(tm);
        }
        ma2 = sqrt(ma2 - sqr(ma));

        m_ma.add(id_type, tngal*ma);

        if (write_cache) {
            fitter_cache.update_elements("e1_eazy_ml", ml.e1,      fits::at(iter));
            fitter_cache.update_elements("e2_eazy_ml", ml.e2,      fits::at(iter));
            fitter_cache.update_elements("r2_eazy_ml", ml.r2,      fits::at(iter));
            fitter_cache.update_elements("e1_eazy_ma", ma.e1,      fits::at(iter));
            fitter_cache.update_elements("e2_eazy_ma", ma.e2,      fits::at(iter));
            fitter_cache.update_elements("r2_eazy_ma", ma.r2,      fits::at(iter));
            fitter_cache.update_elements("e1_eazy_ml_err", ml2.e1, fits::at(iter));
            fitter_cache.update_elements("e2_eazy_ml_err", ml2.e2, fits::at(iter));
            fitter_cache.update_elements("r2_eazy_ml_err", ml2.r2, fits::at(iter));
            fitter_cache.update_elements("e1_eazy_ma_err", ma2.e1, fits::at(iter));
            fitter_cache.update_elements("e2_eazy_ma_err", ma2.e2, fits::at(iter));
            fitter_cache.update_elements("r2_eazy_ma_err", ma2.r2, fits::at(iter));
        }
    }

    void process_cached(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        fitter_cache.read_elements("best_coef", cache_bestc, fits::at(iter,_,_));
        fitter_cache.read_elements("best_z",    cache_bestz, fits::at(iter,_));
        fitter_cache.read_elements("pz",        cache_pz,    fits::at(iter,_));
        fitter_cache.read_elements("pc",        cache_pc,    fits::at(iter,_));

        compute_averages(id_type, tngal);
    }

    void do_fit(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        matrix::mat2d alpha(ntemplate, ntemplate);
        vec1d beta(ntemplate);
        vec1d tcoefs(ntemplate);

        // Simulate redshift measurements
        // We do it in two passes to save computations:
        // 1) Iterate over redshifts, build model matrices, and fit all the MC realizations
        // 2) Iterate over MC simulations and find p(z)

        // Step 1): compute chi2
        for (uint_t i : range(nmc)) {
            chi2_best.safe[i] = finf;
        }

        for (uint_t iz : range(nzfit)) {
            cache_pz.safe[iz] = 0;

            // Create noise-free photometry
            for (uint_t l : range(nband)) {
                double bftot = fdisk.safe[l+1] + fbulge.safe[l+1];
                weight.safe[l] = 1.0/sqrt(phot_err2.safe[l] + sqr(bftot*(rel_err + tpl_err.safe(iz,l))));
                rflux.safe[l] = wflux.safe[l] = bftot*weight.safe[l];
            }

            // Create models
            for (uint_t it : range(ntemplate)) {
                uint_t im = iz*ntemplate + it;
                for (uint_t l : range(nband)) {
                    wmodel.safe(it,l) = tpl_flux.safe(im,l)*weight.safe[l];
                }
            }

            double* wm = &wmodel.safe(0, 0);
            double* wf = &rflux.safe[0];

            // Create alpha matrix
            for (uint_t it0 : range(ntemplate))
            for (uint_t it1 : range(ntemplate)) {
                if (it1 >= it0) {
                    double tmp = 0.0;
                    for (uint_t l : range(nband)) {
                        tmp += wm[it0*nband+l]*wm[it1*nband+l];
                    }

                    alpha.safe(it0,it1) = tmp;
                } else {
                    // Symmetrize
                    alpha.safe(it0,it1) = alpha.safe(it1,it0);
                }
            }

            for (uint_t i : range(nmc)) {
                // Create noisy photometry
                for (uint_t it : range(ntemplate)) {
                    beta.safe[it] = 0.0;
                }

                for (uint_t l : range(nband)) {
                    // In weighted units, the random perturbations have a sigma of unity
                    if (!no_noise) {
                        rflux.safe[l] = wflux.safe[l] + mc.safe(i,l);
                    }

                    // Create beta vector
                    for (uint_t it : range(ntemplate)) {
                        // Ignore negative fluxes (NB: as in EAzY)
                        if (rflux.safe[l] > 0.0) {
                            beta.safe[it] += rflux.safe[l]*wm[it*nband+l];
                        }
                    }
                }

                double tchi2;

                if (limited_set > 0) {
                    tchi2 = finf;

                    // Find best combination using 1, 2, or 3 templates from the library.
                    auto check_chi2 = [&](const vec1d& coefs) {
                        double chi2 = 0.0;
                        for (uint_t l : range(nband)) {
                            double tm = 0.0;
                            for (uint_t it : range(ntemplate)) {
                                tm += coefs.safe[it]*wm[it*nband+l];
                            }
                            chi2 += sqr(wf[l] - tm);
                        }

                        if (chi2 < tchi2) {
                            tchi2 = chi2;
                            tcoefs = coefs;
                        }
                    };

                    auto try_fit_single = [&](uint_t j) {
                        double coef = beta.safe[j]/alpha.safe(j,j);

                        if (coef >= 0) {
                            vec1d coefs(ntemplate);
                            coefs.safe[j] = coef;
                            check_chi2(coefs);
                        }
                    };

                    auto try_fit_two = [&](uint_t j1, uint_t j2) {
                        double det = alpha.safe(j1,j1)*alpha.safe(j2,j2) - sqr(alpha.safe(j1,j2));

                        double coef1 =  beta.safe[j1]*alpha.safe(j2,j2) - beta.safe[j2]*alpha.safe(j1,j2);
                        double coef2 = -beta.safe[j1]*alpha.safe(j1,j2) + beta.safe[j2]*alpha.safe(j1,j1);

                        coef1 /= det;
                        coef2 /= det;

                        if (coef1 >= 0 && coef2 >= 0) {
                            vec1d coefs(ntemplate);
                            coefs.safe[j1] = coef1;
                            coefs.safe[j2] = coef2;
                            check_chi2(coefs);
                        }
                    };

                    auto try_fit_three = [&](uint_t j1, uint_t j2, uint_t j3) {
                        double det1 = alpha.safe(j3,j3)*alpha.safe(j2,j2) - sqr(alpha.safe(j2,j3));
                        double det2 = alpha.safe(j3,j3)*alpha.safe(j1,j2) - alpha.safe(j2,j3)*alpha.safe(j1,j3);
                        double det3 = alpha.safe(j2,j3)*alpha.safe(j1,j2) - alpha.safe(j2,j2)*alpha.safe(j1,j3);
                        double det  = alpha.safe(j1,j1)*det1 - alpha.safe(j1,j2)*det2 + alpha.safe(j1,j3)*det3;

                        double det4 = alpha.safe(j3,j3)*alpha.safe(j1,j1) - sqr(alpha.safe(j1,j3));
                        double det5 = alpha.safe(j2,j3)*alpha.safe(j1,j1) - alpha.safe(j1,j2)*alpha.safe(j1,j3);
                        double det6 = alpha.safe(j2,j2)*alpha.safe(j1,j1) - sqr(alpha.safe(j1,j2));

                        double coef1 =  beta.safe[j1]*det1 - beta.safe[j2]*det2 + beta.safe[j3]*det3;
                        double coef2 = -beta.safe[j1]*det2 + beta.safe[j2]*det4 - beta.safe[j3]*det5;
                        double coef3 =  beta.safe[j1]*det3 - beta.safe[j2]*det5 + beta.safe[j3]*det6;

                        coef1 /= det;
                        coef2 /= det;
                        coef3 /= det;

                        if (coef1 >= 0 && coef2 >= 0 && coef3 >= 0) {
                            vec1d coefs(ntemplate);
                            coefs.safe[j1] = coef1;
                            coefs.safe[j2] = coef2;
                            coefs.safe[j3] = coef3;
                            check_chi2(coefs);
                        }
                    };

                    for (uint_t i1 : range(ntemplate)) {
                        try_fit_single(i1);

                        if (limited_set > 1) {
                            for (uint_t i2 : range(i1+1, ntemplate)) {
                                try_fit_two(i1, i2);

                                if (limited_set > 2) {
                                    for (uint_t i3 : range(i2+1, ntemplate)) {
                                        try_fit_three(i1, i2, i3);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // Find best combination using all templates in the library.
                    // Compute non-negative decomposition using the Sha, Saul & Lee 2006 algorithm.
                    // It is iterative; may not converge and could be non-optimal, but it scales
                    // better than a brute force approach when using large number of templates.

                    // Initialize coefficients (NB: as in EAzY)
                    for (uint_t it : range(ntemplate)) {
                        tcoefs.safe[it] = (beta.safe[it] > 0.0 ? 1.0 : 0.0);
                    }

                    uint_t titer = 0;
                    const uint_t titermax = 10000;

                    double ta, tb;
                    do {
                        ta = 0.0; tb = 0.0;
                        for (uint_t it0 : range(ntemplate)) {
                            double av = 0.0;
                            for (uint_t it1 : range(ntemplate)) {
                                av += alpha.safe(it0,it1)*tcoefs.safe[it1];
                            }

                            // Update coeff
                            double old = tcoefs.safe[it0];
                            tcoefs.safe[it0] *= beta.safe[it0]/av;

                            ta += abs(tcoefs.safe[it0] - old);
                            tb += old;
                        }

                        ++titer;
                    } while (ta/tb > fit_tftol && titer < titermax);

                    // Compute chi2
                    tchi2 = 0.0;
                    for (uint_t l : range(nband)) {
                        double tm = 0.0;
                        for (uint_t it : range(ntemplate)) {
                            tm += tcoefs.safe[it]*wm[it*nband+l];
                        }

                        tchi2 += sqr(wf[l] - tm);
                    }
                }

                // Store results
                pz.safe(i,iz) = tchi2;
                for (uint_t it : range(ntemplate)) {
                    uint_t im = iz*ntemplate + it;
                    cache_pc.safe[im] += tcoefs.safe[it];
                }

                if (tchi2 < chi2_best.safe[i]) {
                    chi2_best.safe[i] = tchi2;
                    cache_bestz.safe[i] = iz;
                    for (uint_t it : range(ntemplate)) {
                        cache_bestc.safe(i,it) = tcoefs.safe[it];
                    }
                }
            }
        }

        // Step 2): compute p(z)
        for (uint_t i : range(nmc)) {
            // Create likelihood, apply prior, and marginalize over SEDs to build p(z)
            double tprob = 0.0;
            for (uint_t iz : range(nzfit)) {
                tpz.safe[iz] = exp(-0.5*(pz.safe(i,iz) - chi2_best.safe[i]))*prior.safe[iz];
                tprob += tpz.safe[iz];
            }

            // Stack p(z)
            for (uint_t iz : range(nzfit)) {
                cache_pz.safe[iz] += tpz.safe[iz]/tprob;
            }
        }

        cache_pz /= nmc;
        cache_pc /= nmc;

        if (write_cache) {
            fitter_cache.update_elements("best_chi2", chi2_best,   fits::at(iter,_));
            fitter_cache.update_elements("best_z",    cache_bestz, fits::at(iter,_));
            if (cache_save_pmodel) {
                fitter_cache.update_elements("best_coef", cache_bestc, fits::at(iter,_,_));
                fitter_cache.update_elements("pz",        cache_pz,    fits::at(iter,_));
                fitter_cache.update_elements("pc",        cache_pc,    fits::at(iter,_));
            }
        }

        compute_averages(id_type, tngal);
    }

    vec1d get_inoue(double z, const vec1d& lam) const {
        // C++ translation of the fortran IGM code from Inoue et al. (2014).
        // Adapted from the C version of EAzY.
        // http://adsabs.harvard.edu/abs/2014arXiv1402.0677I
        // http://www.las.osaka-sandai.ac.jp/~inoue/ANAIGM/ANAIGM.tar.gz

        vec1d tau = replicate(1.0, lam.size());

        // Note: unused in EAzY, don't know why?
        // // Lyman series, LAF
        // auto get_lslaf = [&](double lobs) {
        //     const double z1 = 1.2;
        //     const double z2 = 4.7;
        //     double v = 0.0;
        //     for (uint_t i : range(inoue_nlam)) {
        //         double li = inoue_lam.safe[i];
        //         if (lobs < li*(1.0+z) && lobs > li) {
        //             if (lobs < li*(1.0+z1)) {
        //                 v += inoue_laf1.safe[i]*pow(lobs/li, 1.2);
        //             } else if (lobs < li*(1.0+z2)) {
        //                 v += inoue_laf2.safe[i]*pow(lobs/li, 3.7);
        //             } else {
        //                 v += inoue_laf3.safe[i]*pow(lobs/li, 5.5);
        //             }
        //         }
        //     }

        //     return v;
        // };

        // // Lyman series, DLA
        // auto get_lsdla = [&](double lobs) {
        //     const double z1 = 2.0;
        //     double v = 0.0;
        //     for (uint_t i : range(inoue_nlam)) {
        //         double li = inoue_lam.safe[i];
        //         if (lobs < li*(1.0+z) && lobs > li) {
        //             if (lobs < li*(1.0+z1)) {
        //                 v += inoue_dla1.safe[i]*pow(lobs/li, 2.0);
        //             } else {
        //                 v += inoue_dla2.safe[i]*pow(lobs/li, 3.0);
        //             }
        //         }
        //     }

        //     return v;
        // };

        // Lyman continuum, LAF
        auto get_lclaf = [&](double lobs) {
            const double z1 = 1.2;
            const double z2 = 4.7;
            const double l0 = 911.8;

            if (lobs > l0*(1.0+z)) {
                return 0.0;
            } else if (z < z1) {
                return 0.3248*(pow(lobs/l0, 1.2) - pow(1.0+z, -0.9)*pow(lobs/l0, 2.1));
            } else if (z < z2) {
                if (lobs > l0*(1.0+z1)) {
                    return 2.545e-2*(pow(1.0+z, 1.6)*pow(lobs/l0, 2.1) - pow(lobs/l0, 3.7));
                } else {
                    return 2.545e-2*pow(1.0+z, 1.6)*pow(lobs/l0, 2.1)
                        + 0.3248*pow(lobs/l0, 1.2) - 0.2496*pow(lobs/l0, 2.1);
                }
            } else {
                if (lobs > l0*(1.0+z2)) {
                    return 5.221e-4*(pow(1.0+z, 3.4)*pow(lobs/l0, 2.1) - pow(lobs/l0, 5.5));
                } else if (lobs > l0*(1.0+z1)) {
                    return 5.221e-4*pow(1.0+z, 3.4)*pow(lobs/l0, 2.1)
                        + 0.2182*pow(lobs/l0, 2.1) - 2.545e-2*pow(lobs/l0, 3.7);
                } else {
                    return 5.221e-4*pow(1.0+z, 3.4)*pow(lobs/l0, 2.1)
                        + 0.3248*pow(lobs/l0, 1.2) - 3.140e-2*pow(lobs/l0, 2.1);
                }
            }
        };

        // Lyman continuum, DLA
        auto get_lcdla = [&](double lobs) {
            const double z1 = 2.0;
            const double l0 = 911.8;

            if (lobs > l0*(1.0+z)) {
                return 0.0;
            } else if (z < z1) {
                return 0.2113*pow(1.0+z, 2.0) - 0.07661*pow(1.0+z, 2.3)*pow(lobs/l0, -0.3)
                    - 0.1347*pow(lobs/l0, 2.0);
            } else {
                double tmp = 0.04696*pow(1.0+z, 3.0) - 0.01779*pow(1.0+z, 3.3)*pow(lobs/l0, -0.3);
                if (lobs > l0*(1.0+z1)) {
                    return tmp - 0.02916*pow(lobs/l0, 3.0);
                } else {
                    return tmp - 0.1347*pow(lobs/l0, 2.0) - 0.2905*pow(lobs/l0, -0.3) + 0.6340;
                }
            }
        };

        double da; {
            // Note: same as in FAST
            double l0 = 1026.0*(1.0 + z);
            double l1 = 1216.0*(1.0 + z);
            uint_t nstep = 100;
            vec1d tl = rgen(l0, l1, nstep);
            vec1d ptau = exp(-3.6e-3*pow(tl/1216.0, 3.46));
            da = mean(ptau);
        }

        double db; {
            // Note: more ellaborate than in FAST
            const vec1d ll = {1216.0, 1026.0, 972.8, 950.0, 938.1, 931.0, 926.5, 923.4, 921.2, 919.6,
                918.4, 917.5, 916.7, 916.1, 915.6, 915.2};
            const vec1d aa = {3.6e-03, 1.7e-03, 1.2e-03, 9.4e-04, 8.2e-04, 7.5e-04, 7.1e-04, 6.8e-04,
                6.6e-04, 6.4e-04, 6.3e-04, 6.2e-04, 6.1e-04, 6.0e-04, 6.0e-04, 6.0e-04};

            double l0 = 912.0*(1.0 + z);
            double l1 = 1026.0*(1.0 + z);
            uint_t nstep = 100;
            vec1d tl = rgen(l0, l1, nstep);
            vec1d ptau(nstep);
            for (uint_t i : range(ll)) {
                vec1u idl = where(tl < ll.safe[i]*(1.0+z));
                ptau.safe[idl] += aa.safe[i]*pow(tl.safe[idl]/ll.safe[i], 3.46);
            }

            ptau = exp(-ptau);
            db = mean(ptau);
        }

        uint_t l0 = lower_bound(lam, 912.0);
        uint_t l1 = lower_bound(lam, 1026.0);
        uint_t l2 = lower_bound(lam, 1216.0);

        for (uint_t l : range(l0)) {
            tau.safe[l] = exp(-(get_lclaf(lam.safe[l]*(1.0+z)) + get_lcdla(lam.safe[l]*(1.0+z))));
        }
        for (uint_t l : range(l0, l1)) {
            tau.safe[l] = db;
        }
        for (uint_t l : range(l1, l2)) {
            tau.safe[l] = da;
        }

        return tau;
    }

    void initialize_redshift_bin(uint_t iz) override {
        zml.resize(uzf.size());
        zma.resize(uzf.size());
    }

    std::string make_cache_hash() override {
        return hash(use_noline_library, use_egg_library, egg_sed_step, limited_set,
            prior_file, prior_filter, apply_igm, zfit, eazy_seds, tpl_error_y, fit_tftol);
    }

    void initialize_redshift_slice(uint_t itz) override {
        double zf = uzf[itz];

        // Enforce redshift if asked
        if (force_true_z) {
            zfit = {zf};
            nzfit = 1;
            nmodel = ntemplate;
            prior_table = replicate(1.0, prior_mag.size(), 1);
        } else {
            zfit = zfit_base;
            nzfit = nzfit_base;
            nmodel = nmodel_base;
            prior_table = prior_table_base;
        }

        // Pre-compute EAzY template fluxes and PSF moments
        bool compute_moments = false;
        if (force_true_z || eazy_q11.empty()) {
            eazy_q11.resize(nmodel);
            eazy_q12.resize(nmodel);
            eazy_q22.resize(nmodel);
            eazy_fvis.resize(nmodel);
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
                ascii::read_table(eazy_seds[it], rlam, rsed);
                rsed *= 1e-19;

                for (uint_t iz : range(nzfit)) {
                    // Apply IGM absorption
                    vec1d olam = rlam*(1.0 + zfit[iz]);
                    vec1d osed = rsed;
                    if (apply_igm) {
                        osed *= get_inoue(zfit[iz], rlam);
                    }

                    osed = cgs2uJy(olam, osed);
                    olam *= 1e-4;

                    // Normalize all templates to unit flux at 5500A rest-frame (NB: as in EAzY)
                    osed /= interpolate(osed, rlam, 5500.0);

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
                        eazy_q11.safe[iz*ntemplate+it] = sed2flux(
                            psf_filter.lam, psf_filter.res*mono_q11, olam, osed
                        )/fvis;
                        eazy_q12.safe[iz*ntemplate+it] = sed2flux(
                            psf_filter.lam, psf_filter.res*mono_q12, olam, osed
                        )/fvis;
                        eazy_q22.safe[iz*ntemplate+it] = sed2flux(
                            psf_filter.lam, psf_filter.res*mono_q22, olam, osed
                        )/fvis;
                        eazy_fvis.safe[iz*ntemplate+it] = fvis;
                    }
                }
            }
        }

        // Pre-compute template error function
        tpl_err.resize(nzfit, nband);
        for (uint_t izf : range(nzfit)) {
            for (uint_t l : range(nband)) {
                double lobs = filters[l].rlam/(1.0+zfit[izf]);
                if (lobs >= tpl_error_x[0] && lobs <= tpl_error_x[-1]) {
                    tpl_err.safe(izf,l) = interpolate(tpl_error_y, tpl_error_x, lobs);
                } else {
                    tpl_err.safe(izf,l) = 0.0;
                }
            }
        }

        // Reset averages
        m_ml.reset();
        m_ma.reset();

        // Create workspace
        prior.resize(nzfit);
        pz.resize(nmc, nzfit);
        chi2_best.resize(nmc);
        tpz.resize(nzfit);

        wflux.resize(nband);
        rflux.resize(nband);
        weight.resize(nband);
        wmodel.resize(ntemplate, nband);

        cache_bestz.resize(nmc);
        cache_bestc.resize(nmc,ntemplate);
        cache_pz.resize(nzfit);
        cache_pc.resize(nmodel);

        note("done.");
    }

    void initialize_cache() override {
        // Initialize arrays
        fitter_cache.allocate_column<float>("best_chi2",      niter, nmc);
        fitter_cache.allocate_column<uint_t>("best_z",        niter, nmc);
        if (cache_save_pmodel) {
            fitter_cache.allocate_column<float>("best_coef",  niter, nmc, ntemplate);
            fitter_cache.allocate_column<float>("pz",         niter, nzfit);
            fitter_cache.allocate_column<float>("pc",         niter, nmodel);
        }
        fitter_cache.allocate_column<float>("e1_eazy_ml",     niter);
        fitter_cache.allocate_column<float>("e2_eazy_ml",     niter);
        fitter_cache.allocate_column<float>("r2_eazy_ml",     niter);
        fitter_cache.allocate_column<float>("e1_eazy_ma",     niter);
        fitter_cache.allocate_column<float>("e2_eazy_ma",     niter);
        fitter_cache.allocate_column<float>("r2_eazy_ma",     niter);
        fitter_cache.allocate_column<float>("e1_eazy_ml_err", niter);
        fitter_cache.allocate_column<float>("e2_eazy_ml_err", niter);
        fitter_cache.allocate_column<float>("r2_eazy_ml_err", niter);
        fitter_cache.allocate_column<float>("e1_eazy_ma_err", niter);
        fitter_cache.allocate_column<float>("e2_eazy_ma_err", niter);
        fitter_cache.allocate_column<float>("r2_eazy_ma_err", niter);

        // Save meta data
        fitter_cache.write_columns("z_grid", zfit, "sed_grid", eazy_seds);
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
        to_fits("psf-mean-z"+to_string(iz)+"-eazy-ml.fits", m_ml,
            uzf, dndz, dndz_qu, dndz_sf, zml
        );
        to_fits("psf-mean-z"+to_string(iz)+"-eazy-ma.fits", m_ma,
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
    double min_mag_err = 0.0;
    uint_t nmc = 200;
    double dz = 0.01;
    uint_t seds_step = 5;
    double zfit_max = 7.0;
    double zfit_dz = 0.01;
    std::string prior_filter = "sdss-r";
    std::string prior_file = "/home/cschreib/programming/eazy-photoz/templates/prior_R_extend.dat";
    std::string template_error = "/home/cschreib/programming/eazy-photoz/templates/TEMPLATE_ERROR.eazy_v1.0";
    double template_error_amp = 0.5;
    bool apply_igm = true;
    bool force_true_z = false;
    bool no_noise = false;
    bool use_noline_library = false;
    bool use_egg_library = false;
    uint_t limited_set = 0;
    uint_t egg_sed_step = 1;
    bool write_cache = true;
    bool use_cache = true;
    bool cache_save_pmodel = true;
    uint_t iz = 5;

    read_args(argc, argv, arg_list(
        maglim, selection_band, filters, depths, nmc, min_mag_err, prior_filter, prior_file, dz,
        seds_step, apply_igm, zfit_max, zfit_dz, write_cache, use_cache, iz, template_error,
        template_error_amp, force_true_z, no_noise, use_noline_library, use_egg_library,
        limited_set, egg_sed_step, cache_save_pmodel
    ));

    eazy_averager pavg;
    pavg.write_cache = write_cache;
    pavg.use_cache = use_cache;

    // Setup survey
    egg::generator_options opts;
    opts.share_dir = "/home/cschreib/code/egg/share/";
    opts.filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
    opts.filter_flambda = true; // equivalent to FILTER_FORMAT=1
    opts.filter_photons = true; // equivalent to FILTER_FORMAT=1
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
    mopts.psf_file = "/home/cschreib/code/euclid_psf/psf-averager/psf-mono.fits";
    pavg.configure_mock(mopts);

    // Setup redshift fitting
    fitter_options fopts;
    fopts.prior_filter = prior_filter;
    fopts.prior_file = prior_file;
    fopts.zfit_max = zfit_max;
    fopts.zfit_dz = zfit_dz;
    fopts.template_error = template_error;
    fopts.template_error_amp = template_error_amp;
    fopts.apply_igm = apply_igm;
    fopts.force_true_z = force_true_z;
    fopts.use_noline_library = use_noline_library;
    fopts.use_egg_library = use_egg_library;
    fopts.egg_sed_step = egg_sed_step;
    fopts.limited_set = limited_set;
    fopts.cache_save_pmodel = cache_save_pmodel;
    pavg.configure_fitter(fopts);

    // Average PSF metrics
    // for (uint_t iz : range(pavg.zb)) {
    //     if (!pavg.average_redshift_bin(iz)) continue;
    // }

    if (!pavg.average_redshift_bin(iz)) return 1;

    return 0;
}
