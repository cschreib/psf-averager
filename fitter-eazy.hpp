#ifndef FITTER_EAZY_INCLUDED
#define FITTER_EAZY_INCLUDED

#include "fitter-base.hpp"
#include "common.hpp"

struct fitter_options {
    std::string prior_file;
    std::string prior_filter;
    std::string sed_dir;
    double zfit_max = 7.0;
    double zfit_dz = 0.01;
    std::string template_error;
    double template_error_amp = 0.5;
    double fit_tftol = 1e-4;
    bool apply_igm = true;
    bool use_noline_library = false;
    bool use_egg_library = false;
    bool use_eggpp_library = false;
    bool add_high_ew_template = false;
    uint_t egg_sed_step = 1;
    bool egg_sed_borders = false;
    bool egg_sed_add_center = false;
    bool force_true_z = false;
    uint_t limited_set = 0;
    bool cache_save_pmodel = true;
    bool indiv_save_coefs = false;
    double min_mag_err = 0.0;
};

class eazy_fitter : public fitter_base {
public :
    // Individuals
    vec3f indiv_coefs;
    vec3u indiv_seds;

    // EAzY PSF library
    vec1d eazy_q11, eazy_q12, eazy_q22, eazy_rlam;
    vec1d eazy_fvis;

    // EAzY SEDs
    vec3f save_sed_eazy_fluxes; // [ntemplate,nzfit,nlam]
    vec1f save_sed_lambda;      // [nlam]

    // Fit models
    vec1s eazy_seds;
    vec1d zfit_base;
    uint_t nzfit_base = npos;
    uint_t nmodel_base = npos;
    vec1d zfit;
    uint_t nmodel = npos;
    uint_t nzfit = npos;
    uint_t ntemplate = npos;
    uint_t nsed = npos;

    // Internal variables (constant)
    vec1d tpl_error_x, tpl_error_y;
    vec1d prior_mag;        // [nmag]
    vec2d prior_table_base; // [nmag,nzfit_base]

    // Internal variables (constant per slice)
    vec2d prior_table;      // [nmag,nzfit]
    vec2d tpl_flux;         // [nmodel,nband]
    vec2d tpl_err;          // [nzfit,nband]

    // Internal variables (workspace)
    struct workspace {
        vec1d prior;                 // [nzfit]
        vec2d chi2;                  // [nmc,nzfit]
        vec1d pz;                    // [nzfit]

        vec1d wflux, rflux, weight;  // [nband]
        vec2d wmodel;                // [ntemplate,nband]

        vec1u cache_bestz;           // [nmc]
        vec2d cache_bestc;           // [nmc,ntemplate]
        vec2d cache_pzc;             // [nmc,nmodel] or [nmc,nzfit*limset]
        vec2u cache_pzs;             // [nmc,nmodel] or [nmc,nzfit*limset]
    };

    workspace global;

    // Config
    double dzfit = dnan;
    std::string prior_file, prior_filter;
    uint_t id_prior = npos;
    bool apply_igm = true;
    bool force_true_z = false;
    bool use_noline_library = false;
    bool use_egg_library = false;
    bool cache_save_pmodel = true;
    bool indiv_save_coefs = false;
    uint_t egg_sed_step = 1;
    bool egg_sed_borders = false;
    bool egg_sed_add_center = false;
    bool add_high_ew_template = false;
    uint_t ihew = npos;
    uint_t limited_set = 0;
    double fit_tftol = 1e-4;
    bool save_seds = false;
    std::string sed_dir;
    double rel_err = dnan;


    explicit eazy_fitter(filter_database& db, psf_moments& pm) : fitter_base(db, pm) {
        code_name = "eazy";
    }

    void do_read_options(program_arguments& opts) override {
        fitter_options fopts;

        // Behavior
        fopts.zfit_max = 7.0;
        fopts.zfit_dz = 0.01;
        fopts.prior_filter = "sdss-r";
        fopts.prior_file = "/home/cschreib/programming/eazy-photoz/templates/prior_R_extend.dat";
        fopts.template_error = "/home/cschreib/programming/eazy-photoz/templates/TEMPLATE_ERROR.eazy_v1.0";
        fopts.template_error_amp = 0.5;
        fopts.apply_igm = true;
        fopts.force_true_z = false;
        fopts.min_mag_err = 0.0;

        // Library
        fopts.sed_dir = "/home/cschreib/code/euclid_psf/psf-averager/seds/";
        fopts.use_noline_library = false;
        fopts.use_egg_library = false;
        fopts.use_eggpp_library = false;
        fopts.add_high_ew_template = false;
        fopts.limited_set = 0;
        fopts.egg_sed_step = 1;
        fopts.egg_sed_borders = false;
        fopts.egg_sed_add_center = false;

        // Output
        fopts.cache_save_pmodel = true;
        fopts.indiv_save_coefs = false;

        opts.read(arg_list(
            fopts.prior_filter, fopts.prior_file, fopts.zfit_max, fopts.zfit_dz,
            fopts.template_error, fopts.template_error_amp, fopts.apply_igm, fopts.force_true_z,
            fopts.use_noline_library, fopts.use_egg_library, fopts.use_eggpp_library,
            fopts.egg_sed_step, fopts.limited_set, fopts.cache_save_pmodel, fopts.indiv_save_coefs,
            fopts.add_high_ew_template, fopts.egg_sed_borders, fopts.egg_sed_add_center,
            fopts.sed_dir, save_seds, fopts.min_mag_err
        ));

        if (fopts.use_eggpp_library) {
            fopts.use_egg_library = true;
        }

        initialize(fopts);
    }

    void initialize(const fitter_options& fopts) {
        apply_igm = fopts.apply_igm;
        force_true_z = fopts.force_true_z;
        use_noline_library = fopts.use_noline_library;
        use_egg_library = fopts.use_egg_library;
        egg_sed_step = fopts.egg_sed_step;
        egg_sed_borders = fopts.egg_sed_borders;
        egg_sed_add_center = fopts.egg_sed_add_center;
        limited_set = fopts.limited_set;
        fit_tftol = fopts.fit_tftol;
        cache_save_pmodel = fopts.cache_save_pmodel;
        indiv_save_coefs = fopts.indiv_save_coefs;
        sed_dir = fopts.sed_dir;
        add_high_ew_template = fopts.add_high_ew_template;

        if (limited_set > 3) {
            warning("'limited_set' can be at most equal to 3, set to zero to use all set");
            limited_set = 3;
        }

        if (add_high_ew_template && limited_set == 0) {
            add_high_ew_template = false;
            warning("'add_high_ew_template' disabled when 'limited_set=0'");
        }

        prior_filter = fopts.prior_filter;
        prior_file = fopts.prior_file;
        vif_check(is_any_of(prior_filter, bands),
            "prior filter is not in the filter list ('", prior_filter, "' not found in ", bands, ")");
        id_prior = where_first(prior_filter == bands);

        // Read template error file
        ascii::read_table(fopts.template_error, tpl_error_x, tpl_error_y);
        tpl_error_x /= 1e4;
        tpl_error_y *= fopts.template_error_amp;

        // Relative error on flux, sets minimum uncertainties
        rel_err = fopts.min_mag_err*(log(10.0)/2.5);

        // Setup redshift grid
        nzfit_base = ceil(fopts.zfit_max/fopts.zfit_dz);
        zfit_base = fopts.zfit_dz*indgen<double>(nzfit_base);

        if (save_seds) {
            file::mkdir("seds/");
            save_sed_lambda = rgen_step(0.1,2.0,0.001);
        }

        // Setup EAzY library
        make_sed_library(fopts);

        // Read prior
        read_prior();
    }

    void make_sed_library(const fitter_options& fopts) {
        // List SEDs
        if (use_egg_library) {
            std::string tsed_dir;
            if (fopts.use_eggpp_library) {
                if (fopts.use_noline_library) {
                    tsed_dir = sed_dir+"EGG++/";
                } else {
                    tsed_dir = sed_dir+"EGG++-lines/";
                }
                eazy_seds = tsed_dir+file::list_files(tsed_dir, "*.sed");
            } else {
                tsed_dir = sed_dir+"EGG/";
                eazy_seds = tsed_dir+file::list_files(tsed_dir, "*.dat");
            }

            inplace_sort(eazy_seds);

            vif_check(!eazy_seds.empty(), "no SED found for fitting in ", tsed_dir);

            vec1u suv(eazy_seds.size());
            vec1u svj(eazy_seds.size());
            for (uint_t i : range(eazy_seds)) {
                std::string sid = eazy_seds[i].substr((tsed_dir+"egg-").size(), 5);
                bool converted = from_string(sid.substr(0,2), suv[i]) &&
                                 from_string(sid.substr(3,2), svj[i]);
                vif_check(converted, "could not read UVJ ids from ", sid, " (",
                            sid.substr(0,2), ", ", sid.substr(3,2), ")");
            }

            if (egg_sed_borders) {
                vec1b keep(eazy_seds.size());
                vec1b center(eazy_seds.size());
                for (uint_t i : range(eazy_seds)) {
                    vec1u iduv = where(suv == suv[i]);
                    vec1u idvj = where(svj == svj[i]);
                    if (suv[i] == max(suv[idvj]) || suv[i] == min(suv[idvj]) ||
                        svj[i] == max(svj[iduv]) || svj[i] == min(svj[iduv])) {
                        keep[i] = true;
                    }
                    if (egg_sed_add_center && suv[i] == median(suv[idvj])) {
                        keep[i] = true;
                        center[i] = true;
                    }
                }

                if (egg_sed_step > 1) {
                    // Remove some SEDs to save time

                    // Sort by angle for borders
                    double muv = 0.5*(max(suv) + min(suv));
                    double mvj = 0.5*(max(svj) + min(svj));
                    vec1d th = atan2(suv - muv, svj - mvj);
                    vec1u idc = where(keep && !center);
                    idc = idc[sort(th[idc])];
                    keep[idc] = indgen(idc.size()) % egg_sed_step == 0u;

                    // Sort by VJ for center
                    idc = where(keep && center);
                    idc = idc[sort(svj[idc])];
                    keep[idc] = indgen(idc.size()) % egg_sed_step == 0u;
                }

                eazy_seds = eazy_seds[where(keep)];
                vif_check(!eazy_seds.empty(), "no SED left after keeping border and skipping");

            } else if (egg_sed_step > 1) {
                // Remove some SEDs to save time
                eazy_seds = eazy_seds[where((suv + svj) % egg_sed_step == 0u)];
                vif_check(!eazy_seds.empty(), "no SED left after skipping, reduce 'egg_sed_step'");
            }

            if (fopts.add_high_ew_template) {
                ihew = eazy_seds.size();
                eazy_seds.push_back(sed_dir+"eazy/erb2010_highEW.dat");
            }
        } else {
            if (use_noline_library) {
                eazy_seds = sed_dir+"eazy/" + vec1s{
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
                eazy_seds = sed_dir+"eazy/" + vec1s{
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

        ntemplate = nsed = eazy_seds.size();
        if (limited_set > 0) {
            nsed = limited_set;
            if (add_high_ew_template) {
                ++nsed;
            }
        }

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
            vif_check(spl.size() > 2 && spl[0] == "#" && spl[1] == "z",
                "ill-formed header line for prior file");

            spl = spl[2-_];
            vif_check(count(!from_string(spl, prior_mag)) == 0, "could not read prior magnitudes "
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

    void set_priors(workspace& w, const vec1d& ftot) {
        double mag = -2.5*log10(ftot.safe[id_prior]) + 23.9;

        uint_t imag = lower_bound(prior_mag, mag);
        if (imag == npos) {
            for (uint_t iz : range(nzfit)) {
                w.prior.safe[iz] = prior_table.safe(0,iz);
            }
        } else if (imag == prior_mag.size()-1) {
            for (uint_t iz : range(nzfit)) {
                w.prior.safe[iz] = prior_table.safe(imag,iz);
            }
        } else {
            double x = (mag - prior_mag.safe[imag])/(prior_mag.safe[imag+1] - prior_mag.safe[imag]);
            for (uint_t iz : range(nzfit)) {
                w.prior.safe[iz] = prior_table.safe(imag,iz)*(1-x) + x*prior_table.safe(imag+1,iz);
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
        set_priors(w, ftot);

        // Local variables
        matrix::mat2d alpha(ntemplate, ntemplate);
        vec1d beta(ntemplate);
        vec1d tcoefs(ntemplate);

        vec2f osed;
        if (save_seds) {
            osed.resize(nmc, save_sed_lambda.size());
        }

        // Simulate redshift measurements
        // We do it in two passes to save computations:
        // 1) Iterate over redshifts, build model matrices, and fit all the MC realizations
        // 2) Iterate over MC simulations and find p(z)

        // Step 1): compute chi2
        for (uint_t i : range(nmc)) {
            fr.chi2.safe[i] = finf;
        }

        for (uint_t iz : range(nzfit)) {
            // Create noise-free photometry
            for (uint_t l : range(nband)) {
                double bftot = ftot.safe[l];
                w.weight.safe[l] = 1.0/sqrt(phot_err2.safe[l] + sqr(bftot*(rel_err + tpl_err.safe(iz,l))));
                w.rflux.safe[l] = w.wflux.safe[l] = bftot*w.weight.safe[l];
            }

            // Create models
            for (uint_t it : range(ntemplate)) {
                uint_t im = iz*ntemplate + it;
                for (uint_t l : range(nband)) {
                    w.wmodel.safe(it,l) = tpl_flux.safe(im,l)*w.weight.safe[l];
                }
            }

            const double* wm = &w.wmodel.safe(0, 0);
            const double* wf = &w.rflux.safe[0];

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
                        w.rflux.safe[l] = w.wflux.safe[l] + mc.safe(i,l);
                    }

                    // Create beta vector
                    for (uint_t it : range(ntemplate)) {
                        // Ignore negative fluxes (NB: as in EAzY)
                        if (limited_set != 0 || w.rflux.safe[l] > 0.0) {
                            beta.safe[it] += w.rflux.safe[l]*wm[it*nband+l];
                        }
                    }
                }

                double tchi2;

                if (limited_set > 0) {
                    tchi2 = finf;

                    // Find best combination using 1, 2, or 3 templates from the library.
                    auto try_fit_single = [&](uint_t j) {
                        double coef = beta.safe[j]/alpha.safe(j,j);

                        if (coef >= 0) {
                            double chi2 = 0.0;
                            for (uint_t l : range(nband)) {
                                double tm = coef*wm[j*nband+l];
                                chi2 += sqr(wf[l] - tm);
                            }

                            if (chi2 < tchi2) {
                                tchi2 = chi2;
                                tcoefs = replicate(0.0, ntemplate);
                                tcoefs.safe[j] = coef;
                            }
                        }
                    };

                    auto try_fit_two = [&](uint_t j1, uint_t j2) {
                        double det = alpha.safe(j1,j1)*alpha.safe(j2,j2) - sqr(alpha.safe(j1,j2));

                        double coef1 =  beta.safe[j1]*alpha.safe(j2,j2) - beta.safe[j2]*alpha.safe(j1,j2);
                        double coef2 = -beta.safe[j1]*alpha.safe(j1,j2) + beta.safe[j2]*alpha.safe(j1,j1);

                        coef1 /= det;
                        coef2 /= det;

                        if (coef1 >= 0 && coef2 >= 0) {
                            double chi2 = 0.0;
                            for (uint_t l : range(nband)) {
                                double tm = coef1*wm[j1*nband+l] + coef2*wm[j2*nband+l];
                                chi2 += sqr(wf[l] - tm);
                            }

                            if (chi2 < tchi2) {
                                tchi2 = chi2;
                                tcoefs = replicate(0.0, ntemplate);
                                tcoefs.safe[j1] = coef1;
                                tcoefs.safe[j2] = coef2;
                            }
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
                            double chi2 = 0.0;
                            for (uint_t l : range(nband)) {
                                double tm = coef1*wm[j1*nband+l] + coef2*wm[j2*nband+l] + coef3*wm[j3*nband+l];
                                chi2 += sqr(wf[l] - tm);
                            }

                            if (chi2 < tchi2) {
                                tchi2 = chi2;
                                tcoefs = replicate(0.0, ntemplate);
                                tcoefs.safe[j1] = coef1;
                                tcoefs.safe[j2] = coef2;
                                tcoefs.safe[j3] = coef3;
                            }
                        }
                    };

                    for (uint_t i1 : range(ntemplate)) {
                        if (add_high_ew_template && i1 == ihew) continue;

                        try_fit_single(i1);

                        if (add_high_ew_template) {
                            try_fit_two(i1, ihew);
                        }

                        if (limited_set > 1) {
                            for (uint_t i2 : range(i1+1, ntemplate)) {
                                if (add_high_ew_template && i2 == ihew) continue;

                                try_fit_two(i1, i2);

                                if (add_high_ew_template) {
                                    try_fit_three(i1, i2, ihew);
                                }

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
                w.chi2.safe(i,iz) = tchi2;

                if (tchi2 < fr.chi2.safe[i]) {
                    fr.chi2.safe[i] = tchi2;
                    w.cache_bestz.safe[i] = iz;
                    for (uint_t it : range(ntemplate)) {
                        w.cache_bestc.safe(i,it) = tcoefs.safe[it];
                    }
                }

                if (limited_set > 0) {
                    // Only store the non-zero coefs and SEDs (saves space)
                    uint_t inz = 0;
                    for (uint_t it : range(nsed)) {
                        uint_t im = iz*nsed + it;

                        while (inz < ntemplate && tcoefs.safe[inz] == 0.0) {
                            ++inz;
                        }

                        if (inz == ntemplate) {
                            // Fewer templates used than maximum, set rest to zero
                            w.cache_pzs.safe(i,im) = 0;
                            w.cache_pzc.safe(i,im) = 0.0;
                        } else {
                            // Found template
                            w.cache_pzs.safe(i,im) = inz;
                            w.cache_pzc.safe(i,im) = tcoefs.safe[inz];
                            ++inz;
                        }
                    }
                } else {
                    // Store all coefs
                    for (uint_t it : range(ntemplate)) {
                        uint_t im = iz*ntemplate + it;
                        w.cache_pzc.safe(i,im) = tcoefs.safe[it];
                    }
                }
            }
        }

        // Step 2): compute p(z)
        for (uint_t i : range(nmc)) {
            fr.z_obs.safe[i] = zfit.safe[w.cache_bestz.safe[i]];
            fr.z_obsm.safe[i] = 0.0;

            // Create likelihood, apply prior, and marginalize over SEDs to build p(z)
            double tprob = 0.0;
            for (uint_t iz : range(nzfit)) {
                w.pz.safe[iz] = exp(-0.5*(w.chi2.safe(i,iz) - fr.chi2.safe[i]))*w.prior.safe[iz];
                tprob += w.pz.safe[iz];
            }

            // Apply p(z) to cache_pzc
            for (uint_t iz : range(nzfit)) {
                w.pz.safe[iz] /= tprob;
                fr.z_obsm.safe[i] += zfit.safe[iz]*w.pz.safe[iz];

                for (uint_t it : range(nsed)) {
                    uint_t im = iz*nsed + it;
                    w.cache_pzc.safe(i,im) *= w.pz.safe[iz];
                }
            }

            if (keep_individuals_in_memory && indiv_save_coefs) {
                if (limited_set > 0) {
                    // Only store the non-zero coefs and SEDs (saves space)
                    uint_t inz = 0;
                    for (uint_t it : range(nsed)) {
                        while (inz < ntemplate && w.cache_bestc.safe(i,inz) == 0.0) {
                            ++inz;
                        }

                        if (inz == ntemplate) {
                            // Fewer templates used than maximum, set rest to zero
                            indiv_seds(iter,i,it)       = 0;
                            indiv_coefs.safe(iter,i,it) = 0.0;
                        } else {
                            // Found template
                            indiv_seds(iter,i,it)       = inz;
                            indiv_coefs.safe(iter,i,it) = w.cache_bestc.safe(i,inz);
                            ++inz;
                        }
                    }
                } else {
                    // Store all coefs
                    for (uint_t it : range(ntemplate)) {
                        indiv_coefs.safe(iter,i,it) = w.cache_bestc.safe(i,it);
                    }
                }
            }

            if (save_seds) {
                uint_t iz = w.cache_bestz.safe[i];
                for (uint_t it : range(ntemplate)) {
                    double c = w.cache_bestc.safe(i,it);
                    if (c > 0.0) {
                        for (uint_t il : range(save_sed_lambda)) {
                            osed.safe(i,il) += c*save_sed_eazy_fluxes.safe(it,iz,il);
                        }
                    }
                }
            }

            // Maximum likelihood
            if (!no_psf) {
                uint_t iz = w.cache_bestz.safe[i];

                // Compute flux-weighted average moments for this model
                double q11 = 0.0, q12 = 0.0, q22 = 0.0, rlam = 0.0;
                double wtot = 0.0;
                for (uint_t it : range(ntemplate)) {
                    uint_t im = iz*ntemplate + it;

                    double tw = w.cache_bestc.safe(i,it)*eazy_fvis.safe[im];
                    wtot += tw;

                    q11 += tw*eazy_q11.safe[im];
                    q12 += tw*eazy_q12.safe[im];
                    q22 += tw*eazy_q22.safe[im];
                    rlam += tw*eazy_rlam.safe[im];
                }

                fr.psf_obs.safe[i] = metrics(q11/wtot, q12/wtot, q22/wtot, rlam/wtot);
            }


            // Marginalization
            if (!no_psf) {
                // Compute flux-weighted average moments for this model
                double q11 = 0.0, q12 = 0.0, q22 = 0.0, rlam = 0.0;
                double wtot = 0.0;

                for (uint_t iz : range(nzfit)) {
                    if (limited_set > 0) {
                        for (uint_t it : range(nsed)) {
                            uint_t imm = iz*nsed      + it;
                            uint_t im  = iz*ntemplate + w.cache_pzs.safe(i,imm);

                            double tw = w.cache_pzc.safe(i,imm)*eazy_fvis.safe[im];
                            wtot += tw;

                            q11 += tw*eazy_q11.safe[im];
                            q12 += tw*eazy_q12.safe[im];
                            q22 += tw*eazy_q22.safe[im];
                            rlam += tw*eazy_rlam.safe[im];
                        }
                    } else {
                        for (uint_t it : range(ntemplate)) {
                            uint_t im = iz*ntemplate + it;

                            double tw = w.cache_pzc.safe(i,im)*eazy_fvis.safe[im];
                            wtot += tw;

                            q11 += tw*eazy_q11.safe[im];
                            q12 += tw*eazy_q12.safe[im];
                            q22 += tw*eazy_q22.safe[im];
                            rlam += tw*eazy_rlam.safe[im];
                        }
                    }
                }

                fr.psf_obsm.safe[i] = metrics(q11/wtot, q12/wtot, q22/wtot, rlam/wtot);
            }
        }

        // Save SEDs
        if (save_seds) {
            fits::update_table("seds/sed_"+to_string(iter)+".fits", "sed_obs", osed);
        }

        // Save things in the cache
        if (write_cache) {
            // Could be executed concurrently, use mutex when in multithreading context
            auto lock = (multi_threaded ?
                std::unique_lock<std::mutex>(*cache_mutex_ptr) : std::unique_lock<std::mutex>());

            cache_ptr->update_elements("best_z",   w.cache_bestz, fits::at(iter,_));
            if (cache_save_pmodel) {
                cache_ptr->update_elements("best_coef", w.cache_bestc, fits::at(iter,_,_));
                cache_ptr->update_elements("pz_coef",   w.cache_pzc,   fits::at(iter,_,_));
                if (limited_set > 0) {
                    cache_ptr->update_elements("pz_sed", w.cache_pzs, fits::at(iter,_,_));
                }
            }
        }

        return fr;
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

    std::string make_cache_hash() override {
        return hash(use_noline_library, use_egg_library, egg_sed_step, limited_set,
            prior_file, prior_filter, apply_igm, zfit, eazy_seds, tpl_error_y, fit_tftol,
            bands, phot_err2, nmc, rel_err);
    }

    void reset_workspace(workspace& w) {
        w.prior.resize(nzfit);
        w.chi2.resize(nmc, nzfit);
        w.pz.resize(nzfit);;

        w.wflux.resize(nband);
        w.rflux.resize(nband);
        w.weight.resize(nband);
        w.wmodel.resize(ntemplate, nband);

        w.cache_bestz.resize(nmc);
        w.cache_bestc.resize(nmc,ntemplate);
        w.cache_pzc.resize(nmc,nzfit*nsed);
        if (limited_set > 0) {
            w.cache_pzs.resize(nmc,nzfit*nsed);
        }
    }

    void do_prepare_fit(double zf) override {
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
        if (!no_psf && (force_true_z || eazy_q11.empty())) {
            eazy_q11.resize(nmodel);
            eazy_q12.resize(nmodel);
            eazy_q22.resize(nmodel);
            eazy_fvis.resize(nmodel);
            eazy_rlam.resize(nmodel);
            compute_moments = true;
        }

        bool compute_fluxes = false;
        if (tpl_flux.empty()) {
            tpl_flux.resize(nmodel, nband);
            compute_fluxes = true;
        }

        if (compute_moments || compute_fluxes || save_seds) {
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
                        psf.get_moments(olam, osed,
                            eazy_q11.safe[iz*ntemplate+it], eazy_q12.safe[iz*ntemplate+it],
                            eazy_q22.safe[iz*ntemplate+it], eazy_fvis.safe[iz*ntemplate+it],
                            eazy_rlam.safe[iz*ntemplate+it]
                        );
                    }

                    if (save_seds) {
                        if (save_sed_eazy_fluxes.empty()) {
                            save_sed_eazy_fluxes.resize(ntemplate,nzfit,save_sed_lambda.size());
                        }

                        save_sed_eazy_fluxes(it,iz,_) = resample_sed(osed, olam, save_sed_lambda);
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

        // Create workspace
        if (!multi_threaded) {
            reset_workspace(global);
        }

        // Initialize individual arrays
        if (keep_individuals_in_memory && indiv_save_coefs) {
            indiv_coefs.resize(niter,nmc,nsed);
            if (limited_set > 0) {
                indiv_seds.resize(niter,nmc,nsed);
            }
        }
    }

    void do_initialize_cache() override {
        // Initialize arrays
        cache_ptr->allocate_column<uint_t>("best_z", niter, nmc);
        if (cache_save_pmodel) {
            cache_ptr->allocate_column<float>("best_coef", niter, nmc, ntemplate);
            cache_ptr->allocate_column<float>("pz",        niter, nzfit);
            cache_ptr->allocate_column<float>("pc",        niter, nmodel);
        }

        // Save meta data
        cache_ptr->write_columns("z_grid", zfit, "sed_grid", eazy_seds);
    }

    void save_individuals(const std::string& filename) override {
        // Write to disk the individual measurements
        fits::table otbl(filename);
        otbl.update_column("z_grid",   zfit);
        otbl.update_column("sed_grid", eazy_seds);

        if (indiv_save_coefs) {
            otbl.update_column("coef_obs", indiv_coefs);
            if (limited_set > 0) {
                otbl.update_column("seds_obs", indiv_seds);
            }
        }
    }
};

#endif
