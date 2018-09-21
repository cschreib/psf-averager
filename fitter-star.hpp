#ifndef FITTER_STAR_INCLUDED
#define FITTER_STAR_INCLUDED

#include "fitter-base.hpp"
#include "common.hpp"

class star_fitter : public fitter_base {
public :
    // BPZ PSF library
    vec1d star_q11, star_q12, star_q22;

    // Fit models
    uint_t nmodel = npos;
    uint_t ntemplate = npos;

    // Internal variables (constant per slice)
    vec2d tpl_flux;

    // Internal variables (workspace)
    struct workspace {
        vec1d wflux, rflux, weight;
        vec2d wmodel;
        vec1d wmm;
    };

    workspace global;

    // Config
    double rel_err = 0.0;
    std::string sed_fit;
    bool save_seds = false;

    // Internal variables
    vec1f save_sed_lambda;
    vec2f save_sed_fluxes;


    explicit star_fitter(filter_database& db, psf_moments& pm) : fitter_base(db, pm) {
        code_name = "star";
    }

    void do_read_options(program_arguments& opts) override {
        // Behavior
        double min_mag_err = 0.01;
        // Library
        sed_fit = "/home/cschreib/data/fits/templates/PHOENIX-ACES-AGSS-COND-2011-resample.fits";

        opts.read(arg_list(min_mag_err, sed_fit, save_seds));

        nband = filters.size();

        // Relative error on flux, sets minimum uncertainties
        rel_err = min_mag_err*(log(10.0)/2.5);

        vec1f lam;
        vec2f sed;
        fits::read_table(sed_fit, "lam", lam, "sed", sed);

        nmodel = ntemplate = sed.dims[0];

        if (save_seds) {
            save_sed_lambda = rgen_step(0.1,2.0,0.001);
        }

        tpl_flux.resize(nmodel, nband);
        star_q11.resize(nmodel);
        star_q12.resize(nmodel);
        star_q22.resize(nmodel);
        if (save_seds) {
            save_sed_fluxes.resize(nmodel, save_sed_lambda.size());
        }

        for (uint_t it : range(ntemplate)) {
            vec1d olam = lam;
            vec1d osed = sed(it,_);

            // Compute model fluxes
            for (uint_t l : range(nband)) {
                double flx = sed2flux(filters[l].lam, filters[l].res, olam, osed);
                if (!is_finite(flx)) {
                    // Falling out of the filter, assuming zero flux
                    flx = 0.0;
                }
                tpl_flux.safe(it,l) = flx;
            }

            // Resample SEDs on common grid
            if (save_seds) {
                save_sed_fluxes(it,_) = resample_sed(osed, olam, save_sed_lambda);
            }

            // Compute PSF moments
            psf.get_moments(olam, osed, star_q11.safe[it], star_q12.safe[it], star_q22.safe[it]);
        }
    }

    fit_result do_fit(uint_t iter, const vec1d& ftot) override {
        fit_result fr(nmc);

        workspace local;
        if (multi_threaded) {
            reset_workspace(local);
        }

        workspace& w = (multi_threaded ? local : global);

        vec1f berr;
        vec2f bflx;
        vec2f bflxm;
        vec2f tsed;
        if (save_seds) {
            berr = sqrt(phot_err2 + sqr(ftot*rel_err));
            bflx.resize(nmc,nband);
            bflxm.resize(nmc,nband);
            tsed.resize(nmc,save_sed_lambda.size());

            for (uint_t i : range(nmc))
            for (uint_t l : range(nband)) {
                bflx.safe(i,l) = ftot.safe[l] + (no_noise ? 0.0 : berr.safe[l]*mc.safe(i,l));
            }
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

        // Simulate SED measurements
        for (uint_t i : range(nmc)) {
            // Create noisy photometry
            if (!no_noise) {
                for (uint_t l : range(nband)) {
                    // In weighted units, the random perturbations have a sigma of unity
                    w.rflux.safe[l] = w.wflux.safe[l] + mc.safe(i,l);
                }
            }

            // Find best SED
            uint_t iml = npos;
            double bchi2 = finf;
            double bscale = 0;
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

                if (tchi2 < bchi2) {
                    bchi2 = tchi2;
                    bscale = scale;
                    iml = im;
                }
            }

            fr.chi2.safe[i] = bchi2;
            fr.psf_obs.safe[i] = metrics(star_q11.safe[iml], star_q12.safe[iml], star_q22.safe[iml]);

            if (save_seds) {
                for (uint_t il : range(save_sed_lambda)) {
                    tsed.safe(i,il) = bscale*save_sed_fluxes.safe(iml,il);
                }

                for (uint_t l : range(nband)) {
                    bflxm.safe(i,l) = bscale*tpl_flux.safe(iml,l);
                }
            }
        }

        if (save_seds) {
            fits::update_table("seds/sed_"+to_string(iter)+".fits", "sed_obs", tsed,
                "flux_obs", bflx, "err_obs", berr, "flux_fit", bflxm);
        }

        return fr;
    }

    std::string make_cache_hash() override {
        return hash(sed_fit, rel_err);
    }

    void reset_workspace(workspace& w) {
        w.wflux.resize(nband);
        w.rflux.resize(nband);
        w.weight.resize(nband);
        w.wmodel.resize(nmodel, nband);
        w.wmm.resize(nmodel);
    }

    void do_prepare_fit(double) override {
        // Create workspace
        if (!multi_threaded) {
            reset_workspace(global);
        }
    }

    void do_initialize_cache() override {}

    void save_individuals(const std::string& filename) override {}
};

#endif
