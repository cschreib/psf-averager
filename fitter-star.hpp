#ifndef FITTER_STAR_INCLUDED
#define FITTER_STAR_INCLUDED

#include "fitter-base.hpp"
#include "common.hpp"

class star_fitter : public fitter_base {
public :
    // Star PSF library
    vec<1,metrics> star_psf;
    vec1d star_prob;

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
    bool prior_no_teff = false;

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

        std::string prior_stellar_locus = "/home/cschreib/work_psf/psf-averager/stellar_locus_nep.fits";
        opts.read(arg_list(min_mag_err, sed_fit, save_seds, prior_no_teff, prior_stellar_locus));

        nband = filters.size();

        // Relative error on flux, sets minimum uncertainties
        rel_err = min_mag_err*(log(10.0)/2.5);

        vec1f lam;
        vec2f sed;
        vec1f teff;
        fits::read_table(sed_fit, "lam", lam, "sed", sed);
        if (!prior_no_teff) {
            fits::read_table(sed_fit, "temp", teff);
        }

        nmodel = ntemplate = sed.dims[0];

        if (save_seds) {
            save_sed_lambda = rgen_step(0.1,2.0,0.001);
        }

        tpl_flux.resize(nmodel, nband);
        star_psf.resize(nmodel);
        if (save_seds) {
            save_sed_fluxes.resize(nmodel, save_sed_lambda.size());
        }

        filter_t rest_filter_r, rest_filter_i, rest_filter_z;
        rest_filter_r = filter_db.read_filter("cfht-r");
        rest_filter_i = filter_db.read_filter("cfht-i");
        rest_filter_z = filter_db.read_filter("cfht-z");

        vec1d star_ri(nmodel);
        vec1d star_iz(nmodel);

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
            double q11, q12, q22, rlam, fvis;
            psf.get_moments(olam, osed, q11, q12, q22, fvis, rlam);
            star_psf[it] = metrics(q11, q12, q22, rlam);

            // Compute colors
            double fr = sed2flux(rest_filter_r.lam, rest_filter_r.res, olam, osed);
            double fi = sed2flux(rest_filter_i.lam, rest_filter_i.res, olam, osed);
            double fz = sed2flux(rest_filter_z.lam, rest_filter_z.res, olam, osed);
            star_ri[it] = -2.5*log10(fr/fi);
            star_iz[it] = -2.5*log10(fi/fz);
        }

        // Attribute prior probability
        vec2d teff_min, teff_max, prob;
        vec2u use;
        vec2f bri, biz;
        fits::read_table(prior_stellar_locus, ftable(prob, bri, biz, use, teff_min, teff_max));

        star_prob.resize(nmodel);
        histogram2d(star_ri, star_iz, bri, biz, [&](uint_t iri, uint_t iiz, vec1u ids) {
            if (ids.empty() || !use(iri,iiz)) return;

            if (prior_no_teff) {
                star_prob[ids] = prob(iri,iiz)/ids.size();
            } else {
                vec1u idt = where(teff[ids] >= teff_min(iri,iiz) && teff[ids] <= teff_max(iri,iiz));
                if (idt.empty()) {
                    uint_t mid = min_id(min(
                        abs(teff[ids] - teff_min(iri,iiz)),
                        abs(teff[ids] - teff_max(iri,iiz))
                    ));
                    idt = {mid};
                }

                star_prob[ids[idt]] = prob(iri,iiz)/idt.size();
            }
        });

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
            double tprob = 0.0;
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

                double prob = exp(-0.5*(tchi2 - nband))*star_prob.safe[im];
                fr.psf_obsm.safe[i] += prob*star_psf.safe[im];
                tprob += prob;
            }

            fr.chi2.safe[i] = bchi2;
            fr.psf_obs.safe[i] = star_psf.safe[iml];
            fr.psf_obsm.safe[i] /= tprob;

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
