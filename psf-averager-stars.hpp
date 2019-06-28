#ifndef PSF_AVERAGER_STARS_INCLUDED
#define PSF_AVERAGER_STARS_INCLUDED

#include "fitter-base.hpp"
#include "common.hpp"
#include "metrics.hpp"

class psf_averager_stars {
public :
    // Individuals
    metrics orig_tr;
    vec<1,metrics> indiv_tr;
    vec<2,metrics> indiv_ml;
    vec<2,metrics> indiv_ma;
    vec1u indiv_id;
    vec1f indiv_mag;
    vec1f indiv_ri;
    vec1f indiv_iz;
    vec1f indiv_ngal;
    vec2f indiv_chi2;

    // PSF library
    psf_moments& psf;
    vec1d star_q11, star_q12, star_q22, star_rlam;
    vec1d star_ri, star_iz;
    vec1d star_prob;

    // Stellar locus
    std::string stellar_locus;
    vec2f prob;
    vec2f bri, biz;
    vec2b use;

    // Input SEDs
    vec2d tpl_flux;
    vec1d sel_flux;
    uint_t niter = 0;

    // Config
    std::string sed_gen;
    vec1s bands;
    vec<1,filter_t> filters;
    uint_t nband = npos, nmc = npos;
    filter_database& filter_db;
    fitter_base& fitter;
    bool no_noise = false;
    double mag_min = 18.5, mag_max = 22, mag_step = 0.5;
    double min_mag_err = 0.01;
    bool no_teff = false;
    bool save_seds = false;
    double zp_error = 0.0;
    std::string band_zp_error;
    uint_t id_band_zp_error = npos;
    vec1f mag_grid;
    std::string selection_band = "euclid-vis";
    filter_t selection_filter;
    uint_t nmag;
    uint_t ntemplate;

    // Internal variables
    vec1f blam;
    vec1f save_sed_lambda;
    vec2f save_sed_fluxes;

    // Cache
    std::string force_cache_id;
    std::string cache_dir;

    explicit psf_averager_stars(filter_database& db, psf_moments& pm, fitter_base& f) :
        psf(pm), filter_db(db), fitter(f) {}

    void read_options(program_arguments& opts) {
        bands =        {"sdss-u", "sdss-g", "sdss-r", "sdss-i", "sdss-z", "vista-Y", "vista-J", "vista-H"};
        vec1f depths = {24.5,     24.4,     24.1,     24.1,     23.7,     23.2,     23.2,      23.2}; // 10 sigma

        // Averager options
        nmc = 200;
        min_mag_err = 0.01;
        no_noise = false;
        cache_dir = "cache";
        force_cache_id = "";
        sed_gen = "/home/cschreib/data/fits/templates/PHOENIX-ACES-AGSS-COND-2011-resample.fits";
        stellar_locus = "/home/cschreib/work_psf/psf-averager/stellar_locus_nep.fits";

        // Read custom configuration from command line
        opts.read(arg_list(
            name(bands, "filters"), depths, nmc, no_noise,
            cache_dir, name(force_cache_id, "cache_id"), sed_gen, mag_min, mag_max, mag_step,
            selection_band, min_mag_err, stellar_locus, save_seds, no_teff, band_zp_error, zp_error
        ));

        vif_check(!sed_gen.empty(), "please provide SED generator library in sed_gen=...");

        // Adjust
        if (cache_dir.empty()) {
            cache_dir = file::get_directory(sed_gen);
            opts.write(arg_list(cache_dir));
        }
        if (no_noise && nmc != 1) {
            note("no_noise set, setting nmc=1");
            nmc = 1;
            opts.write(arg_list(nmc));
        }

        // Initialize
        filters.resize(bands.size());
        for (uint_t l : range(bands)) {
            filters[l] = filter_db.read_filter(bands[l]);
        }

        selection_filter = filter_db.read_filter(selection_band);

        // Grid
        vec2d teff_min, teff_max;
        fits::read_table(stellar_locus, ftable(prob, bri, biz, use, teff_min, teff_max));

        if (!band_zp_error.empty()) {
            id_band_zp_error = where_first(bands == band_zp_error);
            vif_check(id_band_zp_error != npos,
                "could not find band '", band_zp_error, "' to apply ZP error");
        }

        mag_grid = rgen_step(mag_min, mag_max, mag_step);
        nmag = mag_grid.size();

        nband = filters.size();
        cache_dir = file::directorize(cache_dir);

        vec1f lam;
        vec2f sed;
        vec1f teff;
        fits::read_table(sed_gen, "lam", lam, "sed", sed, "temp", teff);

        if (save_seds) {
            file::mkdir("seds/");
            save_sed_lambda = rgen_step(0.1,2.0,0.001);

            blam.resize(nband);
            for (uint_t i : range(nband)) {
                blam[i] = filters[i].rlam;
            }
        }

        // Compute colors of SEDs.
        filter_t rest_filter_r, rest_filter_i, rest_filter_z;
        rest_filter_r = filter_db.read_filter("cfht-r");
        rest_filter_i = filter_db.read_filter("cfht-i");
        rest_filter_z = filter_db.read_filter("cfht-z");

        ntemplate = sed.dims[0];

        tpl_flux.resize(ntemplate, nband);
        sel_flux.resize(ntemplate);
        star_q11.resize(ntemplate);
        star_q12.resize(ntemplate);
        star_q22.resize(ntemplate);
        star_rlam.resize(ntemplate);
        star_ri.resize(ntemplate);
        star_iz.resize(ntemplate);
        if (save_seds) {
            save_sed_fluxes.resize(ntemplate, save_sed_lambda.size());
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

            sel_flux[it] = sed2flux(selection_filter.lam, selection_filter.res, olam, osed);

            // Resample SEDs on common grid
            if (save_seds) {
                save_sed_fluxes(it,_) = resample_sed(osed, olam, save_sed_lambda);
            }

            // Compute PSF moments
            double fvis;
            psf.get_moments(olam, osed, star_q11.safe[it], star_q12.safe[it], star_q22.safe[it],
                fvis, star_rlam.safe[it]);

            // Compute colors
            double fr = sed2flux(rest_filter_r.lam, rest_filter_r.res, olam, osed);
            double fi = sed2flux(rest_filter_i.lam, rest_filter_i.res, olam, osed);
            double fz = sed2flux(rest_filter_z.lam, rest_filter_z.res, olam, osed);
            star_ri[it] = -2.5*log10(fr/fi);
            star_iz[it] = -2.5*log10(fi/fz);
        }

        // Attribute probability
        star_prob.resize(ntemplate);
        histogram2d(star_ri, star_iz, bri, biz, [&](uint_t iri, uint_t iiz, vec1u ids) {
            if (ids.empty() || !use(iri,iiz)) return;

            if (no_teff) {
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

        niter = ntemplate*nmag;

        indiv_id.resize(niter);
        indiv_mag.resize(niter);
        indiv_ri.resize(niter);
        indiv_iz.resize(niter);
        indiv_ngal.resize(niter);
        indiv_tr.resize(niter);
        indiv_chi2.resize(niter, nmc);
        indiv_ml.resize(niter, nmc);
        indiv_ma.resize(niter, nmc);
    }

    void average_seds() {
        // Prepare fitter
        fitter.prepare_fit(niter, 1.0);

        // Initialize cache
        std::string cache_id;
        if (force_cache_id.empty()) {
            cache_id = hash(fitter.make_cache_hash(), sed_gen, niter, mag_grid);
        } else {
            cache_id = force_cache_id;
        }

        auto pgi = progress_start(niter);
        uint_t iter = 0;
        for (uint_t it : range(ntemplate))
        for (uint_t im : range(nmag)) {
            double norm = e10(0.4*(23.9 - mag_grid[im]))/sel_flux[it];
            vec1d ftot = tpl_flux(it,_)*norm;

            // Save SEDs if asked
            if (save_seds) {
                vec1f tsed = norm*save_sed_fluxes(it,_);
                fits::update_table("seds/sed_"+to_string(iter)+".fits",
                    "lam", save_sed_lambda, "sed_true", tsed, "blam", blam, "flux_true", ftot
                );
            }

            if (id_band_zp_error != npos) {
                ftot[id_band_zp_error] *= e10(-0.4*zp_error);
            }

            indiv_tr[iter] = metrics(star_q11[it], star_q12[it], star_q22[it], star_rlam[it]);
            indiv_id[iter] = it;
            indiv_mag[iter] = mag_grid[im];
            indiv_ri[iter] = star_ri[it];
            indiv_iz[iter] = star_iz[it];
            indiv_ngal[iter] = star_prob[it]*e10(0.1157*mag_grid[im]);

            // Fit
            fit_result fi = fitter.do_fit(iter, ftot);

            for (uint_t i : range(nmc)) {
                indiv_chi2.safe(iter,i) = fi.chi2.safe[i];
                indiv_ml.safe(iter,i) = fi.psf_obs.safe[i];
                indiv_ma.safe(iter,i) = fi.psf_obsm.safe[i];
            }

            progress(pgi);
            ++iter;
        }

        std::string indiv_filename = cache_dir+"indiv-"+fitter.code_name+"-"+cache_id+".fits";
        note("individuals file: ", indiv_filename);

        vec1d e1bias = partial_mean(1, get_e1(indiv_ml)) - get_e1(indiv_tr);
        double mbias = total(e1bias*indiv_ngal)/total(indiv_ngal);
        double lbias = weighted_percentile(e1bias, indiv_ngal, 0.16);
        double ubias = weighted_percentile(e1bias, indiv_ngal, 0.84);
        note("bias: ", mbias*1e5, " +/- ", (ubias-lbias)/2.0*1e5);

        file::mkdir(cache_dir);
        fits::write_table(indiv_filename,
            "id",        indiv_id,
            "mag",       indiv_mag,
            "ri",        indiv_ri,
            "iz",        indiv_iz,
            "ngal",      indiv_ngal,
            "e1_true",   get_e1(indiv_tr),
            "e2_true",   get_e2(indiv_tr),
            "r2_true",   get_r2(indiv_tr),
            "rlam_true", get_rlam(indiv_tr),
            "chi2_obs",  indiv_chi2,
            "e1_obs",    get_e1(indiv_ml),
            "e2_obs",    get_e2(indiv_ml),
            "r2_obs",    get_r2(indiv_ml),
            "rlam_obs",  get_rlam(indiv_ml),
            "e1_obsm",    get_e1(indiv_ma),
            "e2_obsm",    get_e2(indiv_ma),
            "r2_obsm",    get_r2(indiv_ma),
            "rlam_obsm",  get_rlam(indiv_ma)
        );
    }
};

#endif
