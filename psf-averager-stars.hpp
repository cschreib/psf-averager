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
    vec1u indiv_id;
    vec1f indiv_mag;
    vec1f indiv_uv;
    vec1f indiv_vj;

    // PSF library
    psf_moments& psf;
    vec1d star_q11, star_q12, star_q22;
    vec1d star_uv, star_vj;

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
    vec1f mag_grid;
    std::string selection_band = "euclid-vis";
    filter_t selection_filter;
    uint_t nmag;
    uint_t ntemplate;

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

        // Read custom configuration from command line
        opts.read(arg_list(
            name(bands, "filters"), depths, nmc, no_noise,
            cache_dir, name(force_cache_id, "cache_id"), sed_gen, mag_min, mag_max, mag_step,
            selection_band, min_mag_err
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
        mag_grid = rgen_step(mag_min, mag_max, mag_step);
        nmag = mag_grid.size();

        nband = filters.size();
        cache_dir = file::directorize(cache_dir);

        vec1f lam;
        vec2f sed;
        fits::read_table(sed_gen, "lam", lam, "sed", sed);

        // Compute actual rest-frame colors of EGG SEDs.
        filter_t rest_filter_u, rest_filter_v, rest_filter_j;
        rest_filter_u = filter_db.read_filter("maiz-U");
        rest_filter_v = filter_db.read_filter("maiz-V");
        rest_filter_j = filter_db.read_filter("2mass-J");

        ntemplate = sed.dims[0];

        tpl_flux.resize(ntemplate, nband);
        sel_flux.resize(ntemplate);
        star_q11.resize(ntemplate);
        star_q12.resize(ntemplate);
        star_q22.resize(ntemplate);
        star_uv.resize(ntemplate);
        star_vj.resize(ntemplate);

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

            // Compute PSF moments
            psf.get_moments(olam, osed, star_q11.safe[it], star_q12.safe[it], star_q22.safe[it]);

            // Compute UVJ colors
            double fu = sed2flux(rest_filter_u.lam, rest_filter_u.res, olam, osed);
            double fv = sed2flux(rest_filter_v.lam, rest_filter_v.res, olam, osed);
            double fj = sed2flux(rest_filter_j.lam, rest_filter_j.res, olam, osed);
            star_uv[it] = -2.5*log10(fu/fv);
            star_vj[it] = -2.5*log10(fv/fj);
        }

        niter = ntemplate*nmag;

        indiv_id.resize(niter);
        indiv_mag.resize(niter);
        indiv_uv.resize(niter);
        indiv_vj.resize(niter);
        indiv_tr.resize(niter);
        indiv_ml.resize(niter, nmc);
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
            vec1d ftot = tpl_flux(it,_)*e10(0.4*(23.9 - mag_grid[im]))/sel_flux[it];

            indiv_tr[iter] = metrics(star_q11[it], star_q12[it], star_q22[it]);
            indiv_id[iter] = it;
            indiv_mag[iter] = mag_grid[im];
            indiv_uv[iter] = star_uv[it];
            indiv_vj[iter] = star_vj[it];
            indiv_vj[iter] = mag_grid[im];

            // Fit
            fit_result fr = fitter.do_fit(iter, ftot);

            for (uint_t i : range(nmc)) {
                indiv_ml.safe(iter,i) = fr.psf_obs.safe[i];
            }

            progress(pgi);
            ++iter;
        }

        std::string indiv_filename = cache_dir+"indiv-"+fitter.code_name+"-"+cache_id+".fits";
        note("individuals file: ", indiv_filename);

        file::mkdir(cache_dir);
        fits::write_table(indiv_filename,
            "id",       indiv_id,
            "mag",      indiv_mag,
            "uv",       indiv_uv,
            "vj",       indiv_vj,
            "e1_true",  get_e1(indiv_tr),
            "e2_true",  get_e2(indiv_tr),
            "r2_true",  get_r2(indiv_tr),
            "e1_obs",   get_e1(indiv_ml),
            "e2_obs",   get_e2(indiv_ml),
            "r2_obs",   get_r2(indiv_ml)
        );
    }
};

#endif
