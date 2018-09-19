#ifndef PSF_AVERAGER_FEEDBACK_INCLUDED
#define PSF_AVERAGER_FEEDBACK_INCLUDED

#include "fitter-base.hpp"
#include "common.hpp"
#include "metrics.hpp"

class psf_averager_feedback {
public :
    // Individuals
    metrics orig_tr;
    vec<1,metrics> indiv_tr;
    vec<2,metrics> indiv_ml, indiv_ma;

    // Monochromatic PSF library
    psf_moments& psf;

    // Input SEDs
    vec1f lam;
    vec2f seds;
    vec1f sed_tr;
    uint_t niter = 0;

    // Config
    std::string sed_file;
    vec1s bands;
    vec<1,filter_t> filters;
    uint_t nband = npos, nmc = npos;
    filter_database& filter_db;
    fitter_base& fitter;
    bool no_noise = false;

    // Cache
    std::string force_cache_id;
    std::string cache_dir;

    explicit psf_averager_feedback(filter_database& db, psf_moments& pm, fitter_base& f) :
        psf(pm), filter_db(db), fitter(f) {}

    void read_options(program_arguments& opts) {
        bands =        {"sdss-u", "sdss-g", "sdss-r", "sdss-i", "sdss-z", "vista-Y", "vista-J", "vista-H"};
        vec1f depths = {24.5,     24.4,     24.1,     24.1,     23.7,     23.2,     23.2,      23.2}; // 10 sigma

        // Averager options
        nmc = 200;
        no_noise = false;
        cache_dir = "";
        force_cache_id = "";
        sed_file = "";

        // Read custom configuration from command line
        opts.read(arg_list(
            name(bands, "filters"), depths, nmc, no_noise,
            cache_dir, name(force_cache_id, "cache_id"), sed_file
        ));

        vif_check(!sed_file.empty(), "please provide SED file in sed_file=...");

        // Adjust
        if (cache_dir.empty()) {
            cache_dir = file::get_directory(sed_file);
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

        nband = filters.size();
        cache_dir = file::directorize(cache_dir);

        fits::read_table(sed_file, "lam", lam, "sed_obs", seds, "sed_true", sed_tr);
        niter = seds.dims[0];

        double q11, q12, q22;
        psf.get_moments(lam, sed_tr, q11, q12, q22);
        orig_tr = metrics(q11, q12, q22);

        indiv_tr.resize(niter);
        indiv_ml.resize(niter, nmc);
        indiv_ma.resize(niter, nmc);
    }

    void average_seds() {
        // Prepare fitter
        fitter.prepare_fit(niter, 1.0);

        // Initialize cache
        std::string cache_id;
        if (force_cache_id.empty()) {
            cache_id = fitter.make_cache_hash();
        } else {
            cache_id = force_cache_id;
        }

        auto pgi = progress_start(niter);
        for (uint_t iter : range(niter)) {
            vec1d tlam = lam;
            vec1d tsed = seds(iter,_);

            {
                double q11, q12, q22;
                psf.get_moments(tlam, tsed, q11, q12, q22);
                indiv_tr[iter] = metrics(q11, q12, q22);
            }

            // Compute averages at that redshift
            vec1d ftot(nband);
            for (uint_t l : range(nband)) {
                ftot.safe[l] = sed2flux(filters[l].lam, filters[l].res, tlam, tsed);
            }

            fit_result fr = fitter.do_fit(iter, ftot);

            for (uint_t i : range(nmc)) {
                indiv_ml.safe(iter,i) = fr.psf_obs.safe[i];
                indiv_ma.safe(iter,i) = fr.psf_obsm.safe[i];
            }

            progress(pgi);
        }

        std::string base = file::get_basename(file::remove_extension(sed_file));
        std::string indiv_filename = cache_dir+"indiv-"+fitter.code_name+"-"+
            base+"-"+cache_id+".fits";
        note("individuals file: ", indiv_filename);

        file::mkdir(cache_dir);
        fits::write_table(indiv_filename,
            "e1_orig",  get_e1(orig_tr),
            "e2_orig",  get_e2(orig_tr),
            "r2_orig",  get_r2(orig_tr),
            "e1_true",  get_e1(indiv_tr),
            "e2_true",  get_e2(indiv_tr),
            "r2_true",  get_r2(indiv_tr),
            "e1_obs",   get_e1(indiv_ml),
            "e2_obs",   get_e2(indiv_ml),
            "r2_obs",   get_r2(indiv_ml),
            "e1_obsm",  get_e1(indiv_ma),
            "e2_obsm",  get_e2(indiv_ma),
            "r2_obsm",  get_r2(indiv_ma)
        );
    }
};

#endif
