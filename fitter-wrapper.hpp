#ifndef FITTER_WRAPPER_INCLUDED
#define FITTER_WRAPPER_INCLUDED

#include "fitter-base.hpp"
#include "common.hpp"

class fitter_wrapper {
public :
    // Individuals
    vec1s indiv_id;
    vec2f indiv_chi2;
    vec2f indiv_zml;
    vec2f indiv_zma;

    // Input SEDs
    vec2d obs_flux;
    vec2d obs_err;
    uint_t niter = 0;

    // Config
    uint_t nmc = npos;
    fitter_base& fitter;
    bool no_noise = false;
    double min_mag_err = 0.01;
    std::string catalog;
    std::string out;

    explicit fitter_wrapper(fitter_base& f) : fitter(f) {}

    void read_options(program_arguments& opts) {
        // Averager options
        nmc = 1;
        min_mag_err = 0.01;
        no_noise = true;

        vec1s bands;

        // Read custom configuration from command line
        opts.read(arg_list(nmc, no_noise, min_mag_err, catalog, out));

        vif_check(!catalog.empty(), "please provide flux catalog in catalog=...");
        fits::read_table(catalog, "id", indiv_id,
            "flux", obs_flux, "flux_err", obs_err, "bands", bands);

        if (out.empty()) {
            out = file::remove_extension(catalog)+"-zfit-"+fitter.code_name+".fits";
        }

        niter = obs_flux.dims[0];

        // Overwrite
        {
            vec1f depths = replicate(99, bands.size());
            bool no_psf = true;
            opts.write(arg_list(name(bands, "filters"), depths, no_psf));
        }

        // Adjust
        if (no_noise && nmc != 1) {
            note("no_noise set, setting nmc=1");
            nmc = 1;
            opts.write(arg_list(nmc));
        }

        indiv_id.resize(niter);
        indiv_chi2.resize(niter, nmc);
        indiv_zml.resize(niter, nmc);
        indiv_zma.resize(niter, nmc);
    }

    void fit() {
        // Prepare fitter
        fitter.prepare_fit(niter, 1.0);

        auto pgi = progress_start(niter);
        for (uint_t iter : range(niter)) {
            // Fit
            fit_result fr = fitter.do_fit(iter, obs_flux(iter,_), obs_err(iter,_));

            for (uint_t i : range(nmc)) {
                indiv_chi2.safe(iter,i) = fr.chi2.safe[i];
                indiv_zml.safe(iter,i) = fr.z_obs.safe[i];
                indiv_zma.safe(iter,i) = fr.z_obsm.safe[i];
            }

            progress(pgi);
        }

        file::mkdir(file::get_directory(out));
        fits::write_table(out,
            "id",       indiv_id,
            "chi2_obs", indiv_chi2,
            "z_obs",    indiv_zml,
            "z_obsm",   indiv_zma
        );
    }
};

#endif
