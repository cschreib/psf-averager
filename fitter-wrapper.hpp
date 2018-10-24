#ifndef FITTER_WRAPPER_INCLUDED
#define FITTER_WRAPPER_INCLUDED

#include "fitter-base.hpp"
#include "common.hpp"

class fitter_wrapper {
public :
    // Input SEDs
    uint_t niter = 0;
    vec1f errors;

    // Config
    uint_t nthread = 0;
    uint_t nmc = npos;
    fitter_base& fitter;
    bool no_noise = false;
    bool no_psf = true;
    double min_mag_err = 0.01;
    std::string catalog;
    std::string out;
    vec1u idf;

    explicit fitter_wrapper(fitter_base& f) : fitter(f) {}

    void read_options(program_arguments& opts) {
        // Averager options
        nmc = 1;
        min_mag_err = 0.01;
        no_noise = true;

        vec1s bands;

        // Read custom configuration from command line
        vec1s filters;
        vec1f depths;
        opts.read(arg_list(
            nmc, no_noise, min_mag_err, catalog, out, filters, no_psf, depths, nthread
        ));

        vif_check(!catalog.empty(), "please provide flux catalog in catalog=...");

        fits::input_table itbl(catalog);
        itbl.read_columns(ftable(bands));

        if (!filters.empty()) {
            idf.resize(filters.size());
            for (uint_t i : range(filters)) {
                idf[i] = where_first(bands == filters[i]);
                vif_check(idf[i] != npos, "could not find band '", filters[i], "' in catalog");
            }

            bands = filters;
        }

        fits::column_info ci;
        vif_check(itbl.read_column_info("flux", ci), "missing 'FLUX' column");
        niter = ci.dims[0];

        if (depths.empty()) {
            vif_check(itbl.read_column_info("flux_err", ci),
                "missing 'FLUX_ERR' column (or set depths=...)");
        } else {
            vif_check(depths.size() == bands.size(), "mismatch in depths and filters");
            errors = mag2uJy(depths)/10.0;
        }

        if (out.empty()) {
            out = file::remove_extension(catalog)+"-zfit-"+fitter.code_name+".fits";
        }

        // Overwrite
        {
            if (depths.empty()) {
                depths = replicate(99, bands.size());
                opts.write(arg_list(depths));
            }

            if (filters.empty()) {
                opts.write(arg_list(name(bands, "filters")));
            }
        }

        // Adjust
        if (no_noise && nmc != 1) {
            note("no_noise set, setting nmc=1");
            nmc = 1;
            opts.write(arg_list(nmc));
        }
    }

    void fit() {
        // Prepare fitter
        fitter.prepare_fit(niter, 1.0);

        fits::input_table itbl(catalog);

        {
            file::mkdir(file::get_directory(out));
            fits::output_table otbl(out);

            vec1s id;
            if (itbl.read_column("id", id)) {
                otbl.write_column("id", id);
            }

            if (no_noise) {
                otbl.allocate_column<float>("chi2_obs", niter);
                otbl.allocate_column<float>("z_obs", niter);
                otbl.allocate_column<float>("z_obsm", niter);
            } else {
                otbl.allocate_column<float>("chi2_obs", niter, nmc);
                otbl.allocate_column<float>("z_obs", niter, nmc);
                otbl.allocate_column<float>("z_obsm", niter, nmc);
            }

            if (!no_psf) {
                if (no_noise) {
                    otbl.allocate_column<float>("e1_obs", niter);
                    otbl.allocate_column<float>("e2_obs", niter);
                    otbl.allocate_column<float>("r2_obs", niter);
                    otbl.allocate_column<float>("e1_obsm", niter);
                    otbl.allocate_column<float>("e2_obsm", niter);
                    otbl.allocate_column<float>("r2_obsm", niter);
                } else {
                    otbl.allocate_column<float>("e1_obs", niter, nmc);
                    otbl.allocate_column<float>("e2_obs", niter, nmc);
                    otbl.allocate_column<float>("r2_obs", niter, nmc);
                    otbl.allocate_column<float>("e1_obsm", niter, nmc);
                    otbl.allocate_column<float>("e2_obsm", niter, nmc);
                    otbl.allocate_column<float>("r2_obsm", niter, nmc);
                }
            }
        }

        fits::table otbl(out);
        std::mutex fits_mutex;

        // auto pgi = progress_start(niter);
        // for (uint_t iter : range(niter)) {
        auto do_source = [&](uint_t iter) {
            vec1d flux, flux_err;

            {
                std::unique_lock<std::mutex> l(fits_mutex);

                itbl.read_elements("flux", flux, fits::at(iter,_));
                if (!idf.empty()) {
                    flux = flux[idf];
                }

                if (errors.empty()) {
                    itbl.read_elements("flux_err", flux_err, fits::at(iter,_));
                    flux_err = flux_err[idf];
                } else {
                    flux_err = errors;
                }
            }

            // Fit
            fit_result fr = fitter.do_fit(iter, flux, flux_err);

            {
                std::unique_lock<std::mutex> l(fits_mutex);

                if (no_noise) {
                    otbl.update_elements("chi2_obs", fr.chi2.safe[0], fits::at(iter));
                    otbl.update_elements("z_obs", fr.z_obs.safe[0], fits::at(iter));
                    otbl.update_elements("z_obsm", fr.z_obsm.safe[0], fits::at(iter));
                } else {
                    otbl.update_elements("chi2_obs", fr.chi2, fits::at(iter,_));
                    otbl.update_elements("z_obs", fr.z_obs, fits::at(iter,_));
                    otbl.update_elements("z_obsm", fr.z_obsm, fits::at(iter,_));
                }

                if (!no_psf) {
                    if (no_noise) {
                        otbl.update_elements("e1_obs", fr.psf_obs.safe[0].e1, fits::at(iter));
                        otbl.update_elements("e2_obs", fr.psf_obs.safe[0].e2, fits::at(iter));
                        otbl.update_elements("r2_obs", fr.psf_obs.safe[0].r2, fits::at(iter));
                        otbl.update_elements("e1_obsm", fr.psf_obsm.safe[0].e1, fits::at(iter));
                        otbl.update_elements("e2_obsm", fr.psf_obsm.safe[0].e2, fits::at(iter));
                        otbl.update_elements("r2_obsm", fr.psf_obsm.safe[0].r2, fits::at(iter));
                    } else {
                        otbl.update_elements("e1_obs", get_e1(fr.psf_obs), fits::at(iter,_));
                        otbl.update_elements("e2_obs", get_e2(fr.psf_obs), fits::at(iter,_));
                        otbl.update_elements("r2_obs", get_r2(fr.psf_obs), fits::at(iter,_));
                        otbl.update_elements("e1_obsm", get_e1(fr.psf_obsm), fits::at(iter,_));
                        otbl.update_elements("e2_obsm", get_e2(fr.psf_obsm), fits::at(iter,_));
                        otbl.update_elements("r2_obsm", get_r2(fr.psf_obsm), fits::at(iter,_));
                    }
                }
            }
        };

            // progress(pgi);
        // }

        thread::parallel_for pfor(nthread);
        pfor.verbose = true;
        pfor.progress_step = 1;
        pfor.update_rate = 2.0;
        pfor.chunk_size = 100;
        pfor.execute(do_source, niter);
    }
};

#endif
