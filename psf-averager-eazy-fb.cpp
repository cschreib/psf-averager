#include "fitter-eazy.hpp"
#include "psf-averager-feedback.hpp"

int vif_main(int argc, char* argv[]) {
    uint_t iz = 5;
    bool allbins = false;

    filter_database db;
    psf_moments psf(db);
    eazy_fitter fitter(db, psf);
    psf_averager_feedback pavg(db, psf, fitter);

    // Read setup
    {
        program_arguments opts(argc, argv);

        // Override options
        {
            bool filter_flambda = true; // equivalent to FILTER_FORMAT=1
            bool filter_photons = true; // equivalent to FILTER_FORMAT=1
            bool trim_filters = true;
            opts.write(arg_list(filter_flambda, filter_photons, trim_filters));
        }

        // Default options
        {
            std::string filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
            opts.read(arg_list(filter_db));
        }

        // Setup filter database
        db.read_options(opts);

        // Setup PSF moments
        psf.read_options(opts);

        // Setup survey
        pavg.read_options(opts);

        // Setup fitter
        fitter.read_options(opts);

        opts.read(arg_list(iz, allbins));
    }

    pavg.average_seds();

    return 0;
}
