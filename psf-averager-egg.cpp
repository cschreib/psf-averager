#include "fitter-base.hpp"
#include "psf-averager-photoz.hpp"

class null_fitter : public fitter_base {
public :

    explicit null_fitter(filter_database& db, psf_moments& pm) : fitter_base(db, pm) {
        code_name = "egg";
    }

    fit_result do_fit(uint_t iter, const vec1d& ftot) override {
        fit_result fr(nmc);
        return fr;
    }

    void do_prepare_fit(double zf) override {}

    std::string make_cache_hash() override {
        return "";
    }
};

int vif_main(int argc, char* argv[]) {
    uint_t iz = 5;
    bool allbins = false;

    filter_database db;
    psf_moments psf(db);
    null_fitter fitter(db, psf);
    psf_averager pavg(db, psf, fitter);

    // Read setup
    {
        program_arguments opts(argc, argv);

        // Override options
        {
            bool filter_flambda = true; // equivalent to FILTER_FORMAT=1
            bool filter_photons = true; // equivalent to FILTER_FORMAT=1
            bool trim_filters = true;
            bool write_cache = false;
            bool no_noise = true;
            uint_t nmc = 1;
            opts.write(arg_list(
                filter_flambda, filter_photons, trim_filters, write_cache, no_noise, nmc
            ));
        }

        // Default options
        {
            std::string filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
            bool write_individuals = false;
            bool write_averages = false;
            vec1s filters;
            opts.read(arg_list(filter_db, write_individuals, write_averages, filters));

            vec1d depths = replicate(29.0, filters.size());
            opts.write(arg_list(depths));
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

    // Average PSF metrics
    if (allbins) {
        for (uint_t tiz : range(pavg.zb.size()-1)) {
            if (!pavg.average_redshift_bin(tiz)) return 1;
        }
    } else {
        if (!pavg.average_redshift_bin(iz)) return 1;
    }

    return 0;
}
