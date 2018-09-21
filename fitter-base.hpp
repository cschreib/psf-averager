#ifndef FITTER_BASE_INCLUDED
#define FITTER_BASE_INCLUDED

#include <vif.hpp>
#include "filters.hpp"
#include "psf_moments.hpp"
#include "metrics.hpp"

using namespace vif;
using namespace vif::astro;

struct fit_result {
    explicit fit_result(uint_t n) {
        chi2.resize(n);
        z_obs.resize(n);
        z_obsm.resize(n);
        psf_obs.resize(n);
        psf_obsm.resize(n);
    }

    vec1f chi2;
    vec1f z_obs, z_obsm;
    vec<1,metrics> psf_obs, psf_obsm;
};

class fitter_base {
public :
    // Config
    filter_database& filter_db;
    psf_moments& psf;
    std::string code_name;
    vec1s bands;
    vec1f phot_err2;
    uint_t nband = npos;
    vec<1,filter_t> filters;
    bool multi_threaded = false;
    bool keep_individuals_in_memory = false;
    bool write_cache = false;
    bool no_noise = false;
    bool no_psf = false;
    uint_t nmc = 200;

    // Internal variables
    seed_t       seed = make_seed(42);
    vec2d        mc;
    fits::table* cache_ptr = nullptr;
    std::mutex*  cache_mutex_ptr = nullptr;
    uint_t       niter = 0;

    explicit fitter_base(filter_database& db, psf_moments& pm) : filter_db(db), psf(pm) {}

    virtual ~fitter_base() = default;

    void read_options(program_arguments& opts) {
        uint_t nthread = 0;
        vec1f depths;
        uint_t tseed = 42;
        std::string psf_file;

        opts.read(arg_list(
            name(bands, "filters"), nthread, keep_individuals_in_memory, nmc,
            no_noise, name(tseed, "seed"), write_cache, psf_file, depths, no_psf
        ));

        vif_check(nmc == 1 || nmc % 2 == 0, "nmc must be 1 or an even number, got ", nmc);

        multi_threaded = nthread != 0;

        // Find flux filters
        nband = bands.size();
        filters.resize(nband);
        for (uint_t l : range(bands)) {
            filters[l] = filter_db.read_filter(bands[l]);
        }

        // Square of photometric error (Gaussian additive component)
        phot_err2 = sqr(mag2uJy(depths)/10.0);
        vif_check(phot_err2.size() == nband, "mismatch between filters (", nband, ") and depths (",
            phot_err2.size(), ")");

        // Initialize random seed
        seed = make_seed(seed);

        // Forward configuration to implementation
        do_read_options(opts);
    }

    virtual void do_read_options(program_arguments& opts) {}

    virtual std::string make_cache_hash() {
        return "";
    }

    void reset_random() {
        if (nmc <= 1) {
            mc = randomn(seed, nmc, nband);
        } else {
            // Generate first half
            mc = randomn(seed, nmc/2, nband);
            // Make second half using negative of first half
            // This ensures the noise distribution is symmetric and unbiased
            append<0>(mc, -mc);
        }
    }

    void prepare_fit(uint_t tniter, double zf) {
        // Set number of iterations
        niter = tniter;

        // Initialize random numbers
        reset_random();

        // Initialize implementation
        do_prepare_fit(zf);
    }

    virtual void do_prepare_fit(double zf) = 0;

    void initialize_cache(fits::table& cache, std::mutex& cache_mutex) {
        cache_ptr = &cache;
        cache_mutex_ptr = &cache_mutex;
        do_initialize_cache();
    }

    virtual void do_initialize_cache() {}

    virtual fit_result do_fit(uint_t iter, const vec1d& ftot) = 0;

    virtual void save_individuals(const std::string& filename) {}
};

#endif
