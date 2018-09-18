#ifndef PSF_MOMENTS_INCLUDED
#define PSF_MOMENTS_INCLUDED

#include <vif.hpp>
#include "filters.hpp"

using namespace vif;
using namespace vif::astro;

struct psf_moments {
    filter_database& filter_db;
    filter_t filter;
    vec1d mono_lam, mono_q11, mono_q12, mono_q22, mono_w;

    explicit psf_moments(filter_database& db) : filter_db(db) {}

    void read_options(program_arguments& opts) {
        std::string psf_file = "/home/cschreib/code/euclid_psf/psf-averager/psf-mono.fits";
        opts.read(arg_list(psf_file));

        filter = filter_db.read_filter("euclid-vis");

        // Read monochromatic PSF library
        fits::read_table(psf_file,
            "lambda", mono_lam, "w", mono_w, "q11", mono_q11, "q12", mono_q12, "q22", mono_q22
        );

        // Match it to the PSF filter
        mono_w   = interpolate(mono_w,   mono_lam, filter.lam);
        mono_q11 = interpolate(mono_q11, mono_lam, filter.lam);
        mono_q12 = interpolate(mono_q12, mono_lam, filter.lam);
        mono_q22 = interpolate(mono_q22, mono_lam, filter.lam);
        mono_lam = filter.lam;

        // Ignore data outside of adopted bandpass (450-950)
        mono_w[where(mono_lam < 0.450 || mono_lam > 0.950)] = 0.0;

        // Include PSF weighting flux loss in PSF filter response
        filter.res *= mono_w;
        filter.res /= integrate(filter.lam, filter.res);
    }

    void get_moments(const vec1d& tlam, const vec1d& tsed,
        double& q11, double& q12, double& q22, double& ftot) const {

        vif_check(!filter.res.empty(), "uninitialized PSF moments");

        ftot = sed2flux(filter.lam, filter.res, tlam, tsed);
        q11  = sed2flux(filter.lam, filter.res*mono_q11, tlam, tsed)/ftot;
        q12  = sed2flux(filter.lam, filter.res*mono_q12, tlam, tsed)/ftot;
        q22  = sed2flux(filter.lam, filter.res*mono_q22, tlam, tsed)/ftot;
    }

    void get_moments(const vec1d& tlam, const vec1d& tsed,
        double& q11, double& q12, double& q22) const {

        double ftot;
        get_moments(tlam, tsed, q11, q12, q22, ftot);
    }
};

#endif
