#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    std::string filter_db_file = "filters.dat";
    auto filter_db = read_filter_db(filter_db_file);

    vec1s bands = {"euclid-vis", "sdss-u", "sdss-g", "sdss-r", "sdss-i", "sdss-z", "vista-Y", "vista-J", "vista-H"};
    append(bands, vec1s{"maiz-U","maiz-V","2mass-J"});
    append(bands, "mb-"+to_string_vector(rgen_step(400, 990, 10)));

    file::mkdir("psf-filters");
    std::ofstream odb("psf-filters/db.dat");
    auto write_filter = [&](const std::string& band, const filter_t& f) {
        fits::write_table("psf-filters/"+band+".fits", "lam", f.lam, "res", f.res);
        odb << band << "=" << band << ".fits" << std::endl;
    };

    filter_t vis_filter;
    for (auto b : bands) {
        filter_t filter;
        vif_check(get_filter(filter_db, b, filter),
            "could not find filter '", b, "' in '", filter_db_file, "'");

        // Truncate
        vec1u idg = where(filter.res/max(filter.res) > 1e-3);
        vif_check(!idg.empty(), "filter '", b, "' has no usable data");
        filter.lam = filter.lam[idg[0]-_-idg[-1]];
        filter.res = filter.res[idg[0]-_-idg[-1]];

        // Apply filter definition
        // Filter is defined such that it must integrate f_lambda and not f_nu
        // f_lambda*r(lambda) ~ f_nu*[r(lambda)/lambda^2]
        filter.res /= sqr(filter.lam);
        // Filter is defined such that it integrates photons and not energy
        // n(lambda)*r(lambda) ~ f(lambda)*[r(lambda)*lambda]
        filter.res *= filter.lam;

        // Re-normalize filter
        filter.res /= integrate(filter.lam, filter.res);
        filter.rlam = integrate(filter.lam, filter.lam*filter.res);

        write_filter(b, filter);

        if (b == "euclid-vis") {
            vis_filter = filter;
        }
    }

    vec1s libs = {"psf-mono", "psf-mono-f2", "psf-mono-f4", "psf-mono-f6", "psf-mono-f7", "psf-mono-f9",
        "psf-mono-coated3-0-0", "psf-mono-coated3-3-9", "psf-mono-coated3-5-5", "psf-mono-coated3-9-3",
        "psf-mono-coated3-9-9"};

    for (auto l : libs) {
        // Read monochromatic PSF library
        vec1d mono_lam, mono_w, mono_q11, mono_q12, mono_q22;
        fits::read_table(l+".fits",
            "lambda", mono_lam, "w", mono_w, "q11", mono_q11, "q12", mono_q12, "q22", mono_q22
        );

        // Match it to the PSF filter
        filter_t filter = vis_filter;
        mono_w   = interpolate(mono_w,   mono_lam, filter.lam);
        mono_q11 = interpolate(mono_q11, mono_lam, filter.lam);
        mono_q12 = interpolate(mono_q12, mono_lam, filter.lam);
        mono_q22 = interpolate(mono_q22, mono_lam, filter.lam);
        mono_lam = filter.lam;

        // Ignore data outside of adopted bandpass (450-950)
        {
            vec1u idi = where(mono_lam >= 0.450 && mono_lam <= 0.950);
            mono_w = mono_w[idi];
            mono_q11 = mono_q11[idi];
            mono_q12 = mono_q12[idi];
            mono_q22 = mono_q22[idi];
            mono_lam = mono_lam[idi];
            filter.res = filter.res[idi];
            filter.lam = filter.lam[idi];
        }

        // Include PSF weighting flux loss in PSF filter response
        filter.res *= mono_w;
        filter.res /= integrate(filter.lam, filter.res);

        write_filter("euclid-vis-"+l+"-ftot", filter);

        filter_t tmp = filter;
        tmp.res *= mono_q11;
        write_filter("euclid-vis-"+l+"-q11", tmp);
        tmp = filter;
        tmp.res *= mono_q12;
        write_filter("euclid-vis-"+l+"-q12", tmp);
        tmp = filter;
        tmp.res *= mono_q22;
        write_filter("euclid-vis-"+l+"-q22", tmp);
        tmp = filter;
        tmp.res *= filter.lam;
        write_filter("euclid-vis-"+l+"-lam", tmp);
    }

    return 0;
}
