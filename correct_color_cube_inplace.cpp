#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    std::string input_cube = "color_cube.fits";
    std::string output_cube = "";

    vec1s dbands = {"sdss-r", "sdss-i", "sdss-z", "euclid-vis", "sdss-u", "sdss-g", "vista-Y", "vista-J", "vista-H"};
    vec1f ddepths = {24.1,     24.1,     23.7,     24.5,         24.5,     24.4,     23.2,      23.2,      23.2};

    // vec1d snrs = {5.0, 7.0, 10.0, 12.0, 15.0, 20.0, 40.0, 100.0};
    vec1d snrs = {10.0, 40.0, 100.0};

    // uint_t nmc = 100000;

    read_args(argc, argv, arg_list(input_cube, output_cube));

    uint_t nsnr = snrs.size();

    if (output_cube.empty()) {
        output_cube = file::remove_extension(input_cube)+"_corrected.fits";
    }

    // Copy original cube
    file::copy(input_cube, output_cube);

    // Read original flux data
    vec2d fluxmb_true;
    vec1s bands;
    vec2d flux_true;
    vec1d z_true;
    vec1d weight;
    fits::read_table(input_cube, "fluxmb", fluxmb_true, "bands", bands, "flux", flux_true,
        "weight", weight, "z", z_true);

    vec1u idd(bands.size());
    for (uint_t i : range(bands)) {
        idd[i] = where_first(dbands == bands[i]);
    }

    vec1d depths = mag2uJy(ddepths[idd])/10.0;
    uint_t idi = where_first(bands == "sdss-i");
    depths /= depths[idi];
    double rel_err = 0.01;

    uint_t nband = bands.size();
    uint_t nlib = fluxmb_true.dims[0];
    uint_t nmb = fluxmb_true.dims[1];
    vec3d fluxmb_cor(nsnr, nlib, nmb);
    vec2d z_cor(nsnr, nlib);
    vec2d tweight(nsnr, nlib);

    // auto seed = make_seed(55);

    // vec2d noise = randomn(seed, nmc, nband);
    // vec2d noiser = rel_err*randomn(seed, nmc, nband);
    // for (uint_t i : range(nmc)) {
    //     noise(i,_) *= depths;
    // }

    auto pg = progress_start(nsnr*nlib);
    for (uint_t is : range(snrs))
    for (uint_t i : range(nlib)) {
        double snr = snrs[is];
        double noise_amp = flux_true.safe(i,idi)/snr;
        vec1d error2 = sqr(depths*noise_amp);
        vec1d weight2 = 1.0/(error2 + sqr(rel_err*flux_true.safe(i,_)));

        const double* ff = &flux_true.safe(i,0);
        // const double* fd = &error2.safe[0];
        const double* fw = &weight2.safe[0];

        vec1d chi2_true(nlib);
        // vec1u best_id_mc(nmc);
        // vec1f best_chi2_mc = replicate(finf, nmc);
        for (uint_t j : range(nlib)) {
            const double* fm = &flux_true.safe(j,0);

            // Compute scaling factor and chi2
            double wff = 0.0;
            double wfm = 0.0;
            double wmm = 0.0;
            for (uint_t l : range(nband)) {
                // double w2 = 1.0/(fd[l] + 1e-4*sqr(ff[l]));
                double w2 = fw[l];
                wff += sqr(ff[l])*w2;
                wfm += fm[l]*ff[l]*w2;
                wmm += sqr(fm[l])*w2;
            }

            double scale = wfm/wmm;
            double tchi2 = wff - scale*wfm;

            chi2_true.safe[j] = tchi2;

            // for (uint_t im : range(nmc)) {
            //     const double* nn = &noise.safe(im,0);
            //     const double* nr = &noiser.safe(im,0);

            //     wff = 0.0;
            //     wfm = 0.0;
            //     for (uint_t l : range(nband)) {
            //         double w2 = fw[l];
            //         wff += sqr(ff[l]*(1.0 + nr[l]) + noise_amp*nn[l])*w2;
            //         wfm += fm[l]*(ff[l]*(1.0 + nr[l]) + noise_amp*nn[l])*w2;
            //     }

            //     scale = wfm/wmm;
            //     tchi2 = wff - scale*wfm;

            //     if (tchi2 < best_chi2_mc.safe[im]) {

            //     // double dchi2 = tchi2 - best_chi2_mc.safe[im];
            //     // if (dchi2 > 10) continue;

            //     // double pn = exp(-0.5*dchi2);
            //     // if (dchi2 < -10 || randomu(seed) < pn/(pn + 1.0)) {
            //         best_chi2_mc.safe[im] = tchi2;
            //         best_id_mc.safe[im] = j;
            //     }
            // }

            // progress(pg, 113);
        }

        vec1d p_true = exp(-0.5*(chi2_true - min(chi2_true)));
        p_true /= total(p_true);
        p_true *= weight[i];

        // fits::write_table("tmp.fits", ftable(p_true, best_id_mc));

        for (uint_t j : range(nlib)) {
            fluxmb_cor.safe(is,j,_) += fluxmb_true.safe(i,_)*p_true.safe[j];
            z_cor.safe(is,j) += z_true.safe[i]*p_true.safe[j];
            tweight.safe(is,j) += p_true.safe[j];
        }

        progress(pg, 11);
    }

    for (uint_t is : range(nsnr))
    for (uint_t il : range(nlib)) {
        fluxmb_cor.safe(is,il,_) /= tweight.safe(is,il);
        z_cor.safe(is,il) /= tweight.safe(is,il);
    }

    fits::update_table(output_cube, "fluxmb", fluxmb_cor, "z", z_cor, "snrs", snrs);

    return 0;
}
