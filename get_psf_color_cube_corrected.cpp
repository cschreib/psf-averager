#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

#include "rebin.hpp"
#include "metrics.hpp"
#include "psf_moments.hpp"

int vif_main(int argc, char* argv[]) {
    std::string tinterp = "cst";
    uint_t nthread = 0;
    bool no_noise = false;
    std::string dir = "/home/cschreib/work_psf/psf-averager/";
    std::string model_cube = "color_cube_corrected.fits";
    uint_t step = 1;

    read_args(argc, argv, arg_list(
        name(tinterp, "interp"), nthread, no_noise, dir, model_cube, step
    ));

    dir = file::directorize(dir);

    enum class reconstruct {
        cst,
        mcspline
    } interp;

    if (tinterp == "cst") {
        interp = reconstruct::cst;
    } else if (tinterp == "mcspline") {
        interp = reconstruct::mcspline;
    } else {
        error("unknown reconstruction method '", tinterp, "'");
        return 1;
    }

    // Read Euclid PSF data
    filter_database fdb; {
        filter_options opts;
        opts.filter_db = dir+"filters.dat";
        opts.filter_flambda = true;
        opts.filter_photons = true;
        fdb.initialize(opts);
    }

    psf_moments psf(fdb); {
        program_arguments opts(vec1s{"psf_file="+dir+"psf-mono.fits"});
        psf.read_options(opts);
    }

    auto get_metrics = [&](const vec1d& mblam, const vec1d& mbflx) {
        // Interpolate NB photometry to 1nm resolution
        vec1d flx;
        if (interp == reconstruct::mcspline) {
            flx = rebin_mcspline(mbflx, mblam, psf.filter.lam);
        } else if (interp == reconstruct::cst) {
            flx = rebin_cst(mbflx, mblam, psf.filter.lam);
        }

        // Build PSF
        double q11, q12, q22;
        psf.get_moments_same_grid(flx, q11, q12, q22);

        return metrics(q11, q12, q22);
    };

    // Read trained data cube

    fits::input_table itbl(model_cube);

    uint_t         nmodel = 0;
    uint_t         nband = 0;
    uint_t         nsnr = 0;
    vec1d          model_prior;
    vec<2,metrics> model_psf;
    vec2d          model_fluxes;
    vec3d          model_mbfluxes;
    vec2u          model_id;
    vec1u          model_ri, model_iz;
    vec1s          fit_bands;
    vec1d          model_bcri;
    vec1d          model_bciz;
    vec1d          model_z;
    vec1d          mblam;
    vec1d          snrs;

    uint_t nri = 0, niz = 0;

    {
        note("read models");

        // Read models and broadband fluxes
        vec1s colors;
        itbl.read_columns(
            "id", model_id, "z", model_z, "weight", model_prior,
            "colors", colors, "bands", fit_bands, "snrs", snrs
        );

        model_prior /= median(model_prior);

        nmodel = model_id.dims[0];
        uint_t ncolor = model_id.dims[1];
        nband = ncolor+1;
        nsnr = snrs.size();

        note("found ", nmodel, " models with ", ncolor, " colors and ", nsnr, " SNRs");
        note("bands: ", fit_bands);

        uint_t cri = where_first(colors == "sdss-r/sdss-i");
        uint_t ciz = where_first(colors == "sdss-i/sdss-z");

        model_ri = model_id(_,cri);
        model_iz = model_id(_,ciz);

        vec1d blow, bup;
        itbl.read_column("bin"+to_string(cri)+"_low", blow);
        itbl.read_column("bin"+to_string(cri)+"_up", bup);
        model_bcri = 0.5*(blow + bup);
        nri = model_bcri.size();

        itbl.read_column("bin"+to_string(ciz)+"_low", blow);
        itbl.read_column("bin"+to_string(ciz)+"_up", bup);
        model_bciz = 0.5*(blow + bup);
        niz = model_bciz.size();

        // Read fluxes
        itbl.read_columns("fluxmb", model_mbfluxes, "flux", model_fluxes, "lambda_mb", mblam);

        note("construct PSF map");

        // Construct PSF map
        model_psf.resize(nsnr,nmodel);
        for (uint_t is : range(nsnr))
        for (uint_t i : range(nmodel)) {
            model_psf.safe(is,i) = get_metrics(mblam, model_mbfluxes.safe(is,i,_));
        }
    }


    note("read galaxy data");

    vec1s bands;
    vec1f depths;
    vec2f flux;
    fits::read_table("compiled_noisy.fits", ftable(bands, depths));

    vec1u idf(fit_bands.size());
    for (uint_t b : range(fit_bands)) {
        idf[b] = where_first(bands == fit_bands[b]);
        vif_check(idf[b] != npos, "could not find band '", fit_bands[b], "' in input catalog");
    }

    // depths = depths[idf]; // TODO: uncomment me at some point!

    uint_t br = where_first(fit_bands == "sdss-r");
    uint_t bi = where_first(fit_bands == "sdss-i");
    uint_t bz = where_first(fit_bands == "sdss-z");

    if (!no_noise) {
        fits::read_table("compiled_noisy.fits", ftable(flux));
    } else {
        fits::read_table("compiled.fits", ftable(flux));
    }

    uint_t ngal = flux.dims[0];

    struct package {
        metrics best;
        double z_best;
        uint_t id_best;
    };

    vec1d error2 = sqr(depths);

    auto get_ids = [&](double tri, double tiz, uint_t& iri, uint_t& iiz) {
        if (!is_finite(tri)) {
            iri = nri-1;
        } else {
            double x = clamp(interpolate(indgen<double>(nri), model_bcri, tri), 0, nri-1);
            iri = round(x);
        }

        if (!is_finite(tiz)) {
            iiz = niz-1;
        } else {
            double y = clamp(interpolate(indgen<double>(niz), model_bciz, tiz), 0, niz-1);
            iiz = round(y);
        }
    };

    auto get_psf = [&](const vec1d& flx) {
        package p;

        double chi2_expected = flx.size();
        double tprob = 0.0;
        double best_chi2 = dinf;
        double best_prob = 0.0;
        double best_scale = 0.0;

        const double* ff = &flx.safe[0];
        const double* fd = &error2.safe[0];

        // Estimate border of color range
        double tfr = max(flx[br], depths[br]);
        double tfi = max(flx[bi], depths[bi]);
        double tfz = max(flx[bz], depths[bz]);
        double tri = -2.5*log10(tfr/tfi);
        double tiz = -2.5*log10(tfi/tfz);
        double dri = (2.5/log(10.0))*sqrt(sqr(depths[br]/tfr) + sqr(depths[bi]/tfi));
        double diz = (2.5/log(10.0))*sqrt(sqr(depths[bi]/tfi) + sqr(depths[bz]/tfz));

        uint_t ri0, ri1, iz0, iz1;
        get_ids(tri - 6*dri, tiz - 6*diz, ri0, iz0);
        get_ids(tri + 6*dri, tiz + 6*diz, ri1, iz1);

        if (flx[br] <= depths[br]) {
            ri1 = nri-1;
        }
        if (flx[bi] <= depths[bi]) {
            iz1 = niz-1;
        }

        // Fit each model
        for (uint_t i : range(nmodel)) {
            const double* fm = &model_fluxes.safe(i,0);

            if (model_ri.safe[i] < ri0 || model_ri.safe[i] > ri1 ||
                model_iz.safe[i] < iz0 || model_iz.safe[i] > iz1) {
                // Model is too far away in color space, skip it!
                continue;
            }

            // Compute scaling factor and chi2
            double wff = 0.0;
            double wfm = 0.0;
            double wmm = 0.0;
            // double pe = 1.0;
            for (uint_t l : range(nband)) {
                double w2 = 1.0/(fd[l] + 1e-4*sqr(ff[l]));
                // pe *= w2;
                wff += sqr(ff[l])*w2;
                wfm += fm[l]*ff[l]*w2;
                wmm += sqr(fm[l])*w2;
            }

            // pe = sqrt(pe);

            double scale = wfm/wmm;
            double tchi2 = wff - scale*wfm; // - 2.0*log(pe);

            // Save best chi2
            if (tchi2 < best_chi2) {
                best_chi2 = tchi2;
                p.best = model_psf.safe[i];
                p.z_best = model_z.safe[i];
                p.id_best = i;
                best_scale = scale;
            }
        }

        double snr = best_scale*model_fluxes.safe(p.id_best,id_bsnr)/depths[id_bsnr];
        uint_t isnr = min_id(abs(log(snr/snrs)));
        p.best = model_psf.safe(isnr,p.id_best);
        p.z_best = model_z.safe[p.id_best];

        return p;
    };

    vec1d e1_best(ngal);
    vec1d z_best(ngal);
    vec1u id_best(ngal);

    auto get_source_psf = [&](uint_t j) {
        uint_t i = j*step;
        auto p = get_psf(flux(i,idf));
        e1_best[i] = p.best.e1;
        z_best[i] = p.z_best;
        id_best[i] = p.id_best;
    };

    note("estimating PSFs");

    thread::parallel_for pfor(nthread);
    pfor.verbose = true;
    pfor.progress_step = 1;
    pfor.update_rate = 2.0;
    pfor.chunk_size = 100;
    pfor.execute(get_source_psf, ngal/step);

    note("saving result");

    std::string suffix = "";
    if (!no_noise) suffix = "_noisy";

    fits::write_table("results_n_cor_"+suffix+".fits", ftable(
        e1_best, z_best, id_best
    ));

    return 0;
}
