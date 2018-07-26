#include "egg-analytic.hpp"

struct fitter_options {
    uint_t nmc = 1000;
    uint_t seed = 42;
    vec1f depths;
    std::string prior_filter;
    double min_mag_err = 0.05;
    double gauss_convolve = 0.021; // 0.03/sqrt(2) to match default BPZ
    uint_t ninterp = 2;
};

// Average PSF metrics
struct metrics {
    double q11 = 0, q12 = 0, q22 = 0;
    double e1 = 0, e2 = 0, r2 = 0;

    metrics() = default;

    metrics(double) {} // to enable integrate()

    explicit metrics(double t11, double t12, double t22) :
        q11(t11), q12(t12), q22(t22) {
        get_ellipticities();
    }

    void get_ellipticities() {
        r2 = q11+q22;
        e1 = (q11-q22)/r2;
        e1 = q12/r2;
    }

    metrics& operator *= (double norm) {
        q11 *= norm; q12 *= norm; q22 *= norm;
        e1 *= norm; e2 *= norm; r2 *= norm;
        return *this;
    }

    metrics& operator /= (double norm) {
        q11 /= norm; q12 /= norm; q22 /= norm;
        e1 /= norm; e2 /= norm; r2 /= norm;
        return *this;
    }

    metrics& operator += (const metrics& m) {
        q11 += m.q11; q12 += m.q12; q22 += m.q22;
        e1 += m.e1; e2 += m.e2; r2 += m.r2;
        return *this;
    }

    void reset() {
        operator*=(0.0);
    }
};

metrics operator* (metrics m, double norm) {
    return m *= norm;
}
metrics operator* (double norm, metrics m) {
    return m *= norm;
}
metrics operator+ (metrics m1, const metrics& m2) {
    return m1 += m2; // to enable integrate()
}

namespace std {
    template<>
    struct is_arithmetic<metrics> : std::true_type {}; // to enable integrate()
}

// Set of average PSFs for multiple categories of galaxies
struct metrics_set {
    metrics all, qu, sf;

    void reset() {
        all.reset();
        qu.reset();
        sf.reset();
    }

    void add(uint_t it, const metrics& m) {
        all += m;
        if (it == 0) {
            qu += m;
        } else {
            sf += m;
        }
    }

    void normalize(double ntot, double nqu, double nsf) {
        all /= ntot;
        qu /= nqu;
        sf /= nsf;
    }

    void to_fits(std::string filename,
        const vec1d& z, const vec1d& dndz, const vec1d& dndz_qu, const vec1d& dndz_sf) {

        fits::write_table(filename, "z", z, "dndz", dndz, "dndz_qu", dndz_qu, "dndz_sf", dndz_sf,
            "q11", all.q11, "q12", all.q12, "q22", all.q22,
            "e1", all.e1, "e2", all.e2, "r2", all.r2,
            "q11_qu", qu.q11, "q12_qu", qu.q12, "q22_qu", qu.q22,
            "e1_qu", qu.e1, "e2_qu", qu.e2, "r2_qu", qu.r2,
            "q11_sf", sf.q11, "q12_sf", sf.q12, "q22_sf", sf.q22,
            "e1_sf", sf.e1, "e2_sf", sf.e2, "r2_sf", sf.r2
        );
    }
};

metrics_set integrate(const vec1d& z, const vec<1,metrics_set>& v,
    const vec1d& dndz, const vec1d& dndz_qu, const vec1d& dndz_sf) {

    auto get_all = vectorize_lambda([](const metrics_set& n) {
        return n.all;
    });
    auto get_qu = vectorize_lambda([](const metrics_set& n) {
        return n.qu;
    });
    auto get_sf = vectorize_lambda([](const metrics_set& n) {
        return n.sf;
    });

    metrics_set m;
    m.all = integrate(z, get_all(v)*dndz);
    m.qu  = integrate(z, get_qu(v)*dndz_qu);
    m.sf  = integrate(z, get_sf(v)*dndz_sf);

    return m;
}

class psf_averager : public egg::generator {
public :
    // Total number of galaxies
    double ngal = 0, nqu = 0, nsf = 0;

    // Averages
    metrics_set m_ml; // maximum likelihood
    metrics_set m_ma; // marginalization
    metrics_set m_tr; // EGG (truth)

    // BPZ PSF library
    vec1s bpz_ssed, bpz_zz;
    vec1f bpz_z;
    vec1d bpz_q11, bpz_q12, bpz_q22;
    vec1d zfit;

    // EGG PSF library
    vec1s egg_ssed, egg_zz;
    vec1s egg_sz;
    vec1d egg_z;
    vec1d egg_q11, egg_q12, egg_q22;
    vec2d eggz_q11, eggz_q12, eggz_q22;

    // Internal variables
    vec1d wflux, rflux, weight;
    vec2d wmodel;
    vec2d mc;
    vec2d tpl_flux;
    vec1s bpz_seds;
    uint_t nmodel = npos;
    uint_t ntemplate = npos;
    uint_t nelliptical = npos;
    uint_t nspiral = npos;
    uint_t nirregular = npos;
    uint_t nzfit = npos;
    vec1d pmodel;
    vec1d prior;
    vec1d tt_prior;
    vec1d pz, pzc;
    vec1d wmm;

    vec1d zmeas_ml, zmeas_ma, zstack;

    bool just_count = false;
    uint_t niter = 0;
    uint_t iter = 0;

    progress_t pgi;

    // Config
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    double dz = 0.01;
    uint_t nmc = 1000;
    seed_t seed = make_seed(42);
    uint_t nband = npos;
    vec1f phot_err2;
    double rel_err = dnan;
    double gauss_convolve = dnan;
    double dzfit = dnan;
    vec1d gconv;
    uint_t wconv = npos;
    uint_t id_prior = npos;
    uint_t ninterp = 0;
    bool write_cache = true;

    // Cache
    struct {
        vec1u im, it, idisk, ibulge, ibt;
        vec1d pmag;
        vec2d pmodel, zmeas_ml, zmeas_ma;
    } zmeas_cache;

    psf_averager() : egg::generator() {}

    void configure_fitter(const fitter_options& opts) {
        nband = filters.size();

        phypp_check(is_any_of(opts.prior_filter, filter_names),
            "prior filter is not in the filter list ('", opts.prior_filter, "' not found')");
        id_prior = where_first(opts.prior_filter == filter_names) + 1;

        // Square of photometric error (Gaussian additive component)
        phot_err2 = sqr(mag2uJy(opts.depths)/10.0);
        phypp_check(phot_err2.size() == nband, "mismatch between filters (", nband, ") and depths (",
            phot_err2.size(), ")");

        // Relative error on flux, sets minimum uncertainties
        rel_err = opts.min_mag_err*(log(10.0)/2.5);

        // Cache random noise for re-use (same MC noise will be repeated for all galaxies)
        nmc = opts.nmc;
        seed = make_seed(opts.seed);
        mc = randomn(seed, nmc, nband);

        // Kernel used to convolve the P(z) to wash out too small fluctuations
        gauss_convolve = opts.gauss_convolve;

        // Interpolation of BPZ SEDs
        ninterp = opts.ninterp;

        // Read BPZ PSF library
        std::string bpz_filename = "BPZ/psfs-rebin2-cst.txt";
        ascii::read_table(bpz_filename, 0, bpz_ssed, bpz_zz, _, _, _, bpz_q11, bpz_q12, bpz_q22);

        // Resample library
        resample_library();

        // Make convolver kernel
        {
            uint_t npt = 18*gauss_convolve/dzfit;
            if (npt % 2 == 0) ++npt;
            wconv = npt/2;

            // This is a gaussian trucated between -9*sigma and +9*sigma
            vec1d kz = dzfit*(findgen(npt)-wconv);
            gconv = integrate_gauss(kz-dzfit/2.0, kz+dzfit/2.0, 0.0, gauss_convolve);
            gconv /= total(gconv);
        }
    }

    void set_priors(double mag) {
        // Create prior for three SED classes
        const vec1d alpha = {2.46, 1.81, 0.91};
        const vec1d z0 = {0.431, 0.390, 0.063};
        const vec1d km = {0.091, 0.0636, 0.123};

        const vec1d kt = {0.147, 0.450};
        const vec1d ft = vec1d{0.35, 0.50}/vec1u{nelliptical, nspiral};

        vec1d ptype(3);
        ptype.safe[0] = ft.safe[0]*exp(-kt.safe[0]*(mag - 20.0));
        ptype.safe[1] = ft.safe[1]*exp(-kt.safe[1]*(mag - 20.0));
        ptype.safe[2] = (1.0 - nelliptical*ptype.safe[0] - nspiral*ptype.safe[1])/nirregular;

        for (uint_t tt : range(3)) {
            // For each SED class, evaluate function of redshift
            double tprior = 0.0;
            for (uint_t iz : range(nzfit)) {
                uint_t im = iz*3 + tt;
                double zm = clamp(z0.safe[tt] + km.safe[tt]*(mag - 20.0), 0.01, 15.0);

                tt_prior.safe[im] = pow(zfit.safe[iz], alpha.safe[tt])
                    *exp(-pow(clamp(zfit.safe[iz]/zm, 0.0, 700.0), alpha.safe[tt]));

                tprior += tt_prior.safe[im];
            }

            double tprior2 = 0.0;
            for (uint_t iz : range(nzfit)) {
                uint_t im = iz*3 + tt;

                // Normalize
                tt_prior.safe[im] /= tprior;

                // Clip low probability wings
                if (tt_prior.safe[im] < 1e-2/nzfit) {
                    tt_prior.safe[im] = 0.0;
                } else {
                    tprior2 += tt_prior.safe[im];
                }
            }

            for (uint_t iz : range(nzfit)) {
                uint_t im = iz*3 + tt;

                // Normalize again
                tt_prior.safe[im] /= tprior2;
            }
        }

        // Expand to SED library
        for (uint_t iz : range(nzfit))
        for (uint_t it : range(ntemplate)) {
            uint_t im = iz*ntemplate + it;

            uint_t it0 = it/(ninterp+1);
            uint_t tt0 = (it0 < nelliptical ? 0 : it0 < nelliptical+nspiral ? 1 : 2);
            uint_t it1 = it0 + 1;
            uint_t tt1 = (it1 < nelliptical ? 0 : it1 < nelliptical+nspiral ? 1 : 2);
            uint_t itm = it % (ninterp+1);

            if (itm == 0 || tt0 == tt1) {
                // Non-interpolated SEDs get priors from their class,
                // and so do SEDs inteprolated in between two SEDs of the same class
                prior.safe[im] = tt_prior.safe[iz*3 + tt0];
            } else {
                // SEDs interpolated between different classes get interpolated priors
                double x = (itm+1.0)/(ninterp+1.0);
                prior.safe[im] = tt_prior.safe[iz*3 + tt0]*(1.0-x) - x*tt_prior.safe[iz*3 + tt1];
            }
        }

        prior.safe[_] = 1.0;
    }

    void on_generated(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        const double ftot = fdisk.safe[0] + fbulge.safe[0];

        if (ftot >= flim) {
            // We passed the magnitude cut!

            if (just_count) {
                ++niter;
                return;
            }

            // Compute true PSF moments
            // ------------------------

            // Compute flux-weighted B/T
            const double fbtn = fbulge.safe[0]/ftot;
            const double fbti = 1.0 - fbtn;

            // Add up PSFs
            m_tr.add(id_type, tngal*metrics(
                eggz_q11.safe[id_disk]*fbti + eggz_q11.safe[id_bulge]*fbtn,
                eggz_q12.safe[id_disk]*fbti + eggz_q12.safe[id_bulge]*fbtn,
                eggz_q22.safe[id_disk]*fbti + eggz_q22.safe[id_bulge]*fbtn
            ));

            ngal += tngal;
            if (id_type == 0) {
                nqu += tngal;
            } else {
                nsf += tngal;
            }

            // Now do BPZ
            // ----------

            // Create noise-free photometry
            for (uint_t l : range(nband)) {
                double bftot = fdisk.safe[l+1] + fbulge.safe[l+1];
                weight.safe[l] = 1.0/sqrt(phot_err2.safe[l] + sqr(bftot*rel_err));
                wflux.safe[l] = bftot*weight.safe[l];
            }

            // Create models
            for (uint_t im : range(nmodel)) {
                wmm.safe[im] = 0.0;
                for (uint_t l : range(nband)) {
                    wmodel.safe(im,l) = tpl_flux.safe(im,l)*weight.safe[l];
                    wmm.safe[im] += sqr(wmodel.safe(im,l));
                }
            }

            // Create priors
            double pmag = -2.5*log10(fdisk.safe[id_prior] + fbulge.safe[id_prior]) + 23.9;
            set_priors(pmag);

            // Simulate redshift measurements
            for (uint_t i : range(nmc)) {
                // Create noisy photometry
                for (uint_t l : range(nband)) {
                    // In weighted units, the random perturbations have a sigma of unity
                    rflux.safe[l] = wflux.safe[l] + mc.safe(i,l);
                }

                // Find best template and z
                uint_t iml = npos;
                double chi2_best = finf;
                for (uint_t im : range(nmodel)) {
                    double* wm = &wmodel.safe(im,0);
                    double* wf = &rflux.safe[0];

                    double wfm = 0.0;
                    for (uint_t l : range(nband)) {
                        wfm += wf[l]*wm[l];
                    }

                    double scale = wfm/wmm.safe[im];

                    double tchi2 = 0.0;
                    for (uint_t l : range(nband)) {
                        tchi2 += sqr(wf[l] - scale*wm[l]);
                    }

                    pmodel.safe[im] = tchi2;

                    if (tchi2 < chi2_best) {
                        chi2_best = tchi2;
                        iml = im;
                    }
                }

                // Create likelihood and apply prior
                // and marginalize over SEDs
                for (uint_t iz : range(nzfit)) {
                    pz.safe[iz] = 0.0;
                    for (uint_t it : range(ntemplate)) {
                        uint_t im = iz*ntemplate + it;
                        pmodel.safe[im] = exp(-0.5*(pmodel.safe[im] - chi2_best))*prior.safe[im];
                        pz.safe[iz] += pmodel.safe[im];
                    }
                }

                // Convolve with Gaussian to avoid too many peaks
                double pmax = 0.0;
                for (uint_t iz : range(nzfit)) {
                    // Convolve
                    uint_t i0 = iz > wconv ? iz - wconv : 0;
                    uint_t i1 = iz < nzfit-1-wconv ? iz + wconv : nzfit-1;
                    pzc.safe[iz] = 0.0;
                    for (uint_t tiz : range(i0, i1+1)) {
                        uint_t iconv = (int_t(tiz)-int_t(iz)) + wconv;
                        pzc.safe[iz] += pz.safe[tiz]*gconv.safe[iconv];
                    }

                    if (pzc.safe[iz] > pmax) {
                        pmax = pzc.safe[iz];
                    }
                }

                // Clip too low values, and then apply this back to the p(z,SED)
                double tprob = 0.0;
                for (uint_t iz : range(nzfit)) {
                    // Clip
                    if (pzc.safe[iz] < 1e-2*pmax) {
                        pzc.safe[iz] = 0;
                    }

                    // Update the p(z,SED)
                    for (uint_t it : range(ntemplate)) {
                        if (pz.safe[iz] > 0.0) {
                            uint_t im = iz*ntemplate + it;
                            pmodel.safe[im] *= pzc.safe[iz]/pz.safe[iz];
                            tprob += pmodel.safe[im];
                        }
                    }
                }

                // Maximum likelihood
                // ------------------

                m_ml.add(id_type,
                    (tngal/nmc)*metrics(bpz_q11.safe[iml], bpz_q12.safe[iml], bpz_q22.safe[iml])
                );

                if (write_cache) {
                    pzc /= total(pzc);
                    zmeas_ml.safe[i] = bpz_z.safe[iml];
                    zmeas_ma.safe[i] = total(zfit*pzc);
                    zstack += pmodel/nmc;
                }

                // Marginalization
                // ---------------

                metrics tm;
                for (uint_t im : range(nmodel)) {
                    tm += pmodel.safe[im]*metrics(
                        bpz_q11.safe[im], bpz_q12.safe[im], bpz_q22.safe[im]
                    );
                }

                m_ma.add(id_type, (tngal/tprob/nmc)*tm);

                // static bool saved = false;
                // if (!saved) {
                //     fits::write_table("testfit.fits",
                //         "lam", vec1d{0.35, 0.50, 0.65, 0.75, 0.85, 1.0, 1.2, 1.6},
                //         "flx", rflux/weight, "err", 1.0/weight,
                //         "true_flx", wflux/weight,
                //         "disk_flx", fdisk[1-_], "bulge_flx", fbulge[1-_],
                //         "best_model", scales[iml]*tpl_flux(iml,_),
                //         "iml", iml, "scales", scales,
                //         "pmodel", pmodel, "pzc", pzc, "pz", pz, "zfit", zfit, "prior", prior
                //     );

                //     print("SAVED, ibt=", ibt, " B/T=", bt[ibt]);
                //     saved = true;
                // }
            }

            // double n = now();
            // print(n - t);

            if (write_cache) {
                zmeas_cache.im.safe[iter] = id_mass;
                zmeas_cache.it.safe[iter] = id_type;
                zmeas_cache.idisk.safe[iter] = id_disk;
                zmeas_cache.ibulge.safe[iter] = id_bulge;
                zmeas_cache.ibt.safe[iter] = id_bt;
                zmeas_cache.pmag.safe[iter] = pmag;
                zmeas_cache.pmodel.safe(iter,_) = zstack;
                zmeas_cache.zmeas_ml.safe(iter,_) = zmeas_ml;
                zmeas_cache.zmeas_ma.safe(iter,_) = zmeas_ma;
                ++iter;
            }

            progress(pgi);

            // phypp_check(false, "bla");
        }
    }

    void resample_library() {
        // Sort by z then SED
        {
            vec1u ids = sort(bpz_zz+bpz_ssed);
            bpz_ssed  = bpz_ssed[ids];
            bpz_zz    = bpz_zz[ids];
            bpz_q11   = bpz_q11[ids];
            bpz_q12   = bpz_q12[ids];
            bpz_q22   = bpz_q22[ids];
        }

        // Reduce the size of the BPZ z grid
        vec1s ubz = unique_values_sorted(bpz_zz);
        vec1f tz;
        from_string(replace(ubz, "p", "."), tz);

        {
            uint_t skip_every = 10;
            double zmax_fit = 4.0;
            vec1u idr;
            uint_t k = 0;
            vec1s oubz = ubz;
            ubz.clear();
            for (uint_t i : range(oubz)) {
                if (k % skip_every != 0 || tz[i] > zmax_fit) {
                    append(idr, where(bpz_zz == oubz[i]));
                } else {
                    ubz.push_back(oubz[i]);
                }
                ++k;
            }

            uint_t old_size = bpz_ssed.size();
            inplace_remove(bpz_ssed, idr);
            inplace_remove(bpz_zz,   idr);
            inplace_remove(bpz_q11,  idr);
            inplace_remove(bpz_q12,  idr);
            inplace_remove(bpz_q22,  idr);

            note("shrunk BPZ PSF library from ", old_size, " to ", bpz_ssed.size(), " (",
                ubz.size(), " redshifts, ", bpz_ssed.size()/ubz.size(), " SEDs)");
        }

        // Find redshifts
        from_string(replace(ubz, "p", "."), zfit);
        nzfit = zfit.size();

        dzfit = zfit[1] - zfit[0];
        note("redshift step: ", dzfit);

        // Sort BPZ SEDs by color (red to blue)
        std::string sed_dir = "/home/cschreib/programming/bpz-1.99.3/SED/";
        bpz_seds = file::list_files(sed_dir, "*.sed");

        {
            vec1d color(bpz_seds.size());
            for (uint_t t : range(bpz_seds)) {
                vec1d rlam, rsed;
                ascii::read_table(sed_dir+bpz_seds[t], ascii::auto_skip('#'), rlam, rsed);
                rsed *= 1e-19;
                rsed = cgs2uJy(rlam, rsed);

                vec1u id1 = where(rlam > 2000 && rlam < 3000);
                uint_t id2 = where_first(rlam > 10000);
                color[t] = median(rsed[id1])/rsed[id2];
            }

            vec1u ids = sort(color);
            bpz_seds = bpz_seds[ids];
        }

        ntemplate = bpz_seds.size();
        nelliptical = 1;
        nspiral = 2;
        nirregular = ntemplate - nelliptical - nspiral;

        nmodel = ntemplate*nzfit;

        for (uint_t i : range(bpz_seds)) {
            note(" - ", i, ": ", bpz_seds[i]);
        }

        phypp_check(nmodel == bpz_ssed.size(),
            "mismatch between BPZ SEDs in directory (", ntemplate, ") and PSF library (",
            nzfit, "x", bpz_ssed.size()/nzfit, ")");

        // Find correspondence with PSF library, which was sorted by file name...
        {
            vec1s psf_seds = bpz_seds;
            inplace_sort(psf_seds);

            vec1u id_translate(psf_seds.size());
            for (uint_t it : range(ntemplate)) {
                id_translate[it] = where_first(psf_seds == bpz_seds[it]);
            }

            vec1u ids(nmodel);
            for (uint_t izf : range(nzfit))
            for (uint_t it : range(ntemplate)) {
                ids[izf*ntemplate + it] = where_first(
                    bpz_zz == ubz[izf] && bpz_ssed == to_string(id_translate[it])
                );
            }

            // Reshuffle
            // bpz_ssed is left alone, because we want it to remain sorted and map to
            // SED IDs in our new SED ordering, not the one of the PSF library...
            bpz_zz   = bpz_zz[ids];
            bpz_q11  = bpz_q11[ids];
            bpz_q12  = bpz_q12[ids];
            bpz_q22  = bpz_q22[ids];
        }

        // Interpolate the templates
        bpz_seds = sed_dir+bpz_seds;

        if (ninterp > 0) {
            std::string tmp_dir = "BPZ_interp/";
            file::mkdir(tmp_dir);

            // Save old PSF library
            vec1s oseds = bpz_seds;
            vec1s ossed = bpz_ssed;
            vec1d oq11 = bpz_q11;
            vec1d oq12 = bpz_q12;
            vec1d oq22 = bpz_q22;

            // Clear library
            bpz_seds.clear();
            bpz_ssed.clear();
            bpz_zz.clear();
            bpz_q11.clear();
            bpz_q12.clear();
            bpz_q22.clear();

            // Add first SED
            bpz_seds.push_back(oseds[0]);
            vec1u id0 = where(ossed == "0");
            append(bpz_ssed, replicate("0", nzfit));
            append(bpz_zz,   ubz);
            append(bpz_q11,  oq11[id0]);
            append(bpz_q12,  oq12[id0]);
            append(bpz_q22,  oq22[id0]);

            uint_t ntemplate_old = ntemplate;
            ntemplate = 1;

            for (uint_t it : range(ntemplate_old-1)) {
                // Read two adjacent SEDs
                vec1d rlam1, rsed1, rlam2, rsed2;
                ascii::read_table(oseds[it],   ascii::auto_skip('#'), rlam1, rsed1);
                ascii::read_table(oseds[it+1], ascii::auto_skip('#'), rlam2, rsed2);

                // Convert to uJy (BPZ interpolates fluxes in f_nu)
                rsed1 = cgs2uJy(rlam1, rsed1*1e-19);
                rsed2 = cgs2uJy(rlam2, rsed2*1e-19);

                // Merge wavelength axes in case they don't match and interpolate on merged axis
                vec1d clam = rlam1;
                append(clam, rlam2);
                inplace_sort(clam);
                clam = unique_values_sorted(clam);
                rsed1 = interpolate(rsed1, rlam1, clam);
                rsed2 = interpolate(rsed2, rlam2, clam);

                vec1u id1 = where(ossed == to_string(it+1));

                // Interpolate
                for (uint_t ii : range(ninterp)) {
                    double x = (ii+1.0)/(ninterp+1.0);
                    vec1d tsed1 = (1.0-x)*rsed1;
                    vec1d tsed2 = x*rsed2;
                    vec1d rsed = tsed1 + tsed2;
                    rsed = uJy2cgs(clam, rsed)*1e19;

                    std::string fname = file::remove_extension(file::get_basename(oseds[it]))
                        +"_"+to_string(ii)+"_"+to_string(ninterp)+".sed";
                    ascii::write_table(tmp_dir+fname, 18, clam, rsed);

                    // Insert SED in PSF library
                    bpz_seds.push_back(tmp_dir+fname);
                    std::string istr = align_right(to_string(ii), floor(log10(ninterp)), '0');
                    append(bpz_ssed, replicate(to_string(it)+"."+istr, nzfit));
                    append(bpz_zz,   ubz);
                    ++ntemplate;

                    for (uint_t iz : range(nzfit)) {
                        // We have to interpolate PSFs with a flux-weighted factor inside the VIS passband
                        vec1d tlam = clam*(1e-4*(1.0 + zfit[iz]));
                        double f1 = sed2flux(selection_filter.lam, selection_filter.res, tlam, tsed1);
                        double f2 = sed2flux(selection_filter.lam, selection_filter.res, tlam, tsed2);
                        double xf = f1/(f1 + f2);

                        bpz_q11.push_back(oq11[id0[iz]]*(1.0-xf) + xf*oq11[id1[iz]]);
                        bpz_q12.push_back(oq12[id0[iz]]*(1.0-xf) + xf*oq12[id1[iz]]);
                        bpz_q22.push_back(oq22[id0[iz]]*(1.0-xf) + xf*oq22[id1[iz]]);
                    }
                }

                // Added end SED
                bpz_seds.push_back(oseds[it+1]);
                append(bpz_ssed, replicate(to_string(it+1), nzfit));
                append(bpz_zz,   ubz);
                append(bpz_q11,  oq11[id1]);
                append(bpz_q12,  oq12[id1]);
                append(bpz_q22,  oq22[id1]);
                ++ntemplate;

                std::swap(id0, id1);
            }

            // Re-sort library, presently sorted by SED and z, into z and SED.
            vec1u ids = sort(bpz_zz+bpz_ssed);
            bpz_ssed  = bpz_ssed[ids];
            bpz_zz    = bpz_zz[ids];
            bpz_q11   = bpz_q11[ids];
            bpz_q12   = bpz_q12[ids];
            bpz_q22   = bpz_q22[ids];

            nmodel = nzfit*ntemplate;

            note("expanded PSF library from ", ntemplate_old, " to ", ntemplate, " SEDs");
        }

        // Finalize stuff
        from_string(replace(bpz_zz, "p", "."), bpz_z);

        // Create workspace
        tpl_flux.resize(nmodel, nband);
        wflux.resize(nband);
        rflux.resize(nband);
        weight.resize(nband);
        wmodel.resize(nmodel, nband);
        wmm.resize(nmodel);
        pmodel.resize(nmodel);
        prior.resize(nmodel);
        tt_prior.resize(nzfit*3);
        pz.resize(nzfit);
        pzc.resize(nzfit);

        if (write_cache) {
            zmeas_ml.resize(nmc);
            zmeas_ma.resize(nmc);
            zstack.resize(nmodel);
        }
    }

    bool read_egg_psfs(uint_t iz) {
        std::string zdir = "full_z"+to_string(iz)+"/";
        std::string filename = zdir+"psfs-rebin2-cst.txt";
        if (!file::exists(filename)) return false;

        // Read PSF library
        ascii::read_table(filename, 0, egg_ssed, egg_zz, _, _, _, egg_q11, egg_q12, egg_q22);

        phypp_check(!egg_zz.empty(), "empty PSF file '", filename, "'");

        // Get SED ID
        // Make sure SED ID is of the form xx-yy (i.e. convert 9 to 09)
        for (uint_t i : range(egg_ssed)) {
            vec1s spl = split(egg_ssed.safe[i], "-");
            if (spl[0].size() == 1) {
                spl[0] = '0'+spl[0];
            }
            if (spl[1].size() == 1) {
                spl[1] = '0'+spl[1];
            }

            egg_ssed.safe[i] = spl[0]+'-'+spl[1];
        }

        // Sort by z then SED
        {
            vec1u ids = sort(egg_zz+egg_ssed);
            egg_ssed = egg_ssed[ids];
            egg_zz  = egg_zz[ids];
            egg_q11 = egg_q11[ids];
            egg_q12 = egg_q12[ids];
            egg_q22 = egg_q22[ids];
        }

        // Find redshifts
        egg_sz = unique_values_sorted(egg_zz);
        from_string(replace(egg_sz, "p", "."), egg_z);

        return true;
    }

    bool average_redshift_bin(uint_t iz) {
        if (!read_egg_psfs(iz)) return false;

        uint_t ntz = floor((zb[iz+1]-zb[iz])/dz);
        vec1d uzf = zb[iz] + dz*dindgen(ntz);

        // Start averaging
        vec<1,metrics_set> zml(ntz), zma(ntz), ztr(ntz);
        vec1d dndz(ntz), dndz_qu(ntz), dndz_sf(ntz);
        // for (uint_t itz : range(1)) {
        for (uint_t itz : range(ntz)) {
            print(itz, "/", ntz);

            double zf = uzf.safe[itz];

            // Pre-compute BPZ template fluxes
            for (uint_t t : range(ntemplate)) {
                vec1d rlam, rsed;
                ascii::read_table(bpz_seds[t], ascii::auto_skip('#'), rlam, rsed);
                rsed *= 1e-19;

                for (uint_t izf : range(nzfit)) {
                    vec1d olam = rlam*(1.0 + zfit[izf]);
                    vec1d osed = cgs2uJy(olam, rsed);
                    olam *= 1e-4;

                    for (uint_t l : range(nband)) {
                        tpl_flux.safe(izf*ntemplate+t,l) = sed2flux(
                            filters[l].lam, filters[l].res, olam, osed);
                    }
                }
            }

            // Pre-select redshift slice in EGG PSF library
            {
                vec1u idz = where(egg_zz == egg_sz[min_id(abs(egg_z - zf))]);
                eggz_q11.resize(use.dims);
                eggz_q12.resize(use.dims);
                eggz_q22.resize(use.dims);

                for (uint_t i : idz) {
                    vec1s spl = split(egg_ssed.safe[i], "-");

                    uint_t iuv, ivj;
                    from_string(spl[0], iuv);
                    from_string(spl[1], ivj);

                    eggz_q11.safe(iuv,ivj) = egg_q11.safe[i];
                    eggz_q12.safe(iuv,ivj) = egg_q12.safe[i];
                    eggz_q22.safe(iuv,ivj) = egg_q22.safe[i];
                }
            }

            // Reset
            m_ml.reset();
            m_ma.reset();
            m_tr.reset();
            ngal = 0; nqu = 0; nsf = 0;

            // Pre-compute number of iterations
            niter = 0;
            just_count = true;
            generate(zf, dz);

            // Initialize cache
            if (write_cache) {
                iter = 0;
                zmeas_cache.im.resize(niter);
                zmeas_cache.it.resize(niter);
                zmeas_cache.idisk.resize(niter);
                zmeas_cache.ibulge.resize(niter);
                zmeas_cache.ibt.resize(niter);
                zmeas_cache.pmag.resize(niter);
                zmeas_cache.pmodel.resize(niter, nmodel);
                zmeas_cache.zmeas_ml.resize(niter, nmc);
                zmeas_cache.zmeas_ma.resize(niter, nmc);
            }

            // Compute average at that redshift
            just_count = false;
            pgi = progress_start(niter);
            generate(zf, dz);

            // Average quantities
            m_ml.normalize(ngal, nqu, nsf);
            m_ma.normalize(ngal, nqu, nsf);
            m_tr.normalize(ngal, nqu, nsf);

            // Store
            zml[itz] = m_ml;
            zma[itz] = m_ma;
            ztr[itz] = m_tr;

            dndz[itz]    = ngal/dz;
            dndz_qu[itz] = nqu/dz;
            dndz_sf[itz] = nsf/dz;

            if (write_cache) {
                fits::write_table("cache-bpz-z"+replace(to_string(format::fixed(format::precision(zf, 2))), ".", "p")+".fits",
                    ftable(zmeas_cache.im, zmeas_cache.it, zmeas_cache.idisk, zmeas_cache.ibulge,
                        zmeas_cache.ibt, zmeas_cache.pmag, zmeas_cache.pmodel, zmeas_cache.zmeas_ml,
                        zmeas_cache.zmeas_ma
                    )
                );

                zmeas_cache.im.clear();
                zmeas_cache.it.clear();
                zmeas_cache.idisk.clear();
                zmeas_cache.ibulge.clear();
                zmeas_cache.ibt.clear();
                zmeas_cache.pmag.clear();
                zmeas_cache.pmodel.clear();
                zmeas_cache.zmeas_ml.clear();
                zmeas_cache.zmeas_ma.clear();
            }
        }

        // Average over N(z)
        double ntot = integrate(uzf, dndz);
        double ntot_qu = integrate(uzf, dndz);
        double ntot_sf = integrate(uzf, dndz);
        m_ml = integrate(uzf, zml, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);
        m_ma = integrate(uzf, zma, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);
        m_tr = integrate(uzf, ztr, dndz/ntot, dndz_qu/ntot_qu, dndz_sf/ntot_sf);

        // Save
        m_ml.to_fits("psf-mean-z"+to_string(iz)+"-bpz-ml.fits",
            uzf, dndz, dndz_qu, dndz_sf
        );
        m_ma.to_fits("psf-mean-z"+to_string(iz)+"-bpz-ma.fits",
            uzf, dndz, dndz_qu, dndz_sf
        );
        m_tr.to_fits("psf-mean-z"+to_string(iz)+"-bpz-tr.fits",
            uzf, dndz, dndz_qu, dndz_sf
        );

        // fits::write_table("psf-mean-z"+to_string(iz)+"-bpz-ml.fits",
        //     "z", uzf, "dndz", dndz, "dndz_qu", dndz_qu, "dndz_sf", dndz_sf,
        //     "q11", ml_zq11, "q12", ml_zq12, "q22", ml_zq22,
        //     "e1", (ml_zq11-ml_zq22)/(ml_zq11+ml_zq22), "e2", 2*ml_zq12/(ml_zq11+ml_zq22),
        //     "r2", ml_zq11+ml_zq22
        // );

        // fits::write_table("psf-mean-z"+to_string(iz)+"-bpz-ma.fits",
        //     "z", uzf, "dndz", dndz, "dndz_qu", dndz_qu, "dndz_sf", dndz_sf,
        //     "q11", ma_zq11, "q12", ma_zq12, "q22", ma_zq22,
        //     "e1", (ma_zq11-ma_zq22)/(ma_zq11+ma_zq22), "e2", 2*ma_zq12/(ma_zq11+ma_zq22),
        //     "r2", ma_zq11+ma_zq22
        // );

        // fits::write_table("psf-mean-z"+to_string(iz)+"-bpz-tr.fits",
        //     "z", uzf, "dndz", dndz, "dndz_qu", dndz_qu, "dndz_sf", dndz_sf,
        //     "q11", tr_zq11, "q12", tr_zq12, "q22", tr_zq22,
        //     "e1", (tr_zq11-tr_zq22)/(tr_zq11+tr_zq22), "e2", 2*tr_zq12/(tr_zq11+tr_zq22),
        //     "r2", tr_zq11+tr_zq22
        // );

        return true;
    }
};

int phypp_main(int argc, char* argv[]) {
    // Survey definition
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";

    // Mock photometry
    vec1s filters = {"sdss-u", "sdss-g", "sdss-r", "sdss-i", "sdss-z", "vista-Y", "vista-J", "vista-H"};
    vec1f depths =  {24.5,     24.4,     24.1,     24.1,     23.7,     23.2,     23.2,      23.2}; // 10 sigma
    double min_mag_err = 0.03;
    uint_t nmc = 200;
    uint_t ninterp = 2;
    std::string prior_filter = "sdss-i";

    read_args(argc, argv, arg_list(
        maglim, selection_band, filters, depths, nmc, min_mag_err, prior_filter, ninterp
    ));

    psf_averager pavg;

    // Setup survey
    egg::generator_options opts;
    opts.share_dir = "/home/cschreib/code/egg/share/";
    opts.filter_db = "/home/cschreib/work_psf/scripts/averaging/db_bpz.dat";
    opts.filter_flambda = true;
    opts.filter_photons = true;
    opts.trim_filters = true;
    opts.selection_band = selection_band;
    opts.filters = filters;
    opts.maglim = maglim;
    opts.logmass_steps = 50;
    opts.bt_steps = 5;
    opts.seds_step = 5;
    pavg.initialize(opts);

    // Setup redshift fitting
    fitter_options fopts;
    fopts.depths = depths;
    fopts.nmc = nmc;
    fopts.min_mag_err = min_mag_err;
    fopts.prior_filter = prior_filter;
    fopts.ninterp = ninterp;
    pavg.configure_fitter(fopts);

    // Average PSF metrics
    // for (uint_t iz : range(pavg.zb)) {
    //     if (!pavg.average_redshift_bin(iz)) continue;
    // }

    if (!pavg.average_redshift_bin(5)) return 1;

    return 0;
}
