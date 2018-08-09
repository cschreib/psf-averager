#include "egg-analytic.hpp"
#include "metrics.hpp"

class psf_averager : public egg::generator {
public :
    // Average PSF metrics
    metrics m_all;
    double ngal = 0;

    // PSF library
    vec1s ssed, zz;
    vec1u sed_uv, sed_vj;
    vec1d mq11, mq12, mq22;
    std::array<uint_t,2> bz;

    // Internal variables
    uint_t ipsf_d = 0;
    uint_t ipsf_b = 0;
    double dq11 = dnan, dq12 = dnan, dq22 = dnan;
    double bq11 = dnan, bq12 = dnan, bq22 = dnan;

    psf_averager() : egg::generator() {}

    void update_psf_id(uint_t& ipsf, uint_t ised) {
        uint_t iuv = ised / use.dims[1];
        uint_t ivj = ised % use.dims[1];

        // We have sorted the PSF library by increasing UV and VJ.
        // So when the SED ID is changed, it is very likely that we
        // just went to the next PSF ID.
        while (sed_uv.safe[ipsf] != iuv || sed_vj.safe[ipsf] != ivj) {
            ++ipsf;
            if (ipsf == bz[1]) ipsf = bz[0];
        }
    }

    void on_disk_sed_changed(uint_t ised_d) override {
        update_psf_id(ipsf_d, ised_d);

        dq11 = mq11.safe[ipsf_d];
        dq12 = mq12.safe[ipsf_d];
        dq22 = mq22.safe[ipsf_d];
    }

    void on_bulge_sed_changed(uint_t ised_b) override {
        update_psf_id(ipsf_b, ised_b);

        bq11 = mq11.safe[ipsf_b];
        bq12 = mq12.safe[ipsf_b];
        bq22 = mq22.safe[ipsf_b];
    }

    void on_generated(uint_t im, uint_t it, uint_t ised_d, uint_t ised_b,
        uint_t ibt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {

        const double ftot = fdisk.safe[0] + fbulge.safe[0];

        if (ftot >= flim) {
            // We passed the magnitude cut!

            // Compute flux-weighted B/T
            const double fbtn = fbulge.safe[0]/ftot;
            const double fbti = 1.0 - fbtn;

            // Add up PSFs
            m_all += tngal*metrics(
                dq11*fbti + bq11*fbtn, dq12*fbti + bq12*fbtn, dq22*fbti + bq22*fbtn
            );

            ngal += tngal;
        }
    }

    bool average_redshift_bin(uint_t iz, std::string suffix) {
        std::string zdir = "full_z"+to_string(iz)+"/";
        std::string filename = zdir+"psfs-"+suffix+".txt";
        if (!file::exists(filename)) return false;

        // Read PSF library
        ascii::read_table(filename, ssed, zz, _, _, _, mq11, mq12, mq22);

        phypp_check(!zz.empty(), "empty PSF file '", filename, "'");

        // Get SED ID
        // Make sure SED ID is of the form xx-yy (i.e. convert 9 to 09)
        sed_uv.resize(ssed.size());
        sed_vj.resize(ssed.size());
        for (uint_t i : range(ssed)) {
            vec1s spl = split(ssed.safe[i], "-");
            from_string(spl.safe[0], sed_uv.safe[i]);
            from_string(spl.safe[1], sed_vj.safe[i]);
            if (spl[0].size() == 1) {
                spl[0] = '0'+spl[0];
            }
            if (spl[1].size() == 1) {
                spl[1] = '0'+spl[1];
            }

            ssed.safe[i] = spl[0]+'-'+spl[1];
        }

        // Sort by z then SED
        {
            vec1u ids = sort(zz+ssed);
            ssed = ssed[ids];
            zz  = zz[ids];
            sed_uv = sed_uv[ids];
            sed_vj = sed_vj[ids];
            mq11 = mq11[ids];
            mq12 = mq12[ids];
            mq22 = mq22[ids];
        }

        // Find redshifts
        const vec1s uz = unique_values_sorted(zz);
        vec1f uzf;
        from_string(replace(uz, "p", "."), uzf);

        const uint_t ntz = uz.size();
        uint_t nsed = 0;

        // Start averaging
        auto pg = progress_start(ntz);
        vec<1,metrics> zm(ntz);
        vec1d dndz(ntz);
        for (uint_t itz : range(ntz)) {
            double zf = uzf.safe[itz];
            double dz = (itz == 0 ? uzf.safe[1] - uzf.safe[0] : uzf.safe[itz] - uzf.safe[itz-1]);

            bz = equal_range(zz, uz[itz]);
            if (nsed == 0) {
                nsed = bz[1]-bz[0];
            } else {
                phypp_check(bz[1]-bz[0] == nsed, "incorrect number of SEDs for z", uz[itz],
                    " (got ", bz[1]-bz[0], ", expected ", nsed, ")");
            }

            ipsf_d = bz[0];
            ipsf_b = bz[0];

            // Reset
            m_all.reset();
            ngal = 0;

            // Compute average at that redshift
            generate(zf, dz);

            // Average quantities
            m_all /= ngal;

            // Store
            zm[itz] = m_all;
            dndz[itz] = ngal/dz;

            progress(pg);
        }

        // Average over N(z)
        double ntot = integrate(uzf, dndz);
        vec1d ndndz = dndz/ntot;
        m_all = integrate(uzf, zm*ndndz);

        fits::write_table(zdir+"psf-mean-"+suffix+".fits",
            "z", uzf, "dndz", dndz, "q11", get_q11(zm), "q12", get_q12(zm), "q22", get_q22(zm),
            "e1", get_e1(zm), "e2", get_e2(zm), "r2", get_r2(zm));

        return true;
    }
};

int phypp_main(int argc, char* argv[]) {
    // Redshift bins
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    uint_t nzbin = zb.size()-1;

    std::string suffix = "rebin2";
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";
    uint_t seds_step = 1;

    read_args(argc, argv, arg_list(suffix, maglim, selection_band, seds_step));

    psf_averager pavg;

    // Setup survey
    egg::generator_options opts;
    opts.filter_db = "/home/cschreib/work_psf/scripts/psf-averager/filters.dat";
    opts.share_dir = "/home/cschreib/code/egg/share/";
    opts.selection_band = selection_band;
    opts.filter_flambda = true;
    opts.filter_photons = true;
    opts.logmass_steps = 100;
    opts.bt_steps = 20;
    opts.seds_step = seds_step;
    opts.maglim = maglim;
    pavg.initialize(opts);

    // Average PSF metrics
    vec1d q11 = replicate(dnan, nzbin);
    vec1d q12 = replicate(dnan, nzbin);
    vec1d q22 = replicate(dnan, nzbin);
    vec1d e1 =  replicate(dnan, nzbin);
    vec1d e2 =  replicate(dnan, nzbin);
    vec1d r2 =  replicate(dnan, nzbin);

    for (uint_t iz : range(nzbin)) {
        if (!pavg.average_redshift_bin(iz, suffix)) continue;

        // Save average moments
        q11[iz] = pavg.m_all.q11;
        q12[iz] = pavg.m_all.q12;
        q22[iz] = pavg.m_all.q22;
        e1[iz] = pavg.m_all.e1;
        e2[iz] = pavg.m_all.e2;
        r2[iz] = pavg.m_all.r2;
    }

    // Recompute ellipticity from moments

    fits::write_table("psf-mean-"+suffix+".fits", ftable(zb, e1, e2, r2, q11, q12, q22));

    return 0;
}
