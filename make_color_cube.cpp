#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

#include "tree.hpp"

int vif_main(int argc, char* argv[]) {
    vec1f ngal;
    vec1s bands_mb;
    vec1d lambda_mb;
    vec1f z;

    vec2f colors;
    vec2f fluxes;

    // Define the colors to use
    vec1s cbands = {"sdss-g","sdss-r","sdss-i","sdss-z","vista-Y","vista-J","vista-H"};

    // Define the width of color bins
    uint_t nmag = 6;
    uint_t nmag_riz = 12;
    double dmag = 0.1;
    double dmag_riz = 0.1;
    bool smart_bin = false;

    // Output file name
    std::string write_associations;
    std::string ofile = "color_cube.fits";

    read_args(argc, argv, arg_list(cbands, nmag, nmag_riz, dmag, dmag_riz, ofile, smart_bin,
        write_associations));

    vec1s cname;

    note("reading colors");

    {
        vec1s bands;
        vec2f flux;
        fits::read_table("compiled.fits", ftable(ngal, flux, bands, z));
        fits::read_table("compiled_mb.fits", "bands", bands_mb, "lambda", lambda_mb);

        // Read in broad band fluxes
        fluxes.resize(flux.dims[0], cbands.size());
        for (uint_t b : range(cbands)) {
            fluxes(_,b) = flux(_,where_first(bands == cbands[b]));
        }

        // Normalize to unit flux at reddest band
        for (uint_t i : range(ngal)) {
            fluxes(i,_) /= fluxes(i,cbands.size()-1);
        }

        // Compute colors
        cname.resize(cbands.size()-1);
        colors.resize(flux.dims[0], cname.size());
        for (uint_t b : range(cbands.size()-1)) {
            cname[b] = cbands[b]+"/"+cbands[b+1];
            colors(_,b) = -2.5*log10(fluxes(_,b)/fluxes(_,b+1));
        }
    }

    uint_t nb = bands_mb.size();
    uint_t ncol = cname.size();

    // Train
    note("design bins");

    vec<1,vec2f> bins(ncol);
    vec<1,vec1f> bcenter(ncol);
    vec<1,vec1f> bwidth(ncol);
    for (uint_t i : range(ncol)) {
        // Color is clamped to 3 to avoid long tails toward +inf for dropouts
        vec1f c = colors(_,i);
        c = min(c, 3.0);

        // Select good data points
        vec1u ids = where(is_finite(c));

        vec2f b;

        if (smart_bin) {
            // Sort data
            ids = ids[sort(c[ids])];

            // Define number of bins
            uint_t nm = nmag;
            if (cname[i] == "sdss-r/sdss-i" || cname[i] == "sdss-i/sdss-z") nm = nmag_riz;

            // Define bins with constant number density (refine where there are lots of galaxies)
            vec1d edges = {c[ids[0]]};
            double acc_ngal = 0.0;
            double ngal_step = total(ngal.safe[ids])/nm;
            uint_t ig = 0;
            while (ig < ids.size()) {
                uint_t isg = ids.safe[ig];
                acc_ngal += ngal.safe[isg];

                if (acc_ngal > ngal_step*edges.size()) {
                    edges.push_back(c.safe[isg]);
                }

                ++ig;
            }

            if (edges.size() != nm+1) {
                edges.push_back(c[ids.back()]);
            } else {
                edges.back() = c[ids.back()];
            }

            print(cname[i], ": ", edges.size()-1, " : ", edges);
            b = make_bins_from_edges(edges);
        } else {
            // Get bounds
            double mi = min(c);
            double ma = max(c);

            // Define width of bins
            double dm = dmag;
            if (cname[i] == "sdss-r/sdss-i" || cname[i] == "sdss-i/sdss-z") dm = dmag_riz;

            print(cname[i], ": ", mi, ", ", ma, ", ", b.dims[1]);
            b = make_bins_from_edges(rgen_step(mi, ma, dm));
        }

        bins[i] = b;
        bcenter[i] = bin_center(b);
        bwidth[i] = bin_width(b);
    }

    struct cell {
        double weight = 0.0;
        vec1d flux;
        vec1d fluxmb;
        double z = 0.0, zmin = +finf, zmax = -finf;
    };

    // Setup averaging (using a sorted tree structure for fast lookup and low memory)
    tree<float,cell> mf;
    mf.setup(bins);

    fits::input_table itbl("compiled_mb.fits");

    if (!write_associations.empty()) {
        fits::output_table otbl(write_associations);
        otbl.allocate_column<uint_t>("id_true", ngal.size(), bins.size());
    }

    fits::table atbl;
    if (!write_associations.empty()) {
        atbl.open(write_associations);
    }

    // Training data is loaded in chunks to avoid loading the entire data set in memory
    // at once, which could overflow the memory
    uint_t nchunk = 1e5;
    uint_t nstep = ceil(ngal.size()/float(nchunk));

    note("build tree");

    auto pg = progress_start(nstep);
    for (uint_t k : range(nstep)) {
        uint_t i0 = k*nchunk;
        uint_t i1 = (k+1)*nchunk;
        if (i1 > ngal.size()) i1 = ngal.size();
        --i1;

        vec1u cids = indgen(i1 - i0 + 1) + i0;

        // Grab fluxes from the file and store
        for (uint_t i : cids) {
            vec1f f;
            itbl.read_elements("flux", f, fits::at(i,_));
            f /= mean(f);

            // Remove invalid data (z=0 galaxies typically)
            if (!is_finite(f[0])) continue;

            // Find leaf of the tree with this color
            vec1u ik = mf.bin_key(colors(i,_));
            if (!write_associations.empty()) {
                atbl.update_elements("id_true", ik, fits::at(i,_));
            }

            auto& c = mf.insert_binned(ik);

            // Average stuff using number density as weighting
            double w = ngal[i];
            c.weight += w;

            // Medium band fluxes
            if (c.fluxmb.empty()) {
                c.fluxmb = f*w;
            } else {
                c.fluxmb += f*w;
            }

            // Broad band fluxes
            if (c.flux.empty()) {
                c.flux = fluxes(i,_)*w;
            } else {
                c.flux += fluxes(i,_)*w;
            }

            // Redshift
            c.z += z[i]*w;
            c.zmin = min(c.zmin, z[i]);
            c.zmax = max(c.zmax, z[i]);
        }

        progress(pg);
    }

    note("normalize tree");

    // Finalize averaging, dividing by total weight
    uint_t nflx = 0;
    uint_t nflxmb = 0;
    mf.for_each_leaf([&](const vec1u& k, cell& c) {
        if (nflx == 0) nflx = c.flux.size();
        if (nflxmb == 0) nflxmb = c.fluxmb.size();

        c.flux /= c.weight;
        c.fluxmb /= c.weight;
        c.z /= c.weight;
    });

    uint_t ncell = mf.size();

    // Save binned data
    note("save tree (", ncell, " cells)");


    // First save meta data and allocate space
    {
        fits::output_table otbl(ofile);
        otbl.write_columns("colors", cname, "bands_mb", bands_mb, "lambda_mb", lambda_mb,
            "bands", cbands);

        for (uint_t i : range(ncol)) {
            otbl.write_column("bin"+to_string(i)+"_low", bins[i](0,_));
            otbl.write_column("bin"+to_string(i)+"_up", bins[i](1,_));
        }

        otbl.allocate_column<double>("flux",   ncell, nflx);
        otbl.allocate_column<double>("fluxmb", ncell, nflxmb);
        otbl.allocate_column<double>("weight", ncell);
        otbl.allocate_column<double>("z",      ncell);
        otbl.allocate_column<double>("zmin",   ncell);
        otbl.allocate_column<double>("zmax",   ncell);
        otbl.allocate_column<short>("id",      ncell, ncol);
    }

    // Now save individual cells in the sorted tree one by one
    fits::table otbl(ofile);
    uint_t il = 0;
    mf.for_each_leaf([&](const vec1u& k, cell& c) {
        otbl.update_elements("id",     k,        fits::at(il,_));
        otbl.update_elements("flux",   c.flux,   fits::at(il,_));
        otbl.update_elements("fluxmb", c.fluxmb, fits::at(il,_));
        otbl.update_elements("weight", c.weight, fits::at(il));
        otbl.update_elements("z",      c.z,      fits::at(il));
        otbl.update_elements("zmin",   c.zmin,   fits::at(il));
        otbl.update_elements("zmax",   c.zmax,   fits::at(il));
        ++il;
    });

    return 0;
}
