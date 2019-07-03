#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

#include "tree.hpp"

int vif_main(int argc, char* argv[]) {
    vec1s bands_mb;
    vec1d lambda_mb;
    vec1f z;

    vec2f colors;
    vec2f fluxes;
    vec2f flux_mb;

    auto seed = make_seed(42);

    // Define the colors to use
    vec1s cbands = {"sdss-g","sdss-r","sdss-i","sdss-z","vista-Y","vista-J","vista-H"};
    vec1f cdepths = {27.2,   27.2,    27.0,    26.0,    25.7,     25.4,     24.9};
    cdepths = e10(0.4*(23.9-cdepths))/5.0;

    double mbdepth = 26.0;
    mbdepth = e10(0.4*(23.9-mbdepth))/5.0;

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

    vec1d norm;

    vec1b selected;

    {
        vec1s bands;
        vec1f lambda;
        vec2f flux;
        fits::read_table("egg-catalog.fits", ftable(flux, bands, lambda, z));

        // Read in broad band fluxes
        norm.resize(z.dims[0]);
        fluxes.resize(z.dims[0], cbands.size());
        for (uint_t b : range(cbands)) {
            fluxes(_,b) = flux(_,where_first(bands == cbands[b]));

            // Add noise
            fluxes(_,b) += randomn(seed, z.size())*cdepths[b];
        }

        uint_t bsel = where_first(cbands == "sdss-i");
        double fsel = e10(0.4*(23.9 - 24.5));
        selected = fluxes(_,bsel) > fsel;

        // Normalize to unit flux at reddest band
        for (uint_t i : range(z)) {
            norm[i] = fluxes(i,cbands.size()-1);
            fluxes(i,_) /= norm[i];
        }

        // Compute colors
        cname.resize(cbands.size()-1);
        colors.resize(z.dims[0], cname.size());
        for (uint_t b : range(cbands.size()-1)) {
            cname[b] = cbands[b]+"/"+cbands[b+1];
            colors(_,b) = -2.5*log10(fluxes(_,b)/fluxes(_,b+1));
        }

        vec1u idmb = where(begins_with(bands, "mb30"));
        bands_mb = bands[idmb];
        lambda_mb = lambda[idmb];
        flux_mb = flux(_,idmb);
    }

    uint_t nb = bands_mb.size();
    uint_t ncol = cname.size();

    // Train
    note(count(selected), " galaxies used for training");
    note("design bins");

    vec<1,vec2f> bins(ncol);
    vec<1,vec1f> bcenter(ncol);
    vec<1,vec1f> bwidth(ncol);
    for (uint_t i : range(ncol)) {
        // Color is clamped to 3 to avoid long tails toward +inf for dropouts
        vec1f c = colors(_,i);
        c = min(c, 3.0);

        // Select good data points
        vec1u ids = where(is_finite(c) && selected);

        vec2f b;

        if (smart_bin) {
            // Sort data
            ids = ids[sort(c[ids])];

            // Define number of bins
            uint_t nm = nmag;
            if (cname[i] == "sdss-r/sdss-i" || cname[i] == "sdss-i/sdss-z") nm = nmag_riz;

            // Define bins with constant number density (refine where there are lots of galaxies)
            double ngal_step = ids.size()/float(nm+1);
            vec1d edges = c[ids[max(vec1u{round(ngal_step*indgen(nm+1))}, ids.size()-1)]];

            print(cname[i], ": ", edges.size()-1, " : ", edges);
            b = make_bins_from_edges(edges);
        } else {
            // Get bounds
            double mi = min(c[ids]);
            double ma = max(c[ids]);

            // Define width of bins
            double dm = dmag;
            if (cname[i] == "sdss-r/sdss-i" || cname[i] == "sdss-i/sdss-z") dm = dmag_riz;

            b = make_bins_from_edges(rgen_step(mi, ma, dm));
            print(cname[i], ": ", mi, ", ", ma, ", ", b.dims[1]);
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

    vec2u id_true;
    if (!write_associations.empty()) {
        id_true.resize(z.size(), bins.size());
    }

    // Grab fluxes and store
    auto pg = progress_start(z.size());
    for (uint_t i : range(z)) {
        if (!selected[i]) {
            progress(pg);
            continue;
        }

        // Read MB flux
        vec1f f = flux_mb(i,_);

        // Add noise
        f += randomn(seed, f.size())*mbdepth;

        // Normalize MB flux with same factor as BB
        f /= norm[i];

        // Remove invalid data (z=0 galaxies typically)
        if (!is_finite(f[0])) {
            progress(pg);
            continue;
        }

        // Find leaf of the tree with this color
        vec1u ik = mf.bin_key(colors(i,_));
        if (!write_associations.empty()) {
            id_true(i,_) = ik;
        }

        auto& c = mf.insert_binned(ik);

        // Average stuff using number density as weighting
        double w = 1.0;
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


    // Save associations
    if (!write_associations.empty()) {
        fits::write_table(write_associations, ftable(id_true));
    }

    return 0;
}
