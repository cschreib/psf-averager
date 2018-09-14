#include "psf-averager-photoz.hpp"

class egg_averager : public psf_averager {
public :

    egg_averager() : psf_averager("egg") {
        single_pass = true;
        global_progress_bar = true;
    }

    void process_cached(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {}

    void do_fit(uint_t iter, uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {}

    std::string make_cache_hash() override {
        return "";
    }
};

int vif_main(int argc, char* argv[]) {
    // External data
    std::string share_dir = "/home/cschreib/code/egg-analytic/share/";
    std::string filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
    std::string sed_lib = "/home/cschreib/code/egg-analytic/share/opt_lib_fastpp_hd_noigm.fits";
    std::string sed_imf = "chabrier";
    std::string psf_file  = "/home/cschreib/code/euclid_psf/psf-averager/psf-mono.fits";

    // Survey definition
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";
    double dz = 0.01;
    uint_t seds_step = 5;
    uint_t iz = 5;
    uint_t nthread = 0;
    bool write_individuals = false;
    bool write_averages = true;
    std::string cache_id;
    vec1d zb = {0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500};
    vec1s filters;

    read_args(argc, argv, arg_list(maglim, selection_band, dz, seds_step, iz,
        nthread, write_individuals, write_averages, cache_id, share_dir, filter_db,
        sed_lib, sed_imf, psf_file, zb, filters));

    egg_averager pavg;
    pavg.write_cache = false;
    pavg.use_cache = false;

    // Setup survey
    egg::generator_options opts;
    opts.share_dir = share_dir;
    opts.filter_db = filter_db;
    opts.sed_lib = sed_lib;
    opts.sed_lib_imf = sed_imf;
    opts.filter_flambda = true;
    opts.filter_photons = true;
    opts.selection_band = selection_band;
    opts.maglim = maglim;
    opts.logmass_steps = 50;
    opts.bt_steps = 5;
    opts.logmass_max = 12.0;
    opts.seds_step = seds_step;
    opts.nthread = nthread;
    opts.filters = filters;
    pavg.initialize(opts);

    // Setup mock
    mock_options mopts;
    mopts.nmc = 1;
    mopts.min_mag_err = 0.0;
    mopts.zb = zb;
    mopts.dz = dz;
    mopts.no_noise = true;
    mopts.psf_file = psf_file;
    mopts.write_individuals = write_individuals;
    mopts.write_averages = write_averages;
    mopts.force_cache_id = cache_id;
    mopts.depths = replicate(29.0, filters.size());
    pavg.configure_mock(mopts);

    if (!pavg.average_redshift_bin(iz)) return 1;

    return 0;
}
