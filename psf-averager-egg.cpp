#include "psf-averager-photoz.hpp"

class egg_averager : public psf_averager {
public :

    egg_averager() : psf_averager("egg") {
        single_pass = true;
        global_progress_bar = true;
    }

    void process_cached(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {}

    void do_fit(uint_t id_mass, uint_t id_type, uint_t id_disk, uint_t id_bulge,
        uint_t id_bt, double tngal, const vec1d& fdisk, const vec1d& fbulge) override {}

    void set_priors(const vec1d& fdisk, const vec1d& fbulge) override {}

    std::string make_cache_hash() override {
        return "";
    }
};

int phypp_main(int argc, char* argv[]) {
    // Survey definition
    double maglim = 24.5;
    std::string selection_band = "euclid-vis";
    double dz = 0.01;
    uint_t seds_step = 5;
    uint_t iz = 5;

    read_args(argc, argv, arg_list(maglim, selection_band, dz, seds_step, iz));

    egg_averager pavg;
    pavg.write_cache = false;
    pavg.use_cache = false;

    // Setup survey
    egg::generator_options opts;
    opts.share_dir = "/home/cschreib/code/egg/share/";
    opts.filter_db = "/home/cschreib/code/euclid_psf/psf-averager/filters.dat";
    opts.filter_flambda = true;
    opts.filter_photons = true;
    opts.selection_band = selection_band;
    opts.maglim = maglim;
    opts.logmass_steps = 50;
    opts.bt_steps = 5;
    opts.logmass_max = 12.0;
    opts.seds_step = seds_step;
    pavg.initialize(opts);

    // Setup mock
    mock_options mopts;
    mopts.nmc = 1;
    mopts.min_mag_err = 0.0;
    mopts.dz = dz;
    mopts.no_noise = true;
    mopts.psf_dir = "/home/cschreib/work_psf/psf-library/";
    pavg.configure_mock(mopts);

    if (!pavg.average_redshift_bin(iz)) return 1;

    return 0;
}
