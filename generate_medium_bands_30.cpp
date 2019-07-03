#include <vif.hpp>

using namespace vif;

int vif_main(int argc, char* argv[]) {
    float lmin = 400, lmax = 1000;
    float dl = 30;
    uint_t nband = ceil((lmax - lmin)/dl);

    file::mkdir("mband30");

    for (uint_t i : range(nband)) {
        vec1f lam = {lmin + i*dl, lmin + (i+1)*dl};
        vec1f res = {1.0, 1.0};
        lam /= 1e3;
        res /= integrate(lam, res);
        fits::write_table("mband30/"+to_string(long(lmin + i*dl))+".fits", ftable(lam, res));
    }

    return 0;
}
