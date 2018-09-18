#ifndef COMMON_INCLUDED
#define COMMON_INCLUDED

#include <vif.hpp>

using namespace vif;

vec1d resample_sed(const vec1d& tsed, const vec1d& tlam, const vec1d& blam) {
    double dl = blam[1] - blam[0];
    uint_t npt = blam.size();
    vec1d fobs(npt);
    for (uint_t l : range(npt)) {
        if (blam.safe[l]-dl/2.0 > tlam.front() && blam.safe[l]+dl/2.0 < tlam.back()) {
            fobs.safe[l] = integrate(tlam, tsed, blam.safe[l]-dl/2.0, blam.safe[l]+dl/2.0)/dl;
        }
    }

    return fobs;
}

#endif
