
template<typename TypeX, typename TypeN, typename F>
void rebin_generic(const vec<1,TypeX>& x, const vec<1,TypeN>& xn, F&& func) {
    uint_t n = x.size();
    uint_t nn = xn.size();

    // Build bins
    vec<1,TypeX> xl(n);
    vec<1,TypeX> xu(n);
    vec<1,TypeX> dx(n);
    for (uint_t i : range(x)) {
        xl.safe[i] = (i == 0   ? x.safe[0]   - 0.5*(x.safe[1]-x.safe[0])     : 0.5*(x.safe[i-1] + x.safe[i]));
        xu.safe[i] = (i == n-1 ? x.safe[n-1] + 0.5*(x.safe[n-1]-x.safe[n-2]) : 0.5*(x.safe[i+1] + x.safe[i]));
        dx.safe[i] = xu.safe[i] - xl.safe[i];
    }

    vec<1,TypeN> xedges(nn+1);
    vec1u        iedges(nn+1);

    // Locate integration bounds in integration axis
    uint_t i0 = npos;
    for (uint_t i : range(xedges)) {
        if (i < nn) {
            xedges.safe[i] = (i == 0 ?
                xn.safe[0] - 0.5*(xn.safe[1]-xn.safe[0]) :
                0.5*(xn.safe[i-1]+xn.safe[i]));
        } else {
            xedges.safe[nn] = xn.safe[nn-1] + 0.5*(xn.safe[nn-1]-xn.safe[nn-2]);
        }

        if (i == 0) {
            vif_check(xedges[0] >= xl[0], "requested grid would extrapolate data");
            i0 = lower_bound(xu, xedges[0]);
        }

        while (i0 != xu.size() && xu.safe[i0] < xedges.safe[i]) {
            ++i0;
        }

        vif_check(i0 != xu.size(), "requested grid would extrapolate data");

        iedges.safe[i] = i0;
    }

    // Integrate
    func(xl, xu, dx, xedges, iedges);
}

template<typename TypeY, typename TypeX, typename TypeN>
vec<1,TypeY> rebin_cst(const vec<1,TypeY>& y, const vec<1,TypeX>& x, const vec<1,TypeN>& xn) {
    vec<1,TypeY> data(xn.size());

    rebin_generic(x, xn, [&](
        const vec<1,TypeX>& xl, const vec<1,TypeX>& xu, const vec<1,TypeX>& dx,
        const vec<1,TypeN>& xedges, const vec1u& iedges) {

        for (uint_t i : range(xn)) {
            uint_t iprev = iedges.safe[i];
            uint_t inext = iedges.safe[i+1];
            if (iprev == inext) {
                // Output bin is fully covered by one input bin
                data.safe[i] = y.safe[iprev];
            } else {
                // Add input bins fully covered by output bin
                meta::rtype_t<TypeY> v = 0.0;
                for (uint_t j : range(iprev+1, inext)) {
                    v += y.safe[j]*dx.safe[j];
                }
                // Add left edge
                v += y.safe[iprev]*(xu.safe[iprev]   - xedges.safe[i]);
                // Add right edge
                v += y.safe[inext]*(xedges.safe[i+1] - xl.safe[inext]);
                // Average
                v /= (xedges.safe[i+1] - xedges.safe[i]);
                data.safe[i] = v;
            }
        }
    });


    return data;
}

template<typename TypeY, typename TypeX, typename TypeN>
vec<1,TypeY> rebin_trapz(const vec<1,TypeY>& y, const vec<1,TypeX>& x, const vec<1,TypeN>& xn) {
    vec<1,TypeY> data(xn.size());

    rebin_generic(x, xn, [&](
        const vec<1,TypeX>& xl, const vec<1,TypeX>& xu, const vec<1,TypeX>& dx,
        const vec<1,TypeN>& xedges, const vec1u& iedges) {

        for (uint_t i : range(xn)) {
            // TODO: can optimize this using iedges etc
            data.safe[i] = integrate(x, y, xedges.safe[i], xedges.safe[i+1])/(xedges.safe[i+1] - xedges.safe[i]);
        }
    });

    return data;
}

template<typename TypeY, typename TypeX, typename TypeN>
vec<1,TypeY> rebin_mcspline(const vec<1,TypeY>& y, const vec<1,TypeX>& x, const vec<1,TypeN>& xn) {
    vec<1,TypeY> data(xn.size());

    rebin_generic(x, xn, [&](
        const vec<1,TypeX>& xl, const vec<1,TypeX>& xu, const vec<1,TypeX>& dx,
        const vec<1,TypeN>& xedges, const vec1u& iedges) {

        uint_t gn = x.size();
        vec2d alpha(5,gn);

        matrix::decompose_lu dcc;
        double old_dx = 0.0;

        auto get_coefs = [&](uint_t i0, uint_t i1, uint_t ik0, uint_t ik1) {
            // Build spline
            uint_t n = i1 - i0;
            uint_t nk = ik1 - ik0;
            uint_t kb = ik0 - i0;
            uint_t m = n+1;

            // Try to be smart and catch cases where dx is constant and the same as
            // in the previous call, so we can reuse the decomposed matrix
            double this_dx = dx.safe[i0];
            bool constant_dx = true;
            for (uint_t i : range(1, n)) {
                if (abs(this_dx - dx.safe[i0+i]) > 1e-6) {
                    constant_dx = false;
                    break;
                }
            }

            vec1d s;

            if (constant_dx && abs(this_dx - old_dx) < 1e-6 && dcc.size() == 2*m) {
                // Re-use decomposition!
                vec1d b(2*m);
                for (uint_t i : range(1, n)) {
                    b.safe[3+i] = 60*(y.safe[i0+i] - y.safe[i0+i-1]);
                }
                s = dcc.solve(b);
            } else {
                // Re-compute decomposition
                // Adapted from:
                // "A Spline Interpolation Technique that Preserves Mass Budgets"
                // E.J.M. Delhez, 2003, Applied Mathematics Letters, 16, 17-26

                matrix::mat2d a(2*m, 2*m);
                vec1d b(2*m);

                if (constant_dx) {
                    double tdx = dx.safe[0];
                    double tdx2 = sqr(tdx);

                    // Natural boundary conditions
                    a.safe(0,m) = 1;

                    a.safe(1,m+m-1) = 1;

                    a.safe(2,0)   = +3;
                    a.safe(2,1)   = -3;
                    a.safe(2,m+0) = +2*tdx;
                    a.safe(2,m+1) = +1*tdx;

                    a.safe(3,m-2)   = +3;
                    a.safe(3,m-1)   = -3;
                    a.safe(3,m+m-2) = +1*tdx;
                    a.safe(3,m+m-1) = +2*tdx;

                    // Continuity constraints
                    for (uint_t i : range(1, n)) {
                        a.safe(3+i,i-1)   = +9*tdx;
                        a.safe(3+i,i)     = +21*2*tdx;
                        a.safe(3+i,i+1)   = +9*tdx;

                        a.safe(3+i,m+i-1) = +2*tdx2;
                        a.safe(3+i,m+i+1) = -2*tdx2;

                        b.safe[3+i] = 60*(y.safe[i0+i] - y.safe[i0+i-1]);
                    }

                    for (uint_t i : range(1, n)) {
                        a.safe(1+m+i,i-1)   = +6*tdx;
                        a.safe(1+m+i,i+1)   = -6*tdx;

                        a.safe(1+m+i,m+i-1) = +2*tdx2;
                        a.safe(1+m+i,m+i)   = +4*2*tdx2;
                        a.safe(1+m+i,m+i+1) = +2*tdx2;
                    }
                } else {
                    // Natural boundary conditions
                    a.safe(0,m) = 1;

                    a.safe(1,m+m-1) = 1;

                    a.safe(2,0)   = +3;
                    a.safe(2,1)   = -3;
                    a.safe(2,m+0) = +2*dx.safe[i0];
                    a.safe(2,m+1) = +1*dx.safe[i0];

                    a.safe(3,m-2)   = +3;
                    a.safe(3,m-1)   = -3;
                    a.safe(3,m+m-2) = +1*dx.safe[i0+n-1];
                    a.safe(3,m+m-1) = +2*dx.safe[i0+n-1];

                    // Continuity constraints
                    for (uint_t i : range(1, n)) {
                        a.safe(3+i,i-1)   = +9*dx.safe[i0+i-1];
                        a.safe(3+i,i)     = +21*(dx.safe[i0+i-1] + dx.safe[i0+i]);
                        a.safe(3+i,i+1)   = +9*dx.safe[i0+i];

                        a.safe(3+i,m+i-1) = +2*sqr(dx.safe[i0+i-1]);
                        a.safe(3+i,m+i)   = +3*(sqr(dx.safe[i0+i-1]) - sqr(dx.safe[i0+i]));
                        a.safe(3+i,m+i+1) = -2*sqr(dx.safe[i0+i]);

                        b.safe[3+i] = 60*(y.safe[i0+i] - y.safe[i0+i-1]);
                    }

                    for (uint_t i : range(1, n)) {
                        a.safe(1+m+i,i-1)   = +6*dx.safe[i0+i-1];
                        a.safe(1+m+i,i)     = -6*(dx.safe[i0+i-1] - dx.safe[i0+i]);
                        a.safe(1+m+i,i+1)   = -6*dx.safe[i0+i];

                        a.safe(1+m+i,m+i-1) = +2*sqr(dx.safe[i0+i-1]);
                        a.safe(1+m+i,m+i)   = +4*(sqr(dx.safe[i0+i-1]) + sqr(dx.safe[i0+i]));
                        a.safe(1+m+i,m+i+1) = +2*sqr(dx.safe[i0+i]);
                    }
                }

                // Solve
                matrix::decompose_lu dc;
                vif_check(dc.decompose(a), "could not do LU decomposition of problem");
                s = dc.solve(b);

                if (constant_dx) {
                    // Cache the decomposition so we can re-use it next time
                    dcc = std::move(dc);
                    old_dx = this_dx;
                }
            }

            // Evaluate spline coefficients
            for (uint_t i : range(nk)) {
                uint_t iy1 = kb+i;
                uint_t iy2 = kb+i+m;

                double tdx = dx.safe[ik0+i];

                alpha.safe(0,ik0+i) = y.safe[ik0+i] - (1.0/60.0)*((21*s.safe[iy1] + 9*s.safe[iy1+1]) +
                                                      (3*s.safe[iy2] - 2*s.safe[iy2+1])*tdx)*tdx;
                alpha.safe(1,ik0+i) = s.safe[iy1]*tdx;
                alpha.safe(2,ik0+i) =  (1.0/2.0)*s.safe[iy2]*sqr(tdx);
                alpha.safe(3,ik0+i) = -(1.0/3.0)*(3.0*(s.safe[iy1] - s.safe[iy1+1]) +
                                                 (2*s.safe[iy2] + s.safe[iy2+1])*tdx)*tdx;
                alpha.safe(4,ik0+i) =  (1.0/4.0)*(2.0*(s.safe[iy1] - s.safe[iy1+1]) +
                                                 (s.safe[iy2] + s.safe[iy2+1])*tdx)*tdx;
            }
        };

        // Get spline coefs in chunks:
        //
        //  [-------[======]-------]
        //   nmargin nblock nmargin
        //
        // Compute spline over nmargin+nblock+nmargin, but only keep nblock.

        // Optimized setup for large number of data points
        uint_t nmargin = 20; // lowest value that gives exact results (but why?)
        uint_t nblock = 10; // profiled for speed

        // Optimized for small number of data points
        if (gn <= nblock+2*nmargin) {
            nmargin = 0;
            nblock = gn;
        } else if (gn <= 2*(nblock+nmargin)) {
            nblock = gn/2;
        }

        uint_t nstep = ceil(gn/float(nblock));
        for (uint_t i : range(nstep)) {
            uint_t ik0 = i*nblock;
            uint_t ik1 = min(gn, (i+1)*nblock);
            uint_t i0 = (ik0 > nmargin ? ik0 - nmargin : 0);
            uint_t i1 = (gn-ik1 >= nmargin ? ik1 + nmargin : gn);
            get_coefs(i0, i1, ik0, ik1);
        }

        // Integrate
        auto integrate_spline = [&](uint_t i, double txl, double txu) {
            txl = (txl - xl.safe[i])/dx.safe[i];
            txu = (txu - xl.safe[i])/dx.safe[i];

            double ret = 0.0;
            for (uint_t j : range(5)) {
                ret += alpha.safe(j,i)*(pow(txu, j+1) - pow(txl, j+1))/(j+1.0);
            }

            return dx.safe[i]*ret;
        };

        for (uint_t i : range(xn)) {
            uint_t iprev = iedges.safe[i];
            uint_t inext = iedges.safe[i+1];
            if (iprev == inext) {
                // Output bin is fully covered by one input bin
                data.safe[i] = integrate_spline(iprev, xedges.safe[i], xedges.safe[i+1]);
            } else {
                // Add input bins fully covered by output bin
                meta::rtype_t<TypeY> v = 0.0;
                for (uint_t j : range(iprev+1, inext)) {
                    v += y.safe[j]*dx.safe[j];
                }
                // Add left edge
                v += integrate_spline(iprev, xedges.safe[i], xu.safe[iprev]);
                // Add right edge
                v += integrate_spline(inext, xl.safe[inext], xedges.safe[i+1]);
                data.safe[i] = v;
            }

            // Average
            data.safe[i] /= (xedges.safe[i+1] - xedges.safe[i]);
        }
    });

    return data;
}

// template<typename TypeY, typename TypeX, typename TypeN>
// vec<1,TypeY> rebin_mcspline3(const vec<1,TypeY>& y, const vec<1,TypeX>& x, const vec<1,TypeN>& xn) {
//     vec<1,TypeY> data(xn.size());

//     rebin_generic(x, xn, [&](
//         const vec<1,TypeX>& xl, const vec<1,TypeX>& xu, const vec<1,TypeX>& dx,
//         const vec<1,TypeN>& xedges, const vec1u& iedges) {

//         uint_t gn = x.size();
//         vec2d s(gn);
//         vec2d m(gn);

//         // Solve for the 'm' coefficients
//         {
//             vec2d a(gn,gn);

//             // First derivative boundary conditions
//             a(0,1)       = 1.0;
//             a(gn-1,gn-1) = 1.0;

//             // Guesstimate derivatives from finite difference of integral data (wrong!)
//             double yp0 = (y[1]    - y[0])/dx[0];
//             double yp1 = (y[gn-1] - y[gn-2])/dx[gn-2];

//             b[0]    = yp0;
//             b[gn-1] = yp1;

//             // Continuity constraints
//             for (uint_t i : range(1, gn-1)) {
//                 a(i,i-1) = 1.0;
//                 a(i,i+0) = 10.0;
//                 a(i,i+1) = 1.0;

//                 b[i] = (12.0/sqr(dx[i]))*(y[i+1] - y[i]);
//             }

//             // Solve
//             matrix::decompose_lu dc;
//             vif_check(dc.decompose(a), "could not do LU decomposition of problem");
//             m = dc.solve(b);
//         }

//         // Solve for the 's' coefficients
//         {
//             s(0) =
//         }


//         // Integrate
//         auto integrate_spline = [&](uint_t i, double txl, double txu) {
//             double txll = (txl - xl.safe[i])/dx.safe[i];
//             double txlu = (txu - xl.safe[i])/dx.safe[i];
//             double txul = (txl - xu.safe[i])/dx.safe[i];
//             double txuu = (txu - xu.safe[i])/dx.safe[i];

//             double ret = 0.0;
//             ret += s.safe[i]*(  +(1.0/1.0)*(pow(txuu, 3) - pow(txul, 3))
//                                 +(1.0/2.0)*(pow(txuu, 4) - pow(txul, 4)));
//             ret += s.safe[i+1]*(+(1.0/1.0)*(pow(txlu, 3) - pow(txll, 3))
//                                 -(1.0/2.0)*(pow(txlu, 4) - pow(txll, 4)));
//             ret += m.safe[i]*(  +(1.0/3.0)*(pow(txuu, 3) - pow(txul, 3))
//                                 +(1.0/4.0)*(pow(txuu, 4) - pow(txul, 4)));
//             ret += m.safe[i+1]*(-(1.0/3.0)*(pow(txlu, 3) - pow(txll, 3))
//                                 +(1.0/4.0)*(pow(txlu, 4) - pow(txll, 4)));

//             return dx.safe[i]*ret;
//         };

//         for (uint_t i : range(xn)) {
//             uint_t iprev = iedges.safe[i];
//             uint_t inext = iedges.safe[i+1];
//             if (iprev == inext) {
//                 // Output bin is fully covered by one input bin
//                 data.safe[i] = integrate_spline(iprev, xedges.safe[i], xedges.safe[i+1]);
//             } else {
//                 // Add input bins fully covered by output bin
//                 meta::rtype_t<TypeY> v = 0.0;
//                 for (uint_t j : range(iprev+1, inext)) {
//                     v += y.safe[j]*dx.safe[j];
//                 }
//                 // Add left edge
//                 v += integrate_spline(iprev, xedges.safe[i], xu.safe[iprev]);
//                 // Add right edge
//                 v += integrate_spline(inext, xl.safe[inext], xedges.safe[i+1]);
//                 data.safe[i] = v;
//             }

//             // Average
//             data.safe[i] /= (xedges.safe[i+1] - xedges.safe[i]);
//         }
//     });

//     return data;
// }
