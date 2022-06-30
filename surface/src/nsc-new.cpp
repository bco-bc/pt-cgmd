/*
 *
 *  program NSC
 *  version 1.0  (April 1994)
 *
 *  Author: Frank Eisenhaber
 *
 *  For user notes see file nsc.h !!
 *
 *  Copyright Notice:
 *  All rights reserved, whether the whole or part of the program is
 *  concerned. The software may not be used without specific, prior 
 *  written permission of the author. 
 *
 *  An academic licence agreement for the package ASC/GM or its parts
 *  is granted if you make the following commitments:
 *  1) In using this software, the user will respect the interests of 
 *     the author.
 *  2) The use of the software in commercial activities is not allowed 
 *     without a prior written commercial licence agreement. The program
 *     will not be used in classified research.
 *  3) Other interested research groups will be redirected
 *     to the author. The user will not redistribute the code outside
 *     his immediate research group.
 *  4) The copyright messages will not be modified or suppressed.
 *  5) The reference given below will be cited in any publication
 *     of scientific results based in part or completely on use of the
 *     program.
 *  6) Bugs will be reported to the author.
 *
 *  Permission to use, copy, and modify this software and
 *  its documentation is hereby granted without fee for 
 *  academic use, provided
 *  that the above copyright notices and this permission notice appear in
 *  all copies of the software and related documentation.
 *
 *  THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF
 *  ANY KIND,
 *  EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
 *  WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
 *
 *  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 *  ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
 *  OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 *  WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF
 *  LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
 *  OF THIS SOFTWARE.
 *
 * 
 *  contact address :    European Molecular Biology Laboratory
 *                       Biocomputing Unit
 *                       Meyerhofstr. 1
 *                       Postfach 10.2209
 *                       D-69012 Heidelberg
 *                       Federal Republic of Germany
 *
 *
 *
 *  E-mail : IN%"EISENHABER@EMBL-Heidelberg.DE"
 *  Please send your contact address to get information on updates and
 *  new features. Questions will be answered as soon as possible.
 *
 *
 *  references :
 *  1.F.Eisenhaber, P.Lijnzaad, P.Argos, M.Scharf
 *    "The Double Cubic Lattice Method: Efficient Approaches to
 *    Numerical Integration of Surface Area and Volume and to Dot
 *    Surface Contouring of Molecular Assemblies"
 *    Journal of Computational Chemistry (1994) submitted
 *  2.F.Eisenhaber, P.Argos
 *    "Improved Strategy in Analytic Surface Calculation for Molecular
 *    Systems: Handling of Singularities and Computational Efficiency"
 *    Journal of Computational Chemistry (1993) v.14, N11, pp-1272-1280
 *
 */

#include "simploce/surface/nsc-new.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include "boost/lexical_cast.hpp"

static const int UNSP_ICO_DOD = 9;
static const int UNSP_ICO_ARC = 10;
static const double FOURPI    = 4.0 * M_PI;
static const double DP_TOL    = 0.001;

using double_ptr_t = double *;

static std::vector<double> XPUNSP;
static double del_cube_{0.0};
static std::vector<int> ICO_WK;
static std::vector<int>::iterator ICO_ITER;
static int N_DOT, ISO_CUBE, last_n_dot=0, last_densit=0, last_unsp=0;
static int LAST_CUBUS = 0;

/* routines for dot distributions on the surface of the unit sphere */
static double RG{0.0};
static double RH{0.0};


static double torad_(double A) {
    return A * 0.017453293;
}


static double asin_(double f) {
    if ((std::fabs(f) < 1.00))
        return (std::asin(f));
    if ((std::fabs(f) - 1.00) < DP_TOL) {
        return (M_PI_2);
    } else {
        throw std::domain_error("NSC::asin_: invalid argument: " + boost::lexical_cast<std::string>(f));
    }
}


static std::vector<double> icosaederVertices_() {
    RH = sqrt(1. - 2. * cos(torad_(72.))) / (1. - cos(torad_(72.)));
    RG = cos(torad_(72.)) / (1. - cos(torad_(72.)));
    /* icosaeder vertices */
    std::vector<double> xus(3 * N_DOT, 0.0);
    xus[ 0] = 0.;
    xus[ 1] = 0.;
    xus[ 2] = 1.;
    xus[ 3] = RH * cos(torad_(72.));
    xus[ 4] = RH * sin(torad_(72.));
    xus[ 5] = RG;
    xus[ 6] = RH * cos(torad_(144.));
    xus[ 7] = RH * sin(torad_(144.));
    xus[ 8] = RG;
    xus[ 9] = RH * cos(torad_(216.));
    xus[10] = RH * sin(torad_(216.));
    xus[11] = RG;
    xus[12] = RH * cos(torad_(288.));
    xus[13] = RH * sin(torad_(288.));
    xus[14] = RG;
    xus[15] = RH;
    xus[16] = 0;
    xus[17] = RG;
    xus[18] = RH * cos(torad_(36.));
    xus[19] = RH * sin(torad_(36.));
    xus[20] = -RG;
    xus[21] = RH * cos(torad_(108.));
    xus[22] = RH * sin(torad_(108.));
    xus[23] = -RG;
    xus[24] = -RH;
    xus[25] = 0;
    xus[26] = -RG;
    xus[27] = RH * cos(torad_(252.));
    xus[28] = RH * sin(torad_(252.));
    xus[29] = -RG;
    xus[30] = RH * cos(torad_(324.));
    xus[31] = RH * sin(torad_(324.));
    xus[32] = -RG;
    xus[33] = 0.;
    xus[34] = 0.;
    xus[35] = -1.;
    return std::move(xus);
}


static void divarc_(double x1, double y1, double z1, double x2, double y2, double z2, int div1, int div2,
                    double *xr, double *yr, double *zr) {

    auto xd = y1 * z2 - y2 * z1;
    auto yd = z1 * x2 - z2 * x1;
    auto zd = x1 * y2 - x2 * y1;
    auto dd = sqrt(xd * xd + yd * yd + zd * zd);
    if (dd < DP_TOL) {
      throw std::domain_error("divarc_: rotation axis of length: " + boost::lexical_cast<std::string>(dd));
    }
    auto d1 = x1*x1+y1*y1+z1*z1;
    if (d1 < 0.5) {
      throw std::domain_error("divarc_: vector 1 of sq.length:" + boost::lexical_cast<std::string>(d1));
    }
    auto d2 = x2*x2+y2*y2+z2*z2;
    if (d2 < 0.5) {
        throw std::domain_error("divarc_: vector 2 of sq.length: " + boost::lexical_cast<std::string>(d2));
    }

    auto phi = asin_(dd / sqrt(d1 * d2));
    phi = phi*((double)div1)/((double)div2);
    auto sphi = sin(phi);
    auto cphi = cos(phi);
    auto s  = (x1*xd+y1*yd+z1*zd)/dd;

    auto x = xd*s*(1.-cphi)/dd + x1 * cphi + (yd*z1-y1*zd)*sphi/dd;
    auto y = yd*s*(1.-cphi)/dd + y1 * cphi + (zd*x1-z1*xd)*sphi/dd;
    auto z = zd*s*(1.-cphi)/dd + z1 * cphi + (xd*y1-x1*yd)*sphi/dd;
    dd = sqrt(x*x+y*y+z*z);
    *xr = x/dd;
    *yr = y/dd;
    *zr = z/dd;
}

int ico_dot_arc(int densit) {
    /* densit...required dots per unit sphere */
    /* dot distribution on a unit sphere based on an icosaeder */
    /* great circle average refining of icosahedral face       */

    int i, j, k, tl, tl2, tn, tess;
    double a, d, x, y, z, x2, y2, z2, x3, y3, z3;
    double xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
    xjk, yjk, zjk, xkj, ykj, zkj;

    /* calculate tessalation level */
    a = sqrt((((double) densit)-2.)/10.);
    tess = (int) ceil(a);
    N_DOT = 10 * tess * tess + 2;
    if (N_DOT < densit) {
        //ERROR("ico_dot_arc: error in formula for tessalation level (%d->%d, %d)",
        //  tess, n_dot, densit);
        throw std::domain_error("ico_dot_arc: error in formula for tessalation level.");
    }

    auto xus = icosaederVertices_();
    XPUNSP = xus;

    if (tess > 1) {
        tn = 12;
        a = RH * RH * 2. * (1. - cos(torad_(72.)));
        /* calculate tessalation of icosaeder edges */
        for (i=0; i<11; i++) {
            for (j=i+1; j<12; j++) {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL) continue;

                for (tl=1; tl<tess; tl++) {
                    if (tn >= N_DOT) {
                        //ERROR("ico_dot: tn exceeds dimension of xus");
                        throw std::domain_error("ico_dot: tn exceeds dimension of xus");
                    }
                    divarc_(xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i],
                            xus[3 * j], xus[1 + 3 * j], xus[2 + 3 * j],
                            tl, tess, &xus[3 * tn], &xus[1 + 3 * tn], &xus[2 + 3 * tn]);
                    tn++;
                }
            }
        }
        /* calculate tessalation of icosaeder faces */
        for (i=0; i<10; i++) {
            for (j=i+1; j<11; j++) {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL) continue;

                for (k=j+1; k<12; k++) {
                    x = xus[3*i]-xus[3*k];
                    y = xus[1+3*i]-xus[1+3*k]; z = xus[2+3*i]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if (fabs(a-d) > DP_TOL) continue;

                    x = xus[3*j]-xus[3*k];
                    y = xus[1+3*j]-xus[1+3*k]; z = xus[2+3*j]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if (fabs(a-d) > DP_TOL) continue;

                    for (tl=1; tl<tess-1; tl++) {
                        divarc_(xus[3 * j], xus[1 + 3 * j], xus[2 + 3 * j],
                            xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i],tl, tess, &xji, &yji, &zji);
                        divarc_(xus[3 * k], xus[1 + 3 * k], xus[2 + 3 * k],
                                xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i],
                                tl, tess, &xki, &yki, &zki);

                        for (tl2=1; tl2<tess-tl; tl2++) {
                            divarc_(xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i],
                                    xus[3 * j], xus[1 + 3 * j], xus[2 + 3 * j],
                                tl2, tess, &xij, &yij, &zij);
                            divarc_(xus[3 * k], xus[1 + 3 * k], xus[2 + 3 * k],
                                    xus[3 * j], xus[1 + 3 * j], xus[2 + 3 * j],
                                tl2, tess, &xkj, &ykj, &zkj);
                            divarc_(xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i],
                                xus[3 * k], xus[1 + 3 * k], xus[2 + 3 * k],
                                tess - tl - tl2, tess, &xik, &yik, &zik);
                            divarc_(xus[3 * j], xus[1 + 3 * j], xus[2 + 3 * j],
                                    xus[3 * k], xus[1 + 3 * k], xus[2 + 3 * k],
                                tess - tl - tl2, tess, &xjk, &yjk, &zjk);
                            if (tn >= N_DOT) {
                                //ERROR("ico_dot: tn exceeds dimension of xus");
                                throw std::domain_error("ico_dot: tn exceeds dimension of xus");
                            }
                            divarc_(xki, yki, zki, xji, yji, zji, tl2, tess - tl,
                                    &x, &y, &z);
                            divarc_(xkj, ykj, zkj, xij, yij, zij, tl, tess - tl2,
                                    &x2, &y2, &z2);
                            divarc_(xjk, yjk, zjk, xik, yik, zik, tl, tl + tl2,
                                    &x3, &y3, &z3);
                            x = x+x2+x3;
                            y = y+y2+y3;
                            z = z+z2+z3;
                            d = sqrt(x*x+y*y+z*z);
                            xus[3*tn] = x/d;
                            xus[1+3*tn] = y/d;
                            xus[2+3*tn] = z/d;
                            tn++;
                        }		/* cycle tl2 */
                    }		/* cycle tl */
                }		/* cycle k */
            }		/* cycle j */
        }			/* cycle i */
        if (N_DOT != tn) {
            throw std::domain_error("co_dot: n_dot " + boost::lexical_cast<std::string>(N_DOT) +
                                      " and " + boost::lexical_cast<std::string>(tn) + " differ.");
        }
    }		/* end of if (tess > 1) */
    return N_DOT;
  }		/* end of routine ico_dot_arc */


static int ico_dot_dod(int densit) {
    /* densit...required dots per unit sphere */
    /* dot distribution on a unit sphere based on an icosaeder */
    /* great circle average refining of icosahedral face       */

    int i, j, k, tl, tl2, tn, tess, j1, j2;
    double a, d, x, y, z, x2, y2, z2, x3, y3, z3, ai_d, adod;
    double xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki, xjk, yjk, zjk, xkj, ykj, zkj;

    /* calculate tesselation level */
    a = sqrt((((double) densit)-2.)/30.);
    tess = std::max((int) ceil(a), 1);
    N_DOT = 30 * tess * tess + 2;
    if (N_DOT < densit) {
        //ERROR("ico_dot_dod: error in formula for tessalation level (%d->%d, %d)", tess, n_dot, densit);
        throw std::domain_error("ico_dot_dod: error in formula for tessalation level");
    }

    //xus = (double *) CALLOC(3*n_dot, sizeof(double));
    auto xus = icosaederVertices_();
    XPUNSP = xus;

    tn=12;
    /* square of the edge of an icosaeder */
    a = RH * RH * 2. * (1. - cos(torad_(72.)));
    /* dodecaeder vertices */
    for (i=0; i<10; i++) {
        for (j=i+1; j<11; j++) {
            x = xus[3*i]-xus[3*j];
            y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
            d = x*x+y*y+z*z;
            if (fabs(a-d) > DP_TOL) continue;

            for (k=j+1; k<12; k++) {
                x = xus[3*i]-xus[3*k];
                y = xus[1+3*i]-xus[1+3*k]; z = xus[2+3*i]-xus[2+3*k];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL) continue;

                x = xus[3*j]-xus[3*k];
                y = xus[1+3*j]-xus[1+3*k]; z = xus[2+3*j]-xus[2+3*k];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL) continue;

                x = xus[  3*i]+xus[  3*j]+xus[  3*k];
                y = xus[1+3*i]+xus[1+3*j]+xus[1+3*k];
                z = xus[2+3*i]+xus[2+3*j]+xus[2+3*k];
                d = sqrt(x*x+y*y+z*z);
                xus[3*tn]=x/d; xus[1+3*tn]=y/d; xus[2+3*tn]=z/d;
                tn++;
            }
        }
    }

    if (tess > 1) {
        tn = 32;
        /* square of the edge of an dodecaeder */
        adod = 4. * (cos(torad_(108.)) - cos(torad_(120.))) / (1. - cos(torad_(120.)));
        /* square of the distance of two adjacent vertices of ico- and dodecaeder */
        ai_d = 2.*(1.-sqrt(1.-a/3.));

        /* calculate tessalation of mixed edges */
        for (i=0; i<31; i++) {
            j1 = 12; j2 = 32; a = ai_d;
            if (i>=12) {
                j1=i+1; a = adod;
            }
            for (j=j1; j<j2; j++) {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if (fabs(a-d) > DP_TOL) continue;

                for (tl=1; tl<tess; tl++) {
                    if (tn >= N_DOT) {
                        // ERROR("ico_dot: tn exceeds dimension of xus");
                        throw std::domain_error("ico_dot: tn exceeds dimension of xus");
                    }
                    divarc_(xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i], xus[3 * j], xus[1 + 3 * j], xus[2 + 3 * j],
                            tl, tess, &xus[3 * tn], &xus[1 + 3 * tn], &xus[2 + 3 * tn]);
                    tn++;
                }
            }
        }
        /* calculate tessalation of pentakisdodecahedron faces */
        for (i=0; i<12; i++) {
            for (j=12; j<31; j++) {
                x = xus[3*i]-xus[3*j];
                y = xus[1+3*i]-xus[1+3*j];
                z = xus[2+3*i]-xus[2+3*j];
                d = x*x+y*y+z*z;
                if (fabs(ai_d-d) > DP_TOL) continue;

                for (k=j+1; k<32; k++) {
                    x = xus[3*i]-xus[3*k];
                    y = xus[1+3*i]-xus[1+3*k];
                    z = xus[2+3*i]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if (fabs(ai_d-d) > DP_TOL) continue;

                    x = xus[3*j]-xus[3*k];
                    y = xus[1+3*j]-xus[1+3*k];
                    z = xus[2+3*j]-xus[2+3*k];
                    d = x*x+y*y+z*z;
                    if (fabs(adod-d) > DP_TOL) continue;

                    for (tl=1; tl<tess-1; tl++) {
                        divarc_(xus[3 * j], xus[1 + 3 * j], xus[2 + 3 * j], xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i],
                                tl, tess, &xji, &yji, &zji);
                        divarc_(xus[3 * k], xus[1 + 3 * k], xus[2 + 3 * k], xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i],
                                tl, tess, &xki, &yki, &zki);

                        for (tl2=1; tl2<tess-tl; tl2++) {
                            divarc_(xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i], xus[3 * j], xus[1 + 3 * j],
                                    xus[2 + 3 * j],
                                    tl2, tess, &xij, &yij, &zij);
                            divarc_(xus[3 * k], xus[1 + 3 * k], xus[2 + 3 * k], xus[3 * j], xus[1 + 3 * j],
                                    xus[2 + 3 * j],
                                    tl2, tess, &xkj, &ykj, &zkj);
                            divarc_(xus[3 * i], xus[1 + 3 * i], xus[2 + 3 * i], xus[3 * k], xus[1 + 3 * k],
                                    xus[2 + 3 * k],
                                    tess - tl - tl2, tess, &xik, &yik, &zik);
                            divarc_(xus[3 * j], xus[1 + 3 * j], xus[2 + 3 * j], xus[3 * k], xus[1 + 3 * k],
                                    xus[2 + 3 * k],
                                    tess - tl - tl2, tess, &xjk, &yjk, &zjk);
                            if (tn >= N_DOT) {
                                // ERROR("ico_dot: tn exceeds dimension of xus");
                                throw std::domain_error("ico_dot: tn exceeds dimension of xus");
                            }
                            divarc_(xki, yki, zki, xji, yji, zji, tl2, tess - tl, &x, &y, &z);
                            divarc_(xkj, ykj, zkj, xij, yij, zij, tl, tess - tl2, &x2, &y2, &z2);
                            divarc_(xjk, yjk, zjk, xik, yik, zik, tl, tl + tl2, &x3, &y3, &z3);
                            x = x+x2+x3;
                            y = y+y2+y3;
                            z = z+z2+z3;
                            d = sqrt(x*x+y*y+z*z);
                            xus[3*tn] = x/d;
                            xus[1+3*tn] = y/d;
                            xus[2+3*tn] = z/d;
                            tn++;
                        }		/* cycle tl2 */
                    }		/* cycle tl */
                }		/* cycle k */
            }		/* cycle j */
        }			/* cycle i */
        if (N_DOT != tn) {
            // ERROR("ico_dot: n_dot(%d) and tn(%d) differ", n_dot, tn);
            throw std::domain_error("ico_dot: n_dot (" + boost::lexical_cast<std::string>(N_DOT) + ")" +
                              " and tn (" + boost::lexical_cast<std::string>(tn) + ") differ.");
        }
    }		/* end of if (tess > 1) */
    return N_DOT;
}		/* end of routine ico_dot_dod */


static int unsp_type(int densit) {
    int i1, i2;
    i1 = 1;
    while (10*i1*i1+2 < densit) i1++;
    i2 = 1;
    while (30*i2*i2+2 < densit) i2++;
    if (10*i1*i1-2 < 30*i2*i2-2)
        return UNSP_ICO_ARC;
    else
        return UNSP_ICO_DOD;
}


int make_unsp(int densit, int mode, int * num_dot, int cubus) {
    int ndot, ico_cube_cb, i, j, k, l, ijk, tn, tl, tl2;
    double x, y, z;


    k=1;
    if (mode < 0) {
        k=0;
        mode = -mode;
    }
    if (mode == UNSP_ICO_ARC) {
        ndot = ico_dot_arc(densit);
    } else if (mode == UNSP_ICO_DOD) {
        ndot = ico_dot_dod(densit);
    } else {
        throw std::domain_error(
                "make_unsp: mode '" + boost::lexical_cast<std::string>(mode) + "' not allowed."
                );
    }

    last_n_dot = ndot;
    last_densit = densit;
    last_unsp = mode;
    *num_dot=ndot;
    if (k)
        return 0;

    /* in the following the dots of the unit sphere may be resorted */
    last_unsp = -last_unsp;

    /* determine distribution of points in elementary cubes */
    if (cubus) {
        ISO_CUBE = cubus;
        } else {
        LAST_CUBUS = 0;
            i=1;
            while (i*i*i*2 < ndot) i++;
        ISO_CUBE = std::max(i - 1, 0);
            }
            ico_cube_cb = ISO_CUBE * ISO_CUBE * ISO_CUBE;
            del_cube_= 2.0/((double)ISO_CUBE);
            //work = (int *) CALLOC(ndot, sizeof(int));
            std::vector<int> work(ndot);
            auto xus = XPUNSP;
            for (l=0; l<ndot; l++) {
                i = std::max((int) floor((1.+xus[3*l])/del_cube_), 0);
                if (i >= ISO_CUBE) i = ISO_CUBE - 1;
                j = std::max((int) floor((1.+xus[1+3*l])/del_cube_), 0);
                if (j >= ISO_CUBE) j = ISO_CUBE - 1;
                k = std::max((int) floor((1.+xus[2+3*l])/del_cube_), 0);
                if (k >= ISO_CUBE) k = ISO_CUBE - 1;
                ijk = i + j * ISO_CUBE + k * ISO_CUBE * ISO_CUBE;
                work[l] = ijk;
            }

        // ico_wk = (int *) CALLOC(2*ico_cube_cb+1, sizeof(int));
        ICO_WK.resize(2 * ico_cube_cb + 1, 0);
        //ico_wk = new int[2*ico_cube_cb+1];
        ICO_ITER = ICO_WK.begin() + ico_cube_cb;
        for (l=0; l<ndot; l++) {
            ICO_WK[work[l]]++;   /* dots per elementary cube */
        }

        /* reordering of the coordinate array in accordance with box number */
        tn=0;
        for (i=0; i < ISO_CUBE; i++) {
            for (j=0; j < ISO_CUBE; j++) {
                for (k=0; k < ISO_CUBE; k++) {
                    tl=0;
                    tl2 = tn;
                    ijk = i + ISO_CUBE * j + ISO_CUBE * ISO_CUBE * k;
                    *(ICO_ITER + ijk) = tn;
                    for (l=tl2; l<ndot; l++) {
                        if (ijk == work[l]) {
                            x = xus[3*l]; y = xus[1+3*l]; z = xus[2+3*l];
                            xus[3*l] = xus[3*tn];
                            xus[1+3*l] = xus[1+3*tn]; xus[2+3*l] = xus[2+3*tn];
                            xus[3*tn] = x; xus[1+3*tn] = y; xus[2+3*tn] = z;
                            ijk = work[l]; work[l]=work[tn]; work[tn]=ijk;
                            tn++;
                            tl++;
                        }
                    }
                    ICO_WK[ijk] = tl;
                    //*(ico_wk+ijk) = tl;
                }		/* cycle k */
            }			/* cycle j */
        }			/* cycle i */
    return 0;
}


typedef struct stwknb_ {
    double x;
    double y;
    double z;
    double dot;
  } Neighb;


int NSC(double *co, const double *radius, int nat,
        int densit, int mode,
        double *value_of_area, double **at_area,
        double *value_of_vol,
        double **lidots, int *nu_dots) {

    int iat, i, ii, iii, ix, iy, iz, ixe, ixs, iye, iys, ize, izs, i_ac;
    int jat, j, jj, jjj, jx, jy, jz;
    int distribution;
    int l;
    int maxnei, nnei, last, maxdots;
    //int_ptr_t wkdot= nullptr, wkbox= nullptr, wkat1= nullptr, wkatm= nullptr;
        //Neighb  *wknb, *ctnb;
    std::vector<Neighb> wknb;
    int iii1, iii2, iiat, lfnr, i_at, j_at;
    double dx, dy, dz, dd, ai, aisq, ajsq, aj, as, a;
    double xi, yi, zi, xs=0., ys=0., zs=0.;
    double dotarea, area, vol=0.;
    //double_ptr_t xus, dots=nullptr, atom_area=nullptr;
    std::vector<double> dots{};
    double_ptr_t atom_area=nullptr;

    int    nxbox, nybox, nzbox, nxy, nxyz;
    double xmin, ymin, zmin, xmax, ymax, zmax, ra2max, d, *pco;

    distribution = unsp_type(densit);
    if (distribution != -last_unsp || LAST_CUBUS != 4 || (densit != last_densit && densit != last_n_dot)) {
        if (make_unsp(densit, (-distribution), &N_DOT, 4)) return 1;
    }
    auto xus = XPUNSP;

    dotarea = FOURPI/double(N_DOT);
    area = 0.;


    /* start with neighbour list */
    /* calculate neighbour list with the box algorithm */
    if (nat==0) {
        //WARNING("nsc_dclm: no surface atoms selected");
        //return 1;
        throw std::domain_error("nsc_dclm: no surface atoms selected");
    }
    if (mode & FLAG_VOLUME) vol=0.;
    if (mode & FLAG_DOTS) {
        maxdots = 3 * N_DOT * nat / 10;
        //dots = (double *) CALLOC(maxdots, sizeof(double));
        dots.resize(maxdots);
        lfnr=0;
    }
    if (mode & FLAG_ATOM_AREA) {
        //atom_area = (double *) CALLOC(nat, sizeof(double));
        atom_area = new double[nat];
    }

    /* dimensions of atomic set, cell edge is 2*ra_max */
    xmin = co[0]; 
    xmax = xmin; xs=xmin;
    ymin = co[1]; 
    ymax = ymin; ys=ymin;
    zmin = co[2]; 
    zmax = zmin; zs=zmin;
    ra2max = radius[0];

    for (iat=1; iat<nat; iat++) {
        pco = co+3*iat;
        xmin = std::min(xmin, *pco);     
        xmax = std::max(xmax, *pco);
        ymin = std::min(ymin, *(pco+1)); 
        ymax = std::max(ymax, *(pco+1));
        zmin = std::min(zmin, *(pco+2)); 
        zmax = std::max(zmax, *(pco+2));
        xs= xs+ *pco; 
        ys = ys+ *(pco+1); 
        zs= zs+ *(pco+2);
        ra2max = std::max(ra2max, radius[iat]);
    }
    xs = xs/ (double) nat;
    ys = ys/ (double) nat;
    zs = zs/ (double) nat;
    ra2max = 2.*ra2max;

    d = xmax-xmin; nxbox = (int) std::max(ceil(d/ra2max), 1.);
    d = (((double)nxbox)*ra2max-d)/2.;
    xmin = xmin-d; xmax = xmax+d;
    d = ymax-ymin; nybox = (int) std::max(ceil(d/ra2max), 1.);
    d = (((double)nybox)*ra2max-d)/2.;
    ymin = ymin-d; ymax = ymax+d;
    d = zmax-zmin; nzbox = (int) std::max(ceil(d/ra2max), 1.);
    d = (((double)nzbox)*ra2max-d)/2.;
    zmin = zmin-d; zmax = zmax+d;
    nxy = nxbox*nybox;
    nxyz = nxy*nzbox;

    /* box number of atoms */
    //wkatm = (int *) CALLOC(3*nat, sizeof(int));
    std::vector<int> wkatm(3*nat,0);
    auto wkat1 = wkatm.begin() +nat;
    //wkdot = (int *) CALLOC(n_dot+nxyz+1, sizeof(int));
    std::vector<int> wkdot(N_DOT + nxyz + 1,0);
    auto wkbox = wkdot.begin() + N_DOT;

    for (iat=0; iat<nat; iat++) {
        pco = co+3*iat;
        i = (int) std::max(floor((  *pco  -xmin)/ra2max), 0.0);
        i = std::min(i,nxbox-1);
        j = (int) std::max(floor((*(pco+1)-ymin)/ra2max), 0.0);
        j = std::min(j,nybox-1);
        l = (int) std::max(floor((*(pco+2)-zmin)/ra2max), 0.0);
        l = std::min(l,nzbox-1);
        i = i+j*nxbox+l*nxy;
        wkat1[iat] = i;
        wkbox[i]++;
    }

    /* sorting of atoms in accordance with box numbers */
    j = wkbox[0];
    for (i=1; i<nxyz; i++) j = std::max(wkbox[i], j);
    for (i=1; i<=nxyz; i++) wkbox[i] += wkbox[i-1];
    //maxnei = std::min(nat, 27*j);
    maxnei = std::min(nat, 27*j);
    //wknb = (Neighb *) CALLOC(maxnei, sizeof(Neighb));
    wknb.resize(maxnei);
    for (iat=0; iat<nat; iat++) {
         wkatm[--wkbox[wkat1[iat]]] = iat;
    }


    /* calculate surface for all atoms, step cube-wise */
    for (iz=0; iz<nzbox; iz++) {
        iii = iz*nxy;
        izs = std::max(iz-1,0); ize = std::min(iz+2, nzbox);
        for (iy=0; iy<nybox; iy++) {
            ii = iy*nxbox+iii;
            iys = std::max(iy-1,0); 
            iye = std::min(iy+2, nybox);
            for (ix=0; ix<nxbox; ix++) {
                i = ii+ix;
                iii1=wkbox[i]; 
                iii2=wkbox[i+1];
                if (iii1 >= iii2) continue;
                
                ixs = std::max(ix-1,0); 
                ixe = std::min(ix+2, nxbox);

                iiat = 0;
                /* make intermediate atom list */
                for (jz=izs; jz<ize; jz++) {
                    jjj = jz*nxy;
                    for (jy=iys; jy<iye; jy++) {
                        jj = jy*nxbox+jjj;
                        for (jx=ixs; jx<ixe; jx++) {
                            j = jj+jx;
                            for (jat=wkbox[j]; jat<wkbox[j+1]; jat++) {
                                wkat1[iiat] = wkatm[jat]; iiat++;
                            }     /* end of cycle "jat" */
                        }       /* end of cycle "jx" */
                    }       /* end of cycle "jy" */
                }       /* end of cycle "jz" */
                for (iat=iii1; iat<iii2; iat++) {
                    i_at = wkatm[iat];
                    ai = radius[i_at]; 
                    aisq = ai*ai;
                    pco = co+3*i_at;
                    xi = *pco;
                    yi = *(pco+1);
                    zi = *(pco+2);
                    for (i=0; i < N_DOT; i++)
                        wkdot[i]=0;

                    auto ctnb = wknb.begin(); 
                    nnei = 0;
                    for (j=0; j<iiat; j++) {
                        j_at = *(wkat1+j);
                        if (j_at == i_at) continue;

                        aj = radius[j_at]; 
                        ajsq = aj*aj;
                        pco = co+3*j_at;
                        dx = *pco-xi; 
                        dy = *(pco+1)-yi; 
                        dz = *(pco+2)-zi;
                        dd = dx*dx+dy*dy+dz*dz;

                        as = ai+aj; if (dd > as*as) continue;
        
                        nnei++;
                        ctnb->x = dx; 
                        ctnb->y = dy; 
                        ctnb->z = dz;
                        ctnb->dot = (dd+aisq-ajsq)/(2.*ai); /* reference dot product */
                        ctnb++;
                    }

                    /* check points on accessibility */
                    if (nnei) {
                        last = 0;
                        i_ac = 0;
                        for (l=0; l < N_DOT; l++) {
          if (xus[3*l]*(wknb.begin()+last)->x+
              xus[1+3*l]*(wknb.begin()+last)->y+
              xus[2+3*l]*(wknb.begin()+last)->z <= (wknb.begin()+last)->dot) {
            for (j=0; j<nnei; j++) {
              if (xus[3*l]*(wknb.begin()+j)->x+xus[1+3*l]*(wknb.begin()+j)->y+
                  xus[2+3*l]*(wknb.begin()+j)->z > (wknb.begin()+j)->dot) {
                last = j; break;
                }
              }
            if (j >= nnei) { i_ac++; wkdot[l] = 1; }
            }     /* end of cycle j */
          }       /* end of cycle l */
        } else {
        i_ac  = N_DOT;
        for (l=0; l < N_DOT; l++)
            wkdot[l] = 1;
      }

      a = aisq*dotarea* (double) i_ac;
      area = area + a;
      if (mode & FLAG_ATOM_AREA) {
        atom_area[i_at] = a;
        }
      if (mode & FLAG_DOTS) {
        for (l=0; l < N_DOT; l++) {
          if (wkdot[l]) {
            lfnr++;
            if (maxdots <= 3*lfnr+1) {
              maxdots = maxdots + N_DOT * 3;
              dots.resize(maxdots);
              }
            dots[3*lfnr-3] = ai*xus[3*l]+xi;
            dots[3*lfnr-2] = ai*xus[1+3*l]+yi;
            dots[3*lfnr-1] = ai*xus[2+3*l]+zi;
            }
          }
        }
      if (mode & FLAG_VOLUME) {
        dx=0.; dy=0.; dz=0.;
        for (l=0; l < N_DOT; l++) {
          if (wkdot[l]) {
            dx=dx+xus[3*l];
            dy=dy+xus[1+3*l];
            dz=dz+xus[2+3*l];
            }
          }
        vol = vol+aisq*(dx*(xi-xs)+dy*(yi-ys)+dz*(zi-zs)+ai* (double) i_ac);
        }

      }         /* end of cycle "iat" */
    }           /* end of cycle "ix" */
    }           /* end of cycle "iy" */
    }           /* end of cycle "iz" */


  if (mode & FLAG_VOLUME) {
    vol = vol*FOURPI/(3.* (double) N_DOT);
    *value_of_vol = vol;
    }
  if (mode & FLAG_DOTS) {
    *nu_dots = lfnr;
    *lidots = &dots[0];
    }
  if (mode & FLAG_ATOM_AREA) {
    *at_area = atom_area;
    }
  *value_of_area = area;

  return 0;
  }














