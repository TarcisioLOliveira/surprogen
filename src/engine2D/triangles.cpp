/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of surprogen.
 *
 *   surprogen is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   surprogen is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with surprogen.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <cmath>
#include <random>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "matplotlibcpp.h"
#include "engine2D/triangles.hpp"
#include "roughness/discrete2D.hpp"
#include "solver/MMASolver.hpp"
#include <fenv.h>

namespace engine2D{

inline const double FLOORD = 1e-15;
inline const double SP_MULT = 100;
inline const double MIN_MULT = 50;

inline double softplus(const double x){
    if(x < -1){
        return 0;
    } else if(x < 1){
        return std::log(1 + std::exp(SP_MULT*x))/SP_MULT;
    } else {
        return x;
    }
}

inline double softplus_deriv(const double x){
    if(x < -1){
        return 0;
    } else if(x < 1){
        return 1.0/(1.0 + std::exp(-SP_MULT*x));
    } else {
        return 1.0;
    }
}

inline double smin(const double x1, const double x2){
    const double ex1 = std::exp(-MIN_MULT*x1);
    const double ex2 = std::exp(-MIN_MULT*x2);

    return (x1*ex1 + x2*ex2)/(ex1 + ex2);
}

inline double smin_deriv(const double x1, const double x2, const double dx1, const double dx2){
    const double ex1 = std::exp(x1);
    const double ex2 = std::exp(x2);

    const double sm = (x1*ex1 + x2*ex2)/(ex1 + ex2);

    return (ex1*(1 - MIN_MULT*(x1 - sm))*dx1 + ex2*(1 - MIN_MULT*(x2 - sm))*dx2)/(ex1 + ex2);
}

inline double tri_base(const double x, const double y){
    return smin(softplus(x), softplus(y));
}

inline double tri_base_deriv(const double x, const double y, const double dx, const double dy){
    return smin_deriv(softplus(x), softplus(y), softplus_deriv(x)*dx, softplus_deriv(y)*dy);
}

inline double tri(const double x, const double w, const double h, const double l){
    if(std::abs(x - l) < 1.1*w/2){
      return (h*2/w)*tri_base(x - l + w/2, -x + l + w/2);
    } else {
        return 0;
    }
}

inline double dtri_h(const double x, const double w, const double h, const double l){
    if(std::abs(x - l) < 1.1*w/2){
      return (2/w)*tri_base(x - l + w/2, -x + l + w/2);
    } else {
        return 0;
    }
}
inline double dtri_w(const double x, const double w, const double h, const double l){
    if(std::abs(x - l) < 1.1*w/2){
      return (h*2/w)*tri_base_deriv(x - l + w/2, -x + l + w/2, 1.0/2, 1.0/2) - (h*2/(w*w))*tri_base(x - l + w/2, -x + l + w/2);
    } else {
        return 0;
    }
}
inline double dtri_l(const double x, const double w, const double h, const double l){
    if(std::abs(x - l) < 1.1*w/2){
      return (h*2/w)*tri_base_deriv(x - l + w/2, -x + l + w/2, -1, 1);
    } else {
        return 0;
    }
}

Triangles::Triangles(const double Ra_max, const double Rsk_max, const double Rku_max, const double z_min, const double z_max, const double w_max, const double h_max, const size_t Nz, const size_t Nparams, const double tol):
    Ra_max(Ra_max), Rsk_max(Rsk_max), Rku_max(Rku_max), z_min(z_min), z_max(z_max), z_mid(0.5*(z_min+z_max)), w_max(w_max), h_max(h_max), tol(tol), Nz(Nz), Nparams(Nparams){}

void Triangles::optimize(){

    namespace plt = matplotlibcpp;

    const size_t Mp = 3;

    Vec x(Nz);
    Vec z(Nz);
    std::iota(x.begin(), x.end(), 0);

    Vec dz_avgdp(Nparams*Mp);
    Vec dRadp(Nparams*Mp);
    Vec dRskdp(Nparams*Mp);
    Vec dRkudp(Nparams*Mp);

    Vec pmin(Nparams*Mp, 0);
    Vec pmax(Nparams*Mp, 0);
    Vec params(Nparams*Mp);

    std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> distr_h(-h_max, h_max);
    std::uniform_real_distribution<double> distr_w(2, w_max);
    std::uniform_real_distribution<double> distr_l(0, Nz);
    for(size_t i = 0; i < Nparams; ++i){
        pmin[i            ] = -h_max;
        pmin[i +   Nparams] = 2;
        pmin[i + 2*Nparams] = 0;

        pmax[i            ] = h_max;
        pmax[i +   Nparams] = w_max;
        pmax[i + 2*Nparams] = Nz;

        params[i            ] = distr_h(generator);
        params[i + Nparams  ] = distr_w(generator);
        params[i + 2*Nparams] = distr_l(generator);
    }

    const size_t M = 6;

    MMASolver mma(Nparams*Mp, M, 0, 1e5, 1);
    mma.SetAsymptotes(100, 0.7, 1.2);

    get_z(x, z, params);

    discrete2D::AverageHeight z_avg(z);
    discrete2D::RoughnessAverage Ra(&z_avg, z, 1e-10);
    discrete2D::RoughnessRMS Rq(&z_avg, z);
    discrete2D::RoughnessSkewness Rsk(&z_avg, &Rq, z);
    discrete2D::RoughnessKurtosis Rku(&z_avg, &Rq, z);

    double z_avg_d = z_avg.get() - z_mid;
    double Ra_d = Ra.get() - Ra_max;
    double Rsk_d = Rsk.get() - Rsk_max;
    double Rku_d = Rku.get() - Rku_max;

    const Vec& dz_avg = z_avg.get_deriv();
    const Vec& dRa = Ra.get_deriv();
    const Vec& dRsk = Rsk.get_deriv();
    const Vec& dRku = Rku.get_deriv();

    Vec g(2*M, 0);
    Vec df(Nparams*Mp,0);
    Vec dg(M*Nparams*Mp, 0);

    // Prevent matplotlib from throwing the window to top at every update
    plt::backend("TkAgg");

    plt::figure();

    plt::Plot plot("Profile", x, z);
    plt::xlabel("x [μm]");
    plt::ylabel("z [μm]");
    plt::tight_layout();
    plt::xlim(0.0, (double)Nz);
    plt::ylim(z_min, z_max);

    plt::show(false);
    plt::pause(0.0001);

    std::cout << std::setprecision(15);

    //feenableexcept(FE_ALL_EXCEPT);

    size_t iter = 0;
    do {
        get_dz(x, dz_avg, dz_avgdp, params);
        get_dz(x, dRa, dRadp, params);
        get_dz(x, dRsk, dRskdp, params);
        get_dz(x, dRku, dRkudp, params);

        g[0] = Ra_d ;
        g[1] = Rsk_d;
        g[2] = Rku_d;
        g[3] = -Ra_d ;
        g[4] = -Rsk_d;
        g[5] = -Rku_d;
        
        for(size_t i = 0; i < Nparams*Mp; ++i){
            df[i]       = 2*z_avg_d*dz_avgdp[i];
            dg[i*M + 0] = dRadp[i];
            dg[i*M + 1] = dRskdp[i];
            dg[i*M + 2] = dRkudp[i];
            dg[i*M + 3] = -dRadp[i];
            dg[i*M + 4] = -dRskdp[i];
            dg[i*M + 5] = -dRkudp[i];
        }

        plot.update(x, z);
        plt::draw();
        plt::pause(0.0001);

        std::cout << "z_avg: " << z_avg.get() << std::endl;
        std::cout << "Ra:    " << Ra.get() << std::endl;
        std::cout << "Rsk:   " << Rsk.get() << std::endl;
        std::cout << "Rku:   " << Rku.get() << std::endl;
        std::cout << std::endl;

        if(std::abs(Ra_d) < tol && std::abs(Rsk_d) < tol && std::abs(Rku_d) < tol){
            break;
        }

        // solver.update(z, df, g, dg, z_min, z_max);
        mma.Update(params.data(), df.data(), g.data(), dg.data(), pmin.data(), pmax.data());

        get_z(x, z, params);

        z_avg.update(z);
        Ra.update(&z_avg, z);
        Rq.update(&z_avg, z);
        Rsk.update(&z_avg, &Rq, z);
        Rku.update(&z_avg, &Rq, z);

        z_avg_d = z_avg.get() - z_mid;
        Ra_d = Ra.get() - Ra_max;
        Rsk_d = Rsk.get() - Rsk_max;
        Rku_d = Rku.get() - Rku_max;

        ++iter;
    } while(true);


    std::ofstream file;

    file.open("z.txt");
    for(auto& zi:z){
        file << zi << std::endl;
    }
    file.close();

    file.open("params.txt");
    const double* h = params.data();
    const double* w = h + Nparams;
    const double* l = w + Nparams;
    file << "h:" << std::endl;
    for(size_t i = 0; i < Nparams; ++i){
        file << h[i] << std::endl;
    }
    file << "w:" << std::endl;
    for(size_t i = 0; i < Nparams; ++i){
        file << w[i] << std::endl;
    }
    file << "l:" << std::endl;
    for(size_t i = 0; i < Nparams; ++i){
        file << l[i] << std::endl;
    }
    file.close();

    file.open("roughness.txt");
    file << std::setprecision(15);
    file << "z_avg: " << z_avg.get() << std::endl;
    file << "Ra:    " << Ra.get() << std::endl;
    file << "Rsk:   " << Rsk.get() << std::endl;
    file << "Rku:   " << Rku.get() << std::endl;
    file << std::endl;
    file.close();

    plt::save("plot.pdf");

    plt::close();
}

void Triangles::get_z(const Vec& x, Vec& z, const Vec& params) const{
    const double* h = params.data();
    const double* w = h + Nparams;
    const double* l = w + Nparams;
    #pragma omp parallel for
    for(size_t i = 0; i < Nz; ++i){
        z[i] = 0;
        for(size_t j = 0; j < Nparams; ++j){
            z[i] += tri(x[i], w[j], h[j], l[j]);
        }
    }
}

void Triangles::get_dz(const Vec& x, const Vec& dRdz, Vec& out_dz, const Vec& params) const{
    const double* h = params.data();
    const double* w = h + Nparams;
    const double* l = w + Nparams;
    double* dh = out_dz.data();
    double* dw = dh + Nparams;
    double* dl = dw + Nparams;
    #pragma omp parallel for
    for(size_t j = 0; j < Nparams; ++j){
        dh[j] = 0;
        dw[j] = 0;
        dl[j] = 0;
        for(size_t i = 0; i < Nz; ++i){
            dh[j] += dRdz[i]*dtri_h(x[i], w[j], h[j], l[j]);
            dw[j] += dRdz[i]*dtri_w(x[i], w[j], h[j], l[j]);
            dl[j] += dRdz[i]*dtri_l(x[i], w[j], h[j], l[j]);
        }
    }
}

}
