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

#include <random>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "matplotlibcpp.h"
#include "engine2D/random_discrete.hpp"
#include "roughness/discrete2D.hpp"
#include "solver/slp.hpp"
#include "solver/MMASolver.hpp"

namespace engine2D{

RandomDiscrete::RandomDiscrete(const double Ra_max, const double Rsk_max, const double Rku_max, const double z_min, const double z_max, const size_t N, const double tol):
    Ra_max(Ra_max), Rsk_max(Rsk_max), Rku_max(Rku_max), z_min(z_min), z_max(z_max), z_mid(0.5*(z_min+z_max)), tol(tol), N(N){}

void RandomDiscrete::optimize(){

    Vec z(N);
    Vec x(N);
    std::iota(x.begin(), x.end(), 0);

    const Vec zmin(N, z_min);
    const Vec zmax(N, z_max);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(z_min, z_max);
    for(auto& zi:z){
        zi = distribution(generator);
    }

    const size_t M = 6;

    MMASolver mma(N, M, 0, 1e11, 1);
    mma.SetAsymptotes(1, 0.9, 5);
    // SLP solver;

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

    Vec g{
        Ra_d,
        Rsk_d,
        Rku_d,
        -Ra_d,
        -Rsk_d,
        -Rku_d
    };

    Vec df(N,0);
    Vec dg(M*N, 0);
    for(size_t i = 0; i < N; ++i){
        df[i]       = 2*z_avg_d*dz_avg[i];
        dg[i*M + 0] = dRa[i];
        dg[i*M + 1] = dRsk[i];
        dg[i*M + 2] = dRku[i];
        dg[i*M + 3] = -dRa[i];
        dg[i*M + 4] = -dRsk[i];
        dg[i*M + 5] = -dRku[i];
    }

    std::cout << std::setprecision(15);

    std::cout << "z_avg: " << z_avg.get() << std::endl;
    std::cout << "Ra:    " << Ra.get() << std::endl;
    std::cout << "Rsk:   " << Rsk.get() << std::endl;
    std::cout << "Rku:   " << Rku.get() << std::endl;
    std::cout << std::endl;

    size_t iter = 0;
    do {
        // solver.update(z, df, g, dg, z_min, z_max);
        mma.Update(z.data(), df.data(), g.data(), dg.data(), zmin.data(), zmax.data());

        z_avg.update(z);
        Ra.update(&z_avg, z);
        Rq.update(&z_avg, z);
        Rsk.update(&z_avg, &Rq, z);
        Rku.update(&z_avg, &Rq, z);

        z_avg_d = z_avg.get() - z_mid;
        Ra_d = Ra.get() - Ra_max;
        Rsk_d = Rsk.get() - Rsk_max;
        Rku_d = Rku.get() - Rku_max;

        g[0] = Ra_d ;
        g[1] = Rsk_d;
        g[2] = Rku_d;
        g[3] = -Ra_d ;
        g[4] = -Rsk_d;
        g[5] = -Rku_d;
        
        for(size_t i = 0; i < N; ++i){
            df[i]       = 2*z_avg_d*dz_avg[i];
            dg[i*M + 0] = dRa[i];
            dg[i*M + 1] = dRsk[i];
            dg[i*M + 2] = dRku[i];
            dg[i*M + 3] = -dRa[i];
            dg[i*M + 4] = -dRsk[i];
            dg[i*M + 5] = -dRku[i];
        }

        std::cout << "z_avg: " << z_avg.get() << std::endl;
        std::cout << "Ra:    " << Ra.get() << std::endl;
        std::cout << "Rsk:   " << Rsk.get() << std::endl;
        std::cout << "Rku:   " << Rku.get() << std::endl;
        std::cout << std::endl;

        ++iter;
    } while(std::abs(Ra_d) > tol || std::abs(Rsk_d) > tol || std::abs(Rku_d) > tol);
    
    namespace plt = matplotlibcpp;

    std::ofstream file;

    file.open("z.txt");
    for(auto& zi:z){
        file << zi << std::endl;
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

    plt::figure();

    plt::ion();

    plt::plot(x, z);
    plt::xlabel("x [μm]");
    plt::ylabel("z [μm]");
    plt::tight_layout();

    plt::save("plot.pdf");

    plt::show(true);

    plt::close();
}

}
