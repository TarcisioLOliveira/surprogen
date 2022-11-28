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

#include "solver/MMASolver.hpp"
#include "solver/slp.hpp"
#include <iostream>

void SLP::update(Vec& x, const Vec& dfdx, const Vec& g, const Vec& dgdx, const double xmin, const double xmax) const{
    // x_{i+1} = x_i + l_i*dfdx_i
    // dgdx_i*(x_{i+1}-x_i) <= g_i
    //
    // dgdx_i*(l_i*dfdx_i) <= g_i
    // l_i <= g_i/(dgdx_i*dfdx_i)
    const size_t N = x.size();
    const size_t M = g.size();
    
    // Normalize vectors (currently assuming single constraint)
    const auto norm_dfdx = this->normalized(dfdx, N, 1);
    const auto norm_dgdx = this->normalized(dgdx, N, M);

    MMASolver mma(1, 0, 0, 1e3, 1);
    mma.SetAsymptotes(0.5, 0.7, 1.2);

    // Calculate l_i
    const Vec dg_df = this->dot(norm_dfdx, norm_dgdx, N, M);
    double l_i = 0;//-g[0]/dg_df;
    for(size_t i = 0; i < M; ++i){
        l_i += -g[i]/dg_df[i];
    }
    l_i /= M;

    // Make sure l_i does not make x too low or too high.
    for(size_t i = 0; i < N; ++i){
        if(x[i] + l_i*norm_dfdx[i] < xmin){
            l_i = 0.5*(xmin-x[i])/norm_dfdx[i];
        } else if(x[i] + l_i*norm_dfdx[i] > xmax){
            l_i = 0.5*(xmax-x[i])/norm_dfdx[i];
        }
    }

    // Update x
    for(size_t i = 0; i < N; ++i){
        x[i] += l_i*norm_dfdx[i];
    }

}
