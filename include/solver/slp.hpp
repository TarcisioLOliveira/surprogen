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

#ifndef SOLVER_SLP_HPP
#define SOLVER_SLP_HPP

#include "main.hpp"
#include <cmath>
#include <cstddef>

class SLP{
    public:
    /**
     * Update `x` based on derivatives and constraints.
     */
    void update(Vec& x, const Vec& dfdx, const Vec& g, const Vec& dgdx, const double xmin, const double xmax) const;

    private:
    const double EPS = 1e-7;
    /**
     * Dot product between two vectors.
     */
    inline Vec dot(const Vec& df, const Vec& dg, const size_t N,  const size_t M) const{
        Vec c(M, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < M; ++j){
                c[j] += df[i]*dg[i*M + j];
            }
        }
        return c;
    }
    /**
     * Multiplies a matrix and a vector.
     */
    inline Vec mat_dot_vec(const Vec& m, const Vec& v) const{
        const size_t l = m.size()/v.size();
        const size_t n = v.size();
        Vec r(l,0);
        for(size_t i = 0; i < l; ++i){
            for(size_t j = 0; j < n; ++j){
                r[i] += m[i*n + j]*v[j];
            }
        }
        return r;
    }
    /**
     * Scales a vector or matrix.
     */
    inline Vec scale(const double s, const Vec& v) const{
        auto sv = v;
        for(size_t i = 0; i < v.size(); ++i){
            sv[i] *= s;
        }
        return sv;
    }
    /**
     * Normalizes a vector.
     */
    inline Vec normalized(const Vec& v, const size_t N, const size_t M) const{
        Vec norm(M, 0);
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < M; ++j){
                norm[j] += v[i*M + j]*v[i*M + j];
            }
        }
        for(auto& n : norm){
            n = std::sqrt(n);
        }
        auto svn = v;
        for(size_t i = 0; i < N; ++i){
            for(size_t j = 0; j < M; ++j){
                svn[i*M + j] /= norm[j];
            }
        }
        return svn;
    }
};

#endif
