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

#ifndef ENGINE2D_STEPS_HPP
#define ENGINE2D_STEPS_HPP

#include "engine2D.hpp"

namespace engine2D{

class Steps : public Engine2D{
    public:
    Steps(const double Ra_max, const double Rsk_max, const double Rku_max, const double z_min, const double z_max, const double w_max, const double h_max, const double alpha, const size_t Nz, const size_t Nparams, const double tol);
    virtual ~Steps() = default;

    virtual void optimize() override;
    private:
    const double Ra_max, Rsk_max, Rku_max, z_min, z_max, z_mid, w_max, h_max, alpha, tol;
    const size_t Nz, Nparams;

    void get_z(const Vec& x, Vec& z, const Vec& params) const;
    void get_dz(const Vec& x, const Vec& dRdz, Vec& out_dz, const Vec& params) const;
};

}


#endif
