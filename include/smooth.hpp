/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of texgen.
 *
 *   texgen is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   texgen is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with texgen.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef SMOOTH_HPP
#define SMOOTH_HPP

#include <cmath>

namespace smooth{

inline double abs(const double x, const double EPS){
    return std::sqrt(x*x + EPS);
}

inline double abs_deriv(const double x, const double EPS){
    return x/std::sqrt(x*x + EPS);
}

}

#endif
