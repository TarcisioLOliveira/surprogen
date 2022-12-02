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
#include "engine2D/sinusoids.hpp"

int main(int argc, char* argv[]){

    engine2D::Sinusoids r(1.5, -1, 4, -5, 5, -15.0*8/700, 15.0*8/700, 500, 700, 1e-3);

    r.optimize();

    return 0;
}
