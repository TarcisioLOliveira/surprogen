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

#ifndef TEXTURE2D_HPP
#define TEXTURE2D_HPP

#include "main.hpp"
#include "roughness/discrete2D.hpp"

class Engine2D{
    public:
    virtual ~Engine2D() = default;
    
    virtual void optimize() = 0;
};

#endif
