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

#ifndef ROUGHNESS_DISCRETE2D_HPP
#define ROUGHNESS_DISCRETE2D_HPP

#include <cstddef>
#include "main.hpp"
#include "smooth.hpp"

namespace discrete2D{

class AverageHeight{
    public:
    AverageHeight(const Vec& z):
        z_avg(update_z_avg(z)), dz_avg(init_dz_avg(z))
    {}

    inline void update(const Vec& z){
        z_avg = update_z_avg(z);
    }

    inline double get() const{
        return z_avg;
    }

    inline const Vec& get_deriv() const{
        return dz_avg;
    }

    private:
    double z_avg;
    const Vec dz_avg;

    inline double update_z_avg(const Vec& z){
        double A = 0;
        for(size_t i = 1; i < z.size(); ++i){
            A += (z[i] + z[i-1])/2;
        }
        return A/z.size();
    }
    inline Vec init_dz_avg(const Vec& z){
        Vec dz(z.size(),0);
        for(size_t i = 1; i < z.size(); ++i){
            dz[i-1] += 0.5;
            dz[i]   += 0.5;
        }
        return dz;
    }
};

class RoughnessAverage{
    public:
    RoughnessAverage(const AverageHeight* const z_avg, const Vec& z, const double ABS_EPS):
        ABS_EPS(ABS_EPS), Ra(0), dRa(z.size()){

        this->update(z_avg, z);
    }

    inline void update(const AverageHeight* const z_avg, const Vec& z){
        this->update_Ra(z_avg->get(), z);
        this->update_dRa(z_avg->get(), z_avg->get_deriv(), z);
    }

    inline double get() const{
        return this->Ra;
    }
    inline const Vec& get_deriv() const{
        return this->dRa;
    }

    private:
    const double ABS_EPS;
    double Ra;
    Vec dRa;

    inline void update_Ra(const double z_avg, const Vec& z){
        Ra = 0;
        for(size_t i = 0; i < z.size(); ++i){
            Ra += smooth::abs(z[i] - z_avg, ABS_EPS);
        }
        Ra /= z.size();
    }

    inline void update_dRa(const double z_avg, const Vec& dz_avg, const Vec& z){
        const size_t LEN = z.size();
        for(size_t i = 0; i < LEN; ++i){
            dRa[i] = smooth::abs_deriv(z[i] - z_avg, ABS_EPS)*(1.0 - dz_avg[i])/LEN;
        }
    }
};

class RoughnessRMS{
    public:
    RoughnessRMS(const AverageHeight* const z_avg, const Vec& z):
        Rq(0), dRq(z.size()){

        this->update(z_avg, z);
    }

    inline void update(const AverageHeight* const z_avg, const Vec& z){
        this->update_Rq(z_avg->get(), z);
        this->update_dRq(z_avg->get(), z_avg->get_deriv(), z);
    }

    inline double get() const{
        return this->Rq;
    }
    inline const Vec& get_deriv() const{
        return this->dRq;
    }

    private:
    double Rq;
    Vec dRq;

    inline void update_Rq(const double z_avg, const Vec& z){
        Rq = 0;
        for(size_t i = 0; i < z.size(); ++i){
            const double zz = z[i] - z_avg;
            Rq += zz*zz;
        }
        Rq = std::sqrt(Rq/z.size());
    }

    inline void update_dRq(const double z_avg, const Vec& dz_avg, const Vec& z){
        const size_t LEN = z.size();
        for(size_t i = 0; i < LEN; ++i){
            dRq[i] = (z[i] - z_avg)*(1.0 - dz_avg[i])/(Rq*LEN);
        }
    }
};

class RoughnessSkewness{
    public:
    RoughnessSkewness(const AverageHeight* const z_avg, const RoughnessRMS* Rq, const Vec& z):
        Rsk(0), dRsk(z.size()){

        this->update(z_avg, Rq, z);
    }

    inline void update(const AverageHeight* const z_avg, const RoughnessRMS* Rq, const Vec& z){
        this->update_Rsk(z_avg->get(), Rq->get(), z);
        this->update_dRsk(z_avg->get(), z_avg->get_deriv(), Rq->get(), Rq->get_deriv(), z);
    }

    inline double get() const{
        return this->Rsk;
    }
    inline const Vec& get_deriv() const{
        return this->dRsk;
    }

    private:
    double Rsk;
    Vec dRsk;

    inline void update_Rsk(const double z_avg, const double Rq, const Vec& z){
        Rsk = 0;
        for(size_t i = 0; i < z.size(); ++i){
            const double zz = z[i] - z_avg;
            Rsk += zz*zz*zz;
        }
        Rsk /= Rq*Rq*Rq*z.size();
    }

    inline void update_dRsk(const double z_avg, const Vec& dz_avg, const double Rq, const Vec& dRq, const Vec& z){
        const size_t LEN = z.size();
        for(size_t i = 0; i < LEN; ++i){
            const double zz = z[i] - z_avg;
            dRsk[i] = 3*(zz*zz*(1.0 - dz_avg[i])/(Rq*Rq*Rq*LEN) - Rsk*dRq[i]/Rq);
        }
    }
};

class RoughnessKurtosis{
    public:
    RoughnessKurtosis(const AverageHeight* const z_avg, const RoughnessRMS* Rq, const Vec& z):
        Rku(0), dRku(z.size()){

        this->update(z_avg, Rq, z);
    }

    inline void update(const AverageHeight* const z_avg, const RoughnessRMS* Rq, const Vec& z){
        this->update_Rku(z_avg->get(), Rq->get(), z);
        this->update_dRku(z_avg->get(), z_avg->get_deriv(), Rq->get(), Rq->get_deriv(), z);
    }

    inline double get() const{
        return this->Rku;
    }
    inline const Vec& get_deriv() const{
        return this->dRku;
    }

    private:
    double Rku;
    Vec dRku;

    inline void update_Rku(const double z_avg, const double Rq, const Vec& z){
        Rku = 0;
        for(size_t i = 0; i < z.size(); ++i){
            const double zz = z[i] - z_avg;
            Rku += zz*zz*zz*zz;
        }
        Rku /= Rq*Rq*Rq*Rq*z.size();
    }

    inline void update_dRku(const double z_avg, const Vec& dz_avg, const double Rq, const Vec& dRq, const Vec& z){
        const size_t LEN = z.size();
        for(size_t i = 0; i < LEN; ++i){
            const double zz = z[i] - z_avg;
            dRku[i] = 4*(zz*zz*zz*(1.0 - dz_avg[i])/(Rq*Rq*Rq*Rq*LEN) - Rku*dRq[i]/Rq);
        }
    }
};

}

#endif
