//
// Created by Einar Ólason on 15/09/2022.
//

#ifndef NEXTSIM_INITIALISER_HPP
#define NEXTSIM_INITIALISER_HPP

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>

#include "timer.hpp"
#include "scattergather.hpp"

namespace Nextsim {

    class Initialiser : public ScatterGather
    {
    public:

        Initialiser(Timer::timer* timer)
                : vm(Environment::vm()),
                  ScatterGather(timer) {
            this->initOptAndParams();
        }

    protected:
        setup::IceType M_ice_type;

        void initModelState();

    private:
        void initOptAndParams();

        void initSlabOcean();
        void initIce();

        //ice-init functions
        void constantIce();
        void topazIce();
        void topazIceOsisafIcesat();
        void piomasIce();
        void nemoIce();
        void ciceIce();
        void topazForecastIce();
        void topazForecastAmsr2Ice();
        void topazForecastAmsr2OsisafIce();
        void topazForecastAmsr2OsisafNicIce(bool use_weekly_nic);
        void concBinsNic(double &young_conc_obs_min,double &young_conc_obs_max,double ci,bool use_weekly_nic);
        void cs2SmosIce();
        void cs2SmosAmsr2Ice();
        void smosIce();
        void glorys12Ice();

        //no ice-type option to activate these
        void topazAmsreIce();
        void topazAmsr2Ice();
        void amsr2ConstThickIce();

        void warrenClimatology();

        po::variables_map vm;
    };

}

#endif //NEXTSIM_INITIALISER_HPP
