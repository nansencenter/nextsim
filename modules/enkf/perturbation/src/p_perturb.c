/* main cpp to call perturbation on forcing and ic fields
 * this should be placed in nextsim for online perturbation
 */
#include "ensemble.hpp"

using namespace std;

int main(int argc, char** argv)
{

    /*
     * Here we pass the dimensions of the forcing fields to generate_ensemble
     */
    ensemble perturbation;
    int a = 129600;
    int b = 8;
    int ranstep = 0;

    perturbation.synopticPerturbation(a);
    perturbation.addPerturbation(a);

    cout<< "Size of uwind : \t" << perturbation.synoptic.uwind.data.size()  << '\n';
    cout<< "Size of vwind : \t" << perturbation.synoptic.vwind.data.size()  << '\n';
    cout<< "Size of t2air : \t" << perturbation.synoptic.t2air.data.size()  << '\n';
    cout<< "Size of slp   : \t" << perturbation.synoptic.slp.data.size()    << '\n';
    cout<< "Size of precip: \t" << perturbation.synoptic.precip.data.size() << '\n';
    cout<< "Size of relhum: \t" << perturbation.synoptic.relhum.data.size() << '\n';

    perturbation.computeMinMax( perturbation.synoptic.uwind.data, perturbation.synoptic.uwind.name );
    perturbation.computeVecMean(  perturbation.synoptic.uwind.data, perturbation.synoptic.uwind.name );
    perturbation.computeMinMax( perturbation.synoptic.vwind.data, perturbation.synoptic.vwind.name );
    perturbation.computeVecMean(  perturbation.synoptic.vwind.data, perturbation.synoptic.vwind.name );
    perturbation.computeMinMax( perturbation.synoptic.t2air.data, perturbation.synoptic.t2air.name );
    perturbation.computeVecMean(  perturbation.synoptic.t2air.data, perturbation.synoptic.t2air.name );
    perturbation.computeMinMax( perturbation.synoptic.slp.data  , perturbation.synoptic.slp.name   );
    perturbation.computeVecMean(  perturbation.synoptic.slp.data  , perturbation.synoptic.slp.name   );

    return 0;
}
