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
	Nextsim::ensemble ens;
	int a = 129600;
	int b = 8;

	ens.generate_ensemble(a, b);

	cout<< "Size of uwind : \t" << ens.synoptic.uwind.data.size()  << '\n';
	cout<< "Size of vwind : \t" << ens.synoptic.vwind.data.size()  << '\n';
	cout<< "Size of t2air : \t" << ens.synoptic.t2air.data.size()  << '\n';
	cout<< "Size of slp   : \t" << ens.synoptic.slp.data.size()    << '\n';
	cout<< "Size of precip: \t" << ens.synoptic.precip.data.size() << '\n';
	cout<< "Size of relhum: \t" << ens.synoptic.relhum.data.size() << '\n';

	ens.compute_minmax( ens.synoptic.uwind.data, ens.synoptic.uwind.name );
	ens.compute_vmean(  ens.synoptic.uwind.data, ens.synoptic.uwind.name );
	ens.compute_minmax( ens.synoptic.vwind.data, ens.synoptic.vwind.name );
	ens.compute_vmean(  ens.synoptic.vwind.data, ens.synoptic.vwind.name );
	ens.compute_minmax( ens.synoptic.t2air.data, ens.synoptic.t2air.name );
	ens.compute_vmean(  ens.synoptic.t2air.data, ens.synoptic.t2air.name );
	ens.compute_minmax( ens.synoptic.slp.data  , ens.synoptic.slp.name   );
	ens.compute_vmean(  ens.synoptic.slp.data  , ens.synoptic.slp.name   );

	return 0;
}
