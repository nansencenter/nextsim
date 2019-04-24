#include "ensemble.hpp"

void ensemble::synopticPerturbation(int rdm)
{
    // Call fortran library for synoptic perturbation
    // // A file named ranfld.dat will be written

    p_pseudo2D_fld_sub();

}

void ensemble::addPerturbation(int rdm)
{
    float ranfld[rdm][8];

    // Call fortran library for synoptic perturbation
    // A file called ranfld.dat will be written

    M_ranfile = { "../IO/synforc.00", "../IO/synforc.01" };
    // Find and read the ranfld.dat into struct synoptic
    ifstream franfld;
    franfld.open(M_ranfile[0]);
    int row = 0;
    while(!franfld.eof()){
       std::string str;
       std::getline(franfld, str);
       std::stringstream ss(str);
       int col = 0;
       while(ss >> ranfld[row][col]) col++;

       row++;
    }
    franfld.close();
    synoptic.uwind.data.reserve(rdm);
    synoptic.vwind.data.reserve(rdm);
    synoptic.slp.data.reserve(rdm);
    synoptic.t2air.data.reserve(rdm);
    synoptic.precip.data.reserve(rdm);
    synoptic.relhum.data.reserve(rdm);
    for(int i=0; i<rdm; i++) {
        synoptic.uwind.data.push_back (ranfld[i][synoptic.uwind.id] );
        synoptic.vwind.data.push_back (ranfld[i][synoptic.vwind.id] );
        synoptic.t2air.data.push_back (ranfld[i][synoptic.t2air.id] );
        synoptic.slp.data.push_back   (ranfld[i][synoptic.slp.id]   );
        synoptic.precip.data.push_back(ranfld[i][synoptic.precip.id]);
        synoptic.relhum.data.push_back(ranfld[i][synoptic.relhum.id]);
    }
};


void ensemble::addPerturbation(std::vector<double>& loaded_uwind, std::vector<double>& loaded_vwind, int rdm, int ranid)
{
    float ranfld[rdm][8];

    // Find and read the ranfld.dat into struct synoptic
    M_ranfile = { "synforc.00", "synforc.01" };
    ifstream franfld;
    franfld.open(M_ranfile[ranid]);

    int row = 0;
    while(!franfld.eof()){
        std::string str;
        std::getline(franfld, str);
        std::stringstream ss(str);
        int col = 0;
        while(ss >> ranfld[row][col]) col++;
        row++;
    }
    franfld.close();

    synoptic.uwind.data.reserve (rdm);
    synoptic.vwind.data.reserve (rdm);
    synoptic.slp.data.reserve   (rdm);
    synoptic.t2air.data.reserve (rdm);
    synoptic.precip.data.reserve(rdm);
    synoptic.relhum.data.reserve(rdm);

    for(int i = 0; i < rdm; i++) {
        synoptic.uwind.data.push_back (ranfld[i][synoptic.uwind.id] );
        synoptic.vwind.data.push_back (ranfld[i][synoptic.vwind.id] );
        synoptic.t2air.data.push_back (ranfld[i][synoptic.t2air.id] );
        synoptic.slp.data.push_back   (ranfld[i][synoptic.slp.id]   );
        synoptic.precip.data.push_back(ranfld[i][synoptic.precip.id]);
        synoptic.relhum.data.push_back(ranfld[i][synoptic.relhum.id]);
    }

    for(int i = 0; i < rdm; i++) {
        loaded_uwind[i] += synoptic.uwind.data[i];
        loaded_vwind[i] += synoptic.vwind.data[i];
    }

};


void ensemble::computeMinMax(const std::vector<double> &ivector, const char* iname){
    double max = *max_element(ivector.begin(), ivector.end());
    double min = *min_element(ivector.begin(), ivector.end());
};

void ensemble::computeVecMean(const std::vector<double> &ivector, const char* iname){
    float average = std::accumulate(ivector.begin(), ivector.end(), 0.0)/ivector.size();
};


void ensemble::getpath(std::string iopath){
     M_ranpath = iopath;
};
