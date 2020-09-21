#include "ensemble.hpp"

void ensemble::synopticPerturbation(int const& xdm, int const& ydm, std::vector<std::vector<double> > &synforc,std::vector<std::vector<double> > &randfld) 
{
    // Call fortran library for synoptic perturbation
    // // A file named ranfld.dat will be written   
    p_pseudo2D_fld_sub(&xdm, &ydm, &synforc[0][0], &randfld[0][0]);
}

// void ensemble::addPerturbation(int rdm)
// {
//     std::vector<std::vector<double> > ranfld(rdm, std::vector<double>(8));

//     // Call fortran library for synoptic perturbation
//     // A file called ranfld.dat will be written

//     M_ranfile = { "../IO/synforc.00", "../IO/synforc.01" };
//     // Find and read the ranfld.dat into struct synoptic
//     ifstream franfld;
//     franfld.open(M_ranfile[0]);
//     int row = 0;
//     while(!franfld.eof()){
//         std::string str;
//         std::getline(franfld, str);
//         std::stringstream ss(str);   
//         int col = 0;     
//         while(ss >> ranfld[row][col]) col++;

//        row++;
//     }
//     franfld.close();
//     synoptic.uwind.data.reserve(rdm);
//     synoptic.vwind.data.reserve(rdm);
//     synoptic.slp.data.reserve(rdm);
//     synoptic.t2air.data.reserve(rdm);
//     synoptic.precip.data.reserve(rdm);
//     synoptic.relhum.data.reserve(rdm);   
//     for(int i=0; i<rdm; i++) {
//         synoptic.uwind.data.push_back (ranfld[i][synoptic.uwind.id] );
//         synoptic.vwind.data.push_back (ranfld[i][synoptic.vwind.id] );
//         synoptic.t2air.data.push_back (ranfld[i][synoptic.t2air.id] );
//         synoptic.slp.data.push_back   (ranfld[i][synoptic.slp.id]   );
//         synoptic.precip.data.push_back(ranfld[i][synoptic.precip.id]);
//         synoptic.relhum.data.push_back(ranfld[i][synoptic.relhum.id]);
//     }
// };


// void ensemble::addPerturbation(std::vector<double>& velocity_u, std::vector<double>& velocity_v, int rdm, int ranid)
// {
//     std::vector<std::vector<double> > ranfld(rdm, std::vector<double>(8,0.0)); 

//     // Find and read the ranfld.dat into struct synoptic
//     M_ranfile = { "synforc.00", "synforc.01" };
//     ifstream franfld;
//     franfld.open(M_ranfile[ranid]);

//     int row = 0;    
//     while(!franfld.eof()){
//         std::string str;
//         std::getline(franfld, str);
//         std::stringstream ss(str);
//         int col = 0;
//         while(ss >> ranfld[row][col]) col++;
//         row++;
//     }
//     franfld.close();

//     synoptic.uwind.data.reserve (rdm);
//     synoptic.vwind.data.reserve (rdm);
//     synoptic.slp.data.reserve   (rdm);
//     synoptic.t2air.data.reserve (rdm);
//     synoptic.precip.data.reserve(rdm);
//     synoptic.relhum.data.reserve(rdm);

//     for(int i = 0; i < rdm; i++) {
//         synoptic.uwind.data.push_back (ranfld[i][synoptic.uwind.id] );
//         synoptic.vwind.data.push_back (ranfld[i][synoptic.vwind.id] );
//         synoptic.t2air.data.push_back (ranfld[i][synoptic.t2air.id] );
//         synoptic.slp.data.push_back   (ranfld[i][synoptic.slp.id]   );
//         synoptic.precip.data.push_back(ranfld[i][synoptic.precip.id]);
//         synoptic.relhum.data.push_back(ranfld[i][synoptic.relhum.id]);
//     }

//     for(int i = 0; i < rdm; i++) {
//         velocity_u[i] += synoptic.uwind.data[i];
//         velocity_v[i] += synoptic.vwind.data[i];
//     }

// };


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


//void ensemble::loadPerturbation(std::vector<double>& uwind, std::vector<double>& vwind, int rdm, int ranid)
// todo: remove this function
void ensemble::loadPerturbation(std::vector<std::vector<float> > &synforc,int rdm, int ranid)
{
    std::vector<std::vector<float> > ranfld(rdm, std::vector<float>(8,0.0)); 

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
    for(int i = 0; i < rdm; i++) {
        synforc[0][i] = ranfld[i][3];  //see index in synforc_wr()
        synforc[1][i] = ranfld[i][4];
    }
};

// todo:  ensure the start,count are correct, check with others.
void ensemble::addPerturbation(std::vector<double>& velocity_u, std::vector<double>& velocity_v, std::vector<double>& synforc_u, std::vector<double>& synforc_v, int x_start, int y_start, int x_count, int y_count)
{    
    //int count = x_count*y_count; =velocity_u.size() //
    int count = velocity_u.size();
    int start = x_count*y_start + x_start; // double check
    
    //std::cout<<"subdomain_size="<<synforc_v.size()<<",start="<<start<<", end="<<start+count -1<<", length="<<count<<", x_start="<<x_start<<", y_start="<<y_start<<",x_count="<<x_count<<", y_count="<<y_count<<"\n";
    for(int i = start; i <=count; i++) {
        velocity_u[i] += synforc_u[i];
        velocity_v[i] += synforc_v[i];                
    }
};


// //----------------------------------------------// how to declare type Dataset in this file?
// void ensemble::addPerturbation(Dataset *M_dataset, \   
//     std::vector<std::vector<double> >& ranfld00,\
//     std::vector<std::vector<double> >& ranfld01)
// {
//     int y_start = M_dataset->grid.dimension_y_start;
//     int x_start = M_dataset->grid.dimension_x_start;
//     int y_count = M_dataset->grid.dimension_y_count;
//     int x_count = M_dataset->grid.dimension_x_count;
    
//     int start = x_count*(y_start-1) + x_start; // double check
//     int n = 0;
//     for(int i = start; i < x_count*y_count; i++) {
//         // 
//         M_dataset->variables[0].loaded_data[0][n] += ranfld00[0][i];  // wind u
//         M_dataset->variables[1].loaded_data[0][n] += ranfld00[1][i];  // wind v 
//         //
//         M_dataset->variables[0].loaded_data[1][n] += ranfld01[0][i];  // wind u
//         M_dataset->variables[1].loaded_data[1][n] += ranfld01[1][i];  // wind v 
//         n++;
//     }
// };
