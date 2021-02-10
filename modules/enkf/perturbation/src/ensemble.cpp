#include "ensemble.hpp"

void ensemble::synopticPerturbation(std::vector<double> &synforc,std::vector<double> &randfld, int const& ydim, int const& xdim, int const& synforc_exist) 
{
    // M_full = ydim is size of y domain, N_full = xdim is size of x domain
    // Call fortran library for synoptic perturbation
    p_pseudo2D_fld_sub(&xdim, &ydim, &synforc[0], &randfld[0],&synforc_exist);    
}    


// void ensemble::addPerturbation(std::vector<double>& velocity_u, std::vector<double>& velocity_v, std::vector<double>& synforc, int M_full, int N_full, int x_start, int y_start, int x_count, int y_count)
// {   // submesh size is y_count*x_count; =velocity_u.size()     
//     // full mesh size M_full(y domain)*N_full(x domain)
//     // outer loop rows indicated by y direction. j is index in submesh
//     // interior loop read a row of data for submesh using starting from index x_start_tmp in the full mesh, i is index in the full mesh.
//     int start_tmp,n = 0;
//     for(int j = 0; j <y_count; j++) {  
//         start_tmp = x_start + (y_start + j)*N_full;
//         for(int i = start_tmp; i < start_tmp + x_count; i++ ) {    
//             velocity_u[n] += synforc[i];
//             velocity_v[n] += synforc[i+M_full*N_full];       
//             n++; 
//         }
//     }        
// };


void ensemble::addPerturbation(std::vector<double>& perturbed_field, std::vector<double>& synforc, int M_full, int N_full, int x_start, int y_start, int x_count, int y_count,int index)
{   // submesh size is y_count*x_count; =velocity_u.size()     
    // full mesh size M_full(y domain)*N_full(x domain)
    // outer loop rows indicated by y direction. j is index in submesh
    // interior loop read a row of data for submesh using starting from index x_start_tmp in the full mesh, i is index in the full mesh.
    
    // index of forcing field defined in save_randfld_synforc() of mod_random_forcing.f90, starting from 1
    int start_tmp,n = 0;
    for(int j = 0; j <y_count; j++) {  
        start_tmp = x_start + (y_start + j)*N_full;
        for(int i = start_tmp; i < start_tmp + x_count; i++ ) {    
            perturbed_field[n] += synforc[i+(index-1)*M_full*N_full];       
            n++; 
        }
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

// void ensemble::synopticPerturbation(int const& ydim, int const& xdim, std::vector<std::vector<double> > &synforc,std::vector<std::vector<double> > &randfld, int const& synforc_exist, double *synforc_p) 
// {
//     double *randfld_p = (double *)malloc(randfld.size()*randfld[0].size()*sizeof(double));
//     //
//     int Ncol = xdim*ydim;    // = synforc[0].size() 
//     for(int row = 0; row < synforc.size(); row++){
//         for(int col = 0; col < Ncol; col++) {
//             synforc_p[row*Ncol+col] = synforc[row][col];
//         }
//     }
//     for(int row = 0; row < randfld.size(); row++){
//         for(int col = 0; col < Ncol; col++) {
//             randfld_p[row*Ncol+col] = randfld[row][col];
//         }
//     } 
//     //
//     p_pseudo2D_fld_sub(&xdim, &ydim, &synforc_p[0], &randfld_p[0],&synforc_exist);
//     //
    
//     for(int col = 0; col < Ncol; col++) {
//         for(int row = 0; row < synforc.size(); row++){
//             synforc[row][col] = synforc_p[row*Ncol+col];
//         }
//         //std::cout<< col<< ",  "<<synforc[0][col]<<", "<<synforc[1][col]<<"\n";
//     }           
//     for(int row = 0; row < randfld.size(); row++){
//         for(int col = 0; col < Ncol; col++) {
//             randfld[row][col] = randfld_p[row*Ncol+col];
//         }
//     } 
// }    

// void ensemble::synopticPerturbation(int const& ydim, int const& xdim, std::vector<std::vector<double> > &synforc,std::vector<std::vector<double> > &randfld, int const& synforc_exist) 
// {
//     std::cout<< "t1\n";
//     int rows,cols, id;
//     // rows = synforc.size();
//     // cols = synforc[0].size();
// //     std::cout<< "r,c:"<<rows<<", "<<cols<<"t1\n";
//     // Call fortran library for synoptic perturbation
//     // // A file named ranfld.dat will be written 
 
//     // for(int i = 0; i < xdim; i++) {
//     //     for(int j = 0; j < ydim; j++) {
//     //         id = i*ydim + j;
//     //         synforc[0][id] = i+1;
//     //         synforc[1][id] = j+1;
//     //     }    
//     // }
//     //std::cout<< "t2\n";
//     //return;
    
//     // pass 2D vector from c++ to fortran. Because fortran only recogizes array allocated contiguous in memory, while multidimentional array in c++ is not always contiguous automatically. Thus, a contiguous array is defined here to transfer data.

//     // define two temporally continuous 1d array for transfering data with fortran routine and mpi broadcast       
//     rows = synforc.size();
//     cols = synforc[0].size();
//     double *data1 = (double *)malloc(rows*cols*sizeof(double));
//     double **synforc_p = (double **)malloc(rows*sizeof(double*));
//     for(int i = 0; i < rows; i++) 
//         synforc_p[i] = &(data1[cols*i]);
//     //
//     rows = randfld.size();
//     cols = randfld[0].size();
//     double *data2 = (double *)malloc(rows*cols*sizeof(double));
//     double **randfld_p = (double **)malloc(rows*sizeof(double*));
//     for(int i = 0; i < rows; i++) 
//         randfld_p[i] = &(data2[cols*i]);

//     //
//     for(int i = 0; i < synforc[0].size(); i++) {
//         for(int k = 0; k < synforc.size(); k++) {
//             synforc_p[k][i] = synforc[k][i];
//         }
//     } 
//     for(int i = 0; i < randfld[0].size(); i++) {
//         for(int k = 0; k < randfld.size(); k++) {
//             randfld_p[k][id] = randfld[k][id];
//         }    
//     }    
//     p_pseudo2D_fld_sub(&xdim, &ydim, &synforc_p[0][0], &randfld_p[0][0],&synforc_exist);
//     //
//    // std::cout<<ydim<<", "<<xdim<<","<<synforc_exist<<"\n";
//  //   return;
//     //std::cout<< "t5\n";
    
//     for(int i = 0; i < synforc[0].size(); i++) {
//         for(int k = 0; k < synforc.size(); k++) {
//             synforc[k][i] = synforc_p[k][i];
//         }
//         std::cout<<i+1<< ",  "<< id << ",  "<<synforc[0][id]<<", "<<synforc[1][id]<<"\n";
//     } 
//     for(int i = 0; i < randfld[0].size(); i++) {
//         for(int k = 0; k < randfld.size(); k++) {
//             randfld[k][id] = randfld_p[k][id];
//         }    
//     }    
//     std::cout<< "t6\n";
//  //   return; */  
// }    



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



// // todo: use this function for restart case
// void ensemble::loadPerturbation(double *synforc,int rdm, int ranid)
// {
//     // Find and read the ranfld.dat into struct synoptic
//     std::vector<std::vector<double> > ranfld(rdm, std::vector<double>(8,0.0)); 
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
//     for(int i = 0; i < rdm; i++) {
//         synforc[i]     = ranfld[i][3];  //see index in synforc_wr()
//         synforc[i+rdm] = ranfld[i][4];
//     }
// };
