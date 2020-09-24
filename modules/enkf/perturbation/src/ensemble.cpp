#include "ensemble.hpp"

double *ensemble::synopticPerturbation(int const& ydim, int const& xdim, std::vector<std::vector<double> > &synforc,std::vector<std::vector<double> > &randfld, int const& perturbation_count) 
{
    double *synforc_p = (double *)malloc(synforc.size()*synforc[0].size()*sizeof(double));
    double *randfld_p = (double *)malloc(randfld.size()*randfld[0].size()*sizeof(double));
    //
    int Ncol = xdim*ydim;    // = synforc[0].size() 
    for(int row = 0; row < synforc.size(); row++){
        for(int col = 0; col < Ncol; col++) {
            synforc_p[row*Ncol+col] = synforc[row][col];
        }
    }
    for(int row = 0; row < randfld.size(); row++){
        for(int col = 0; col < Ncol; col++) {
            randfld_p[row*Ncol+col] = randfld[row][col];
        }
    } 
    //
    p_pseudo2D_fld_sub(&xdim, &ydim, &synforc_p[0], &randfld_p[0],&perturbation_count);
    //
    
    for(int col = 0; col < Ncol; col++) {
        for(int row = 0; row < synforc.size(); row++){
            synforc[row][col] = synforc_p[row*Ncol+col];
        }
        //std::cout<< col<< ",  "<<synforc[0][col]<<", "<<synforc[1][col]<<"\n";
    }           
    for(int row = 0; row < randfld.size(); row++){
        for(int col = 0; col < Ncol; col++) {
            randfld[row][col] = randfld_p[row*Ncol+col];
        }
    } 
    return(synforc_p);   
}    

void ensemble::synopticPerturbation(int const& ydim, int const& xdim, std::vector<std::vector<double> > &synforc,std::vector<std::vector<double> > &randfld, int const& perturbation_count) 
{
    std::cout<< "t1\n";
    int rows,cols, id;
    // rows = synforc.size();
    // cols = synforc[0].size();
//     std::cout<< "r,c:"<<rows<<", "<<cols<<"t1\n";
    // Call fortran library for synoptic perturbation
    // // A file named ranfld.dat will be written 
 
    // for(int i = 0; i < xdim; i++) {
    //     for(int j = 0; j < ydim; j++) {
    //         id = i*ydim + j;
    //         synforc[0][id] = i+1;
    //         synforc[1][id] = j+1;
    //     }    
    // }
    //std::cout<< "t2\n";
    //return;
    
    // pass 2D vector from c++ to fortran. Because fortran only recogizes array allocated contiguous in memory, while multidimentional array in c++ is not always contiguous automatically. Thus, a contiguous array is defined here to transfer data.

    // define two temporally continuous 1d array for transfering data with fortran routine and mpi broadcast       
    rows = synforc.size();
    cols = synforc[0].size();
    double *data1 = (double *)malloc(rows*cols*sizeof(double));
    double **synforc_p = (double **)malloc(rows*sizeof(double*));
    for(int i = 0; i < rows; i++) 
        synforc_p[i] = &(data1[cols*i]);
    //
    rows = randfld.size();
    cols = randfld[0].size();
    double *data2 = (double *)malloc(rows*cols*sizeof(double));
    double **randfld_p = (double **)malloc(rows*sizeof(double*));
    for(int i = 0; i < rows; i++) 
        randfld_p[i] = &(data2[cols*i]);

    //
    for(int i = 0; i < synforc[0].size(); i++) {
        for(int k = 0; k < synforc.size(); k++) {
            synforc_p[k][i] = synforc[k][i];
        }
    } 
    for(int i = 0; i < randfld[0].size(); i++) {
        for(int k = 0; k < randfld.size(); k++) {
            randfld_p[k][id] = randfld[k][id];
        }    
    }    
    p_pseudo2D_fld_sub(&xdim, &ydim, &synforc_p[0][0], &randfld_p[0][0],&perturbation_count);
    //
   // std::cout<<ydim<<", "<<xdim<<","<<perturbation_count<<"\n";
 //   return;
    //std::cout<< "t5\n";
    
    for(int i = 0; i < synforc[0].size(); i++) {
        for(int k = 0; k < synforc.size(); k++) {
            synforc[k][i] = synforc_p[k][i];
        }
        std::cout<<j+1<<", "<<i+1<< ",  "<< id << ",  "<<synforc[0][id]<<", "<<synforc[1][id]<<"\n";
    } 
    for(int i = 0; i < randfld[0].size(); i++) {
        for(int k = 0; k < randfld.size(); k++) {
            randfld[k][id] = randfld_p[k][id];
        }    
    }    
    std::cout<< "t6\n";
 //   return; */  
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





// // define two temporally continuous 1d array for transfering data with fortran routine and mpi broadcast       
//                 int rows = M_dataset->synforc.size();
//                 int cols = MN_full;
//                 double *temp1 = (double *)malloc(rows*cols*sizeof(double));
//                 double **synforc_p = (double **)malloc(rows*sizeof(double*));
//                 for(int i = 0; i < rows; i++) 
//                     synforc_p[i] = &(temp1[cols*i]);
                
//                 rows = M_dataset->randfld.size();
//                 cols = MN_full;
//                 double *temp2 = (double *)malloc(rows*cols*sizeof(double));
//                 double **randfld_p = (double **)malloc(rows*sizeof(double*));
//                 for(int i = 0; i < rows; i++) 
//                     randfld_p[i] = &(temp2[cols*i]);
//                 if (M_comm.rank() == 0) {                    
//                     LOG(DEBUG) << "### Generate perturbations based on the loaded wind inputs\n";        
//                     for(int k = 0; k < M_dataset->synforc.size(); k++)
//                          M_dataset->synforc[k].resize(MN_full);
//                          for(int i = 0; i < MN_full; i++)
//                             synforc_p[k][i] = synforc[k][i];
        
//                     for(int k = 0; k < M_dataset->randfld.size(); k++)
//                         M_dataset->randfld[k].resize(MN_full);
//                         for(int i = 0; i < MN_full; i++)
//                             randfld_p[k][i] = randfld[k][i];
//                     }
//                     // 
//                     perturbation.synopticPerturbation(M_full,N_full,&synforc_p[0][0],&randfld_p[0][0],M_dataset->perturbation_count); 
//                     M_dataset->perturbation_count++;
//                     //
//                     for(int k = 0; k < M_dataset->synforc.size(); k++)
//                          M_dataset->synforc[k].resize(MN_full);
//                          for(int i = 0; i < MN_full; i++)
//                             synforc[k][i] = synforc_p[k][i];
                    
//                     for(int k = 0; k < M_dataset->randfld.size(); k++)
//                         M_dataset->randfld[k].resize(MN_full);
//                         for(int i = 0; i < MN_full; i++)
//                             randfld[k][i] = randfld_p[k][i];
//                     //perturbation.synopticPerturbation(M_full,N_full,synforc00,synforc01,randfld01);  // may not need  synforc00
                    
//                     // LOG(DEBUG) << "### Load perturbations\n";
//                     // perturbation.loadPerturbation(synforc00,MN_full,0); 
//                    // perturbation.loadPerturbation(synforc01,MN_full,1); //todo: be replaced by returning variables from synopticPerturbation, but need to consider at the inital condition.                                                
//                 }