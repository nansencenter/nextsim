#include <oasis_cpp_interface.h>

extern "C"
{
  int oasis3_init_comp(int *, const char *, int *);
  int oasis3_terminate();

  int oasis3_def_partition(int *, int *, int *);
  int oasis3_def_var(int *, const char *, int *, int *, int *, int *, int *, int *);
  int oasis3_enddef();

  int oasis3_put_1d(int *, int *, double *, int *);
  int oasis3_put_2d(int *, int *, double *, int *, int *);
  int oasis3_get_1d(int *, int *, double *, int *);
  int oasis3_get_2d(int *, int *, double *, int *, int *);
  // tm 19/10/2015
  int oasis3_abort(int *, const char *, const char *, int *, int *);
  int oasis3_get_localcomm(int *);
  int oasis3_create_couplcomm(int *, int *, int *);
  int oasis3_start_grids_writing(int *);
  int oasis3_write_grid(const char *, int *, int *, double *, double *, int*);
  int oasis3_write_corner(const char *, int *, int *, int *, double *, double *, int*);
  int oasis3_write_area(const char *, int *, int *, double *, int*);
  int oasis3_write_mask(const char *, int *, int *, int *, int*);
  int oasis3_terminate_grids_writing();
};


int OASIS3::init_comp(int * compid, std::string model_name)
{
    int lname = model_name.length();
    return oasis3_init_comp(compid, model_name.c_str(), &lname);
};

int OASIS3::terminate()
{
    return oasis3_terminate();
};

int OASIS3::def_partition(int * il_part_id, int * ig_paral, int ig_paral_len)
{
    return oasis3_def_partition(il_part_id,ig_paral,&ig_paral_len);
};

int OASIS3::def_var(int * var_id,std::string name, int il_part_id, int * var_nodims, int kinout, int * var_actual_shape, int var_type)
{
    int lname = name.length();
    return oasis3_def_var(var_id,name.c_str(),&il_part_id,var_nodims,&kinout,var_actual_shape,&var_type, &lname);
};

int OASIS3::enddef()
{
    return oasis3_enddef();
};

int OASIS3::put_1d(int var_id, int date, double * field_array, int field_array_len)
{
    return oasis3_put_1d(&var_id,&date,field_array,&field_array_len);
};

int OASIS3::put_2d(int var_id, int date, double * field_array, int field_array_len_x, int field_array_len_y)
{
    return oasis3_put_2d(&var_id,&date,field_array,&field_array_len_x,&field_array_len_y);
};

int OASIS3::get_1d(int var_id, int date, double * field_array, int field_array_len)
{
    return oasis3_get_1d(&var_id,&date,field_array,&field_array_len);
};

int OASIS3::get_2d(int var_id, int date, double * field_array, int field_array_len_x, int field_array_len_y)
{
    return oasis3_get_2d(&var_id,&date,field_array,&field_array_len_x,&field_array_len_y);
};

// tm 19/10/2015
int OASIS3::abort(int comp_id, std::string routine_name, std::string abort_message)
{
  int lname, lmess;
  lname = routine_name.length();
  lmess = abort_message.length();
  return oasis3_abort(&comp_id, routine_name.c_str(), abort_message.c_str(), &lname, &lmess);

}

int OASIS3::get_localcomm(MPI_Comm *localcomm)
{
  int fortran_localcomm;
  int ierror;
  // communicateur fortran donne par OASIS
  ierror = oasis3_get_localcomm(&fortran_localcomm);

  // conversion du communicateur MPI fortran (int) en un communicateur C
  *localcomm = MPI_Comm_f2c(fortran_localcomm);
  return ierror;
}

int OASIS3::create_couplcomm(bool active, MPI_Comm *localcomm, MPI_Comm *cplcomm)
{
  int fortran_cplcomm;
  MPI_Fint fortran_localcomm;
  int ierror;
  int icpl = int(active);

  // conversion du communicateur MPI fortran (int) en un communicateur C
  fortran_localcomm = MPI_Comm_c2f(*localcomm);

  // communicateur fortran donne par OASIS
  ierror = oasis3_create_couplcomm(&icpl, &fortran_localcomm, &fortran_cplcomm);

  // conversion du communicateur MPI fortran (int) en un communicateur C
  *cplcomm   = MPI_Comm_f2c(fortran_cplcomm);
  return ierror;
}

int OASIS3::write_grid(std::string gridname, int nx, int ny, double *lon, double *lat)
{
  int lgridname = gridname.length();
  return oasis3_write_grid(gridname.c_str(), &nx, &ny, lon, lat, &lgridname);
}

int OASIS3::write_corner(std::string gridname, int nx, int ny, int nc, double *clon, double *clat)
{
  int lgridname = gridname.length();
  return oasis3_write_corner(gridname.c_str(), &nx, &ny,&nc, clon, clat, &lgridname);
}

int OASIS3::write_area(std::string gridname, int nx, int ny, double *area)
{
  int lgridname = gridname.length();
  return oasis3_write_area(gridname.c_str(), &nx, &ny, area,  &lgridname);
}

int OASIS3::write_mask(std::string gridname, int nx, int ny, int *mask)
{
  int lgridname = gridname.length();
  return oasis3_write_mask(gridname.c_str(), &nx, &ny, mask, &lgridname);
}

int OASIS3::start_grids_writing(int flag)
{
  return oasis3_start_grids_writing(&flag);
}


int OASIS3::terminate_grids_writing()
{
  return oasis3_terminate_grids_writing();
}
