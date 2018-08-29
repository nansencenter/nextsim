#include "mpi.h"
#include <string>

// Class OASIS3 provides access to all OASIS3 procedures via the C++-Fortran interface
class OASIS3
{
 public : static const int CLIM_Apple          =    1;
 public : static const int CLIM_BadName        =   -2;
 public : static const int CLIM_BadPort        =   -3;
 public : static const int CLIM_BadTaskId      =  -15;
 public : static const int CLIM_BadType        =   -4;
 public : static const int CLIM_Box            =    2;
 public : static const int CLIM_Clength        =   32;
 public : static const int CLIM_ContPvm        =    0;
 public : static const int CLIM_Double         =    8;
 public : static const int CLIM_Down           =  -20;
 public : static const int CLIM_FastExit       =   -1;
 public : static const int CLIM_FirstCall      =  -12;
 public : static const int CLIM_Group          =  -14;
 public : static const int CLIM_In             =    1;
 public : static const int CLIM_IncSize        =   -8;
 public : static const int CLIM_IncStep        =   -7;
 public : static const int CLIM_InitBuff       =  -17;
 public : static const int CLIM_InOut          =    2;
 public : static const int CLIM_Integer        =    1;
 public : static const int CLIM_LdX            =    5;
 public : static const int CLIM_Length         =    3;
 public : static const int CLIM_MaxCodes       =  -22;
 public : static const int CLIM_MaxSegments    =   338;
 public : static const int CLIM_Mpi            =  -22;
 public : static const int CLIM_NoTask         =  -16;
 public : static const int CLIM_NotClim        =   -9;
 public : static const int CLIM_NotStep        =   -6;
 public : static const int CLIM_Offset         =    2;
 public : static const int CLIM_Ok             =    0;
 public : static const int CLIM_Orange         =    3;
 public : static const int CLIM_Out            =    0;
 public : static const int CLIM_Pack           =  -18;
 public : static const int CLIM_ParSize        =    2*CLIM_MaxSegments+2;
 public : static const int CLIM_PbRoute        =  -13;
 public : static const int CLIM_Pvm            =  -11;
 public : static const int CLIM_PvmExit        =  -21;
 public : static const int CLIM_Real           =    4;
 public : static const int CLIM_Segments       =    2;
 public : static const int CLIM_Serial         =    0;
 public : static const int CLIM_SizeX          =    3;
 public : static const int CLIM_SizeY          =    4;
 public : static const int CLIM_StopPvm        =    1;
 public : static const int CLIM_Strategy       =    1;
 public : static const int CLIM_TimeOut        =  -10;
 public : static const int CLIM_Unpack         =  -19;
 public : static const int CLIM_Void           =    0;
  
 public : static const int ip_accumul          =    3;
 public : static const int ip_auxilary         =    7;
 public : static const int ip_average          =    2;
 public : static const int ip_exported         =    1;
 public : static const int ip_expout           =    5;
 public : static const int ip_ignored          =    2;
 public : static const int ip_ignout           =    6;
 public : static const int ip_input            =    3;
 public : static const int ip_instant          =    1;
 public : static const int ip_max              =    5;
 public : static const int ip_min              =    4;
 public : static const int ip_output           =    4;
  
 public : static const int PRISM_Double        =    8;
 public : static const int PRISM_DoubleDef     =   -5;
 public : static const int PRISM_FromRest      =   10;
 public : static const int PRISM_FromRestOut   =   13;
 public : static const int PRISM_In            =   21;
 public : static const int PRISM_InOut         =    2;
 public : static const int PRISM_Input         =   11;
 public : static const int PRISM_LocTrans      =    5;
 public : static const int PRISM_NotDef        =   -2;
 public : static const int PRISM_NotFreq       =  -23;
 public : static const int PRISM_Ok            =    0;
 public : static const int PRISM_Out           =   20;
 public : static const int PRISM_Output        =    7;
 public : static const int PRISM_Real          =    4;
 public : static const int PRISM_Recvd         =    3;
 public : static const int PRISM_RecvOut       =   12;
 public : static const int PRISM_Sent          =    4;
 public : static const int PRISM_SentOut       =    8;
 public : static const int PRISM_ToRest        =    6;
 public : static const int PRISM_ToRestOut     =    9;
  
 public : static const int OASIS_Double        = PRISM_Double;
 public : static const int OASIS_DoubleDef     = PRISM_DoubleDef;
 public : static const int OASIS_FromRest      = PRISM_FromRest;
 public : static const int OASIS_FromRestOut   = PRISM_FromRestOut;
 public : static const int OASIS_In            = PRISM_In;
 public : static const int OASIS_InOut         = PRISM_InOut;
 public : static const int OASIS_Input         = PRISM_Input;
 public : static const int OASIS_LocTrans      = PRISM_LocTrans;
 public : static const int OASIS_NotDef        = PRISM_NotDef;
 public : static const int OASIS_NotFreq       = PRISM_NotFreq;
 public : static const int OASIS_Ok            = PRISM_Ok;
 public : static const int OASIS_Out           = PRISM_Out;
 public : static const int OASIS_Output        = PRISM_Output;
 public : static const int OASIS_Real          = PRISM_Real;
 public : static const int OASIS_Recvd         = PRISM_Recvd;
 public : static const int OASIS_RecvOut       = PRISM_RecvOut;
 public : static const int OASIS_Sent          = PRISM_Sent;
 public : static const int OASIS_SentOut       = PRISM_SentOut;
 public : static const int OASIS_ToRest        = PRISM_ToRest;
 public : static const int OASIS_ToRestOut     = PRISM_ToRestOut;
  
 public : static int init_comp(int *, std::string);
 public : static int terminate();
  
 public : static int def_partition(int *, int *, int);
 public : static int def_var(int *, std::string, int, int *, int, int *, int);
 public : static int enddef();
  
 public : static int put_1d(int, int, double *, int);
 public : static int put_2d(int, int, double *, int, int);
 public : static int get_1d(int, int, double *, int);
 public : static int get_2d(int, int, double *, int, int);
  
    // ajout  tm 19/10/2015
 public : static int abort(int ,  std::string, std::string);
 public : static int get_localcomm(MPI_Comm *);
 public : static int create_couplcomm(bool, MPI_Comm *);
 public : static int start_grids_writing(int *);
 public : static int start_grids_writing(int);
 public : static int terminate_grids_writing();
 public : static int write_grid(std::string, int, int, double *, double *);
 public : static int write_corner(std::string, int, int, int, double *, double *);
 public : static int write_area(std::string, int, int, double *);
 public : static int write_mask(std::string, int, int, int *);

};
