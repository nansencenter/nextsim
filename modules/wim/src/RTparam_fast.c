/* Standard files; */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* My header files; */
#include <RTparam_fast.h>
#include <RTparam_hardcoded.h>

#define PI M_PI
#define EPS 1.0e-8
#define MAXIT 100

/*********************************************************************/
/* Main routine here; */
int RTparam_fast(double *ac,double *modT,double *argR,double *argT,
                 double alp_nd,double hnd,double int_adm) {

  double t_a,t_h;            /* inputs to RTparam_inner; */
  double alp_lin3,alp_lin4;
  int    OPT,LOW;
  /* */
  double alp_nd_lims[6] = {1.0e-6,0.005,0.3,1.5,2.5,3.5};
  double mc_alplin[2]   = {-3.323529252398524,3.119943407349375};
  double y0_ll          = 40.0;
  double dy_ll          = 120.0;
  int    n_ll           = 3;
  double h1_ll          = .4;

  double hnd_lims[3]    = {1.0e-2,0.2,0.4};
  int    LOG_A_vec[5]   = {1,1,1,0,1};
  int    LOG_A;
  double a0,a1,h0,h1,l0,l1,dtmp;

  /* start of approximately linear regime; */
  alp_lin3       = mc_alplin[1]+mc_alplin[0]*log(hnd);
  alp_nd_lims[4] = alp_lin3;

  /* end of calculatable results; */
  dtmp      = cos(hnd/h1_ll*PI/2.0);
  alp_lin4  = y0_ll+dy_ll*exp(n_ll*log(dtmp));

  alp_nd_lims[5] = alp_lin4;

  /************************************************/
  /* determine if low/high thickness regime; */
  if(hnd<=hnd_lims[0]) {
     hnd    = hnd_lims[0];
     LOW    = 1;
  }
  else if(hnd>=hnd_lims[2]) {
     hnd    = hnd_lims[2];
     LOW    = 0;
  }

  /* now between hnd_lims[0] and hnd_lims[1]; */
  else if(hnd<hnd_lims[1]) {
     LOW = 1;
  }
  else {
     LOW = 0;
  }

  /* Interpolation value (-1<=t<=1) to use for hnd; */
  h0  = hnd_lims[1-LOW];
  h1  = hnd_lims[2-LOW];

  if(LOW==0) {
	  /* linear interpolation wrt hnd; */
	  t_h   = -1.0+2*(hnd-h0)/(h1-h0);
  }
  else {
	  /* interpolation wrt log(hnd); */
	  /* printf("\n(h0,h1)=(%f,%f)\n",h0,h1); */
	  l0  = log(h0);
	  l1  = log(h1);
	  t_h = -1.0+2*(log(hnd)-l0)/(l1-l0);
  }
  /* printf("\n(t_a,t_h)=(%f,%f)\n",t_a,t_h); */
  /************************************************/

  /************************************************/
  /* determine frequency regime; */
  if(alp_nd<=alp_nd_lims[0]) {
     alp_nd = alp_nd_lims[0];
     OPT    = 1;
  }
  else if(alp_nd>=alp_nd_lims[5]) {
     alp_nd = alp_nd_lims[5];
     OPT    = 5;
  }
  /* */
  else if(alp_nd<alp_nd_lims[1]) {
     OPT = 1;
  }
  else if(alp_nd<alp_nd_lims[2]) {
     OPT = 2;
  }
  else if(alp_nd<alp_nd_lims[3]) {
     OPT = 3;
  }
  else if(alp_nd<alp_nd_lims[4]) {
     OPT = 4;
  }
  else {
     OPT = 5;
  }
  /************************************************/

  a0     = alp_nd_lims[OPT-1];
  a1     = alp_nd_lims[OPT];
  LOG_A  = LOG_A_vec[OPT-1];
  if(LOG_A==1) {
	  /* interpolation wrt log(alp_nd); */
	  l0  = log(a0);
	  l1  = log(a1);
	  t_a = -1.0+2*(log(alp_nd)-l0)/(l1-l0);
  }
  else {
	  /* interpolation wrt alp_nd; */
	  t_a = -1.0+2*(alp_nd-a0)/(a1-a0);
  }

  /* call main interpolation routine; */
  /* (ac,modT,argR,argT) are pointers already; */
  RTparam_inner(ac,modT,argR,argT,
                t_a,t_h,OPT,LOW,int_adm);
}

int RTparam_inner(double *ac2,double *modT2,double *argR2,double *argT2,
                  double tx,double ty,int OPT,int LOW,double int_adm) {

  /***********************************************************/
  /* declare "input/output variables"; */
  double z;
  double ac,modR,modT,argR,argT;
  double Rr,Ri,Tr,Ti;

  /* auxiliary variables; */
  double *chebys;/* use malloc/free later */
  double *Zthings;/* use malloc/free later */
  int    Ncx,Ncy,ncol;
  int    i,j,r;
  int    Nthings;
  int    INTERP_MODE[5] = {1,1,3,2,1};
  int    IM;

  /* Order of chebyshev expansion for lower thicknesses; */
  int Ncx_L[5] = {10,10,10,10,3}; /* for alp_nd; */
  int Ncy_L[5] = {10,10,10,10,10};/* for hnd; */

  /* Order of chebyshev expansion for higher thicknesses; */
  int Ncx_H[5] = {10,10,10,10,4};/* for alp_nd; */
  int Ncy_H[5] = {10,10,10,10,10};/* for hnd; */
  /***********************************************************/

  /* Mode of interpolation; */
  IM  =INTERP_MODE[OPT-1];
  if(IM==1) {
	  /* Interpolate log(ac), Arg[R,T]; */
	  Nthings   = 3;
  }
  else if(IM==2) {
	  /* Interpolate (ac), Arg[R,T]; */
	  Nthings   = 3;
  }
  else if(IM==3) {
	  /* Interpolate Re[R,T], Im[R,T]; */
	  Nthings   = 4;
  }

  /* printf("LOW=%d\n",LOW); */

  if (LOW==1) {
	 /* Lower non-dimensional thicknesses; */
	  Ncx = Ncx_L[OPT-1];
	  Ncy = Ncy_L[OPT-1];

	  /*******************************************************/
	  if(OPT==1) {
        /**************************************************/
        /* Actual interpolation; */

        /* get vector of chebyshev coefficients into; */
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        /* get interpolated quantities; */
        Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L1(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }
     /*******************************************************/
     else if(OPT==2) {
	     /**************************************************/
        /* Actual interpolation; */

        /* get vector of chebyshev coefficients into; */
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        /* get interpolated quantities; */
        Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L2(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }
     /*******************************************************/
     else if(OPT==3) {
        /**************************************************/
        /* Actual interpolation; */

        /* get vector of chebyshev coefficients into; */
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        /* get interpolated quantities; */
        Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L3(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }
     /*******************************************************/
     else if(OPT==4) {
        /**************************************************/
        /* Actual interpolation; */

        /* get vector of chebyshev coefficients into; */
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        /* get interpolated quantities; */
        Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L4(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }

     /*******************************************************/
     else if(OPT==5) {
        /**************************************************/
        /* Actual interpolation; */

        /* get vector of chebyshev coefficients into; */
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        /* get interpolated quantities; */
        Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L5(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }
  } else {
	 /* Higher non-dimensional thicknesses; */
	  Ncx = Ncx_H[OPT-1];
	  Ncy = Ncy_H[OPT-1];

     /*******************************************************/
     if(OPT==1) {
        /**************************************************/
        /* Actual interpolation; */

        /* get vector of chebyshev coefficients into; */
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        /* get interpolated quantities; */
        Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H1(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }
     /*******************************************************/
     else if(OPT==2) {
	     /**************************************************/
	     /* Actual interpolation; */

	     /* get vector of chebyshev coefficients into; */
	     chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

	     /* get interpolated quantities; */
	     Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H2(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }
     /*******************************************************/
     else if(OPT==3) {
        /**************************************************/
        /* Actual interpolation; */

        /* get vector of chebyshev coefficients into; */
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        /* get interpolated quantities; */
        Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H3(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }
     /*******************************************************/
     else if(OPT==4) {
	     /**************************************************/
	     /* Actual interpolation; */

	     /* get vector of chebyshev coefficients into; */
	     chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

	     /* get interpolated quantities; */
	     Zthings   = malloc( Nthings*sizeof(double) );

	     /* printf("OPT=%d\n",OPT); */
	     /* printf("IM=%d\n",IM); */
	     /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

	     for(ncol=0;ncol<Nthings;ncol++) {
		     Amn_fxn_H4(chebys,ncol); /* NB this is specific to (OPT,LOW) */
		     OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
		     Zthings[ncol]   = z;
	     }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }

     /*******************************************************/
     else if(OPT==5) {
        /**************************************************/
        /* Actual interpolation; */

        /* get vector of chebyshev coefficients into; */
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        /* get interpolated quantities; */
        Zthings   = malloc( Nthings*sizeof(double) );

        /* printf("OPT=%d\n",OPT); */
        /* printf("IM=%d\n",IM); */
        /* printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy); */

        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H5(chebys,ncol); /* NB this is specific to (OPT,LOW) */
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        /* free memory from chebys; */
        free(chebys);
        /**************************************************/
     }
  }

  /********************************************************/
  /* interpret interpolated quantities; */
  RTparam_get_ac(ac2,modT2,argR2,argT2,Zthings,IM,int_adm);
  /********************************************************/

  /* free memory from Zthings; */
  free(Zthings);
}

/*********************************************************
*********************************************************/
int OP_chebinterp2d(double *z,double tx,double ty,
                     double *chebys,int Ncx,int Ncy) {

  double Amn, am,z0;
  double Tm0_x, Tm1_x, Tm_x;
  double Tn0_y, Tn1_y, Tn_y;
  int    s,nx,ny;

  /* printf("\n(in OP_chebinterp2d)(tx,ty)=(%f,%f)\n",tx,ty); */
  z0     = 0.0;
  Tm0_x  = 1.0;
  Tm1_x  = tx;
  /* printf("\ntx,Tm1_x,Tm0_x=%f,%f,%f\n",tx,Tm1_x,Tm0_x); */

  /**************************************************
  * Need to calculate z=\sum_{m=0}^Ncx a_m*T_m(tx),
  *  where a_m=\sum_{n=0}^Ncy A_{mn}*T_n(ty);
  **************************************************/

  /* m=0 term; */
  nx    = 0;
  Tm_x  = Tm0_x;
  am    = 0.0;
  /* */
  ny    = 0;
  Tn0_y = 1.0;
  s     = nx + ny*(Ncx+1); /* NB s starts from 0; */
  Amn   = chebys[s];
  am    = am+Amn*Tn0_y;
  /* printf("\n{nx,ny,s}=%d,%d,%d\n", nx,ny,s); */
  /* */
  /* printf("\ns=%d, Amn=%f,Tn=%f,am=%f\n", s,Amn,Tn0_y,am); */
  /* */
  ny    = 1;
  Tn1_y = ty;
  s     = nx + ny*(Ncx+1); /* NB s starts from 0; */
  Amn   = chebys[s];
  am    = am+Amn*Tn1_y;
  /* printf("\n{nx,ny,s}=%d,%d,%d\n", nx,ny,s); */
  /* */
  /* printf("\ns=%d, Amn=%f,Tn=%f,am=%f\n", s,Amn,Tn1_y,am); */
  /* */
  for(ny=2;ny<=Ncy;ny++) {
     Tn_y  = 2*ty*Tn1_y-Tn0_y;
     /* printf("\nny,ty,Tn_y,Tn1_y,Tn0_y=%d,%f,%f,%f,%f\n",ny,ty,Tn_y,Tn1_y,Tn0_y); */
     /* */
     Tn0_y = Tn1_y;
     Tn1_y = Tn_y;
     /* */
     s     = nx + ny*(Ncx+1); /* NB s starts from 0; */
     /* printf("\n{nx,ny,s}=%d,%d,%d\n", nx,ny,s); */
     /* */
     Amn   = chebys[s];
     am    = am+Amn*Tn_y;
     /* printf("\ns=%d, Amn=%f,Tn=%f,am=%f\n", s,Amn,Tn_y,am); */
  }

  /* update z; */
  z0 = z0+am*Tm_x;
  /* printf("z (m=0): %f\n",z0); */

  /*******************************************************
  * finished m=0 term;
  * now do m=1 term;
  *******************************************************/

  nx    = 1;
  Tm_x  = Tm1_x;
  am    = 0.0;
  /* */
  ny    = 0;
  Tn0_y = 1.0;
  s     = nx + ny*(Ncx+1);
  Amn   = chebys[s];
  am    = am+Amn*Tn0_y;
  /* */
  ny    = 1;
  Tn1_y = ty;
  s     = nx + ny*(Ncx+1);
  Amn   = chebys[s];
  am    = am+Amn*Tn1_y;

  for(ny=2;ny<=Ncy;ny++) {
     Tn_y  = 2*ty*Tn1_y-Tn0_y;
     Tn0_y = Tn1_y;
     Tn1_y = Tn_y;
     /* */
     s     = nx + ny*(Ncx+1); /* NB s starts from 0; */
     Amn   = chebys[s];
     am    = am+Amn*Tn1_y;
  }

  /* update z; */
  z0 = z0+am*Tm_x;
  /* printf("z (m=1): %f\n",z0); */

  /*******************************************************
  * Finished m=1 term;
  * Now sum over rest of terms (m=2 to m=Ncx);
  *******************************************************/
  for(nx=2;nx<=Ncx;nx++) {
     Tm_x  = 2*tx*Tm1_x-Tm0_x;
     Tm0_x = Tm1_x;
     Tm1_x = Tm_x;
     am    = 0.0;
     /* */
     ny    = 0;
     Tn0_y = 1.0;
     s     = nx + ny*(Ncx+1);
     Amn   = chebys[s];
     am    = am+Amn*Tn0_y;
     /* */
     ny    = 1;
     Tn1_y = ty;
     s     = nx + ny*(Ncx+1);
     Amn   = chebys[s];
     am    = am+Amn*Tn1_y;
     /* */
     for(ny=2;ny<=Ncy;ny++) {
        Tn_y  = 2*ty*Tn1_y-Tn0_y;
        Tn0_y = Tn1_y;
        Tn1_y = Tn_y;
        /* */
        s     = nx + ny*(Ncx+1); /* NB s starts from 0; */
        Amn   = chebys[s];
        am    = am+Amn*Tn1_y;
     }

     /* update z; */
     z0 = z0+am*Tm_x;
     /* printf("z (m=%d): %f\n",nx,z0); */
   }

   *z = z0;
   /* printf("\n(in OP_chebinterp2d) output (z) = %f\n",*z); */

   return 0;
}

/***********************************************/
int RTparam_get_ac(double *ac2,double *modT2,double *argR2,double *argT2,
                   double *Zthings,int IM,double int_adm) {

  double ac,modR,modT,argR,argT;
  double Rr,Ri,Tr,Ti;

  if(IM==1) {
     ac     = exp(Zthings[0]);
     argR   = Zthings[1];
     argT   = Zthings[2];
     modT   = sqrt(exp(-ac/2.0)/int_adm);
  }
  else if(IM==2) {
     ac     = Zthings[0];
     argR   = Zthings[1];
     argT   = Zthings[2];
     modT   = sqrt(exp(-ac/2.0)/int_adm);
  }
  else if(IM==3) {
     Rr     = Zthings[0];
     Ri     = Zthings[1];
     Tr     = Zthings[2];
     Ti     = Zthings[3];
     /* */
     argR   = atan2(Ri,Rr);
     argT   = atan2(Ti,Tr);
     modR   = sqrt(Rr*Rr+Ri*Ri);
     modT   = sqrt(Tr*Tr+Ti*Ti);
     ac     = -2*log(1-modR*modR);
  }

  *ac2   = ac;
  *modT2 = modT;
  *argR2 = argR;
  *argT2 = argT;

   return 0;
}
/***********************************************/
