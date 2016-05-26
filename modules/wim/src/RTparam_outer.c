/* Standard files; */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* My header files; */
#include <RTparam_fast.h>
#include <RTparam_hardcoded.h>
#include <RTparam_outer.h>

#define PI M_PI
#define EPS 1.0e-12
#define MAXIT 100


int RTparam_outer(double outputs[],double h,double om,double visc_rp,double guess,double params[])
{

   double alp_nd,h_nd,zeta_nd;
   double Hw_nd,varpi,L;
   double damping_nd,visc_rp_nd;
   double H_nd = 4.0;/* infinite depth used in scattering calculation; */
   /* */
   double ki,kw,avc,BG1,BG2;
   double E,g,rhow,rhoi,nu,rho;
   double D;

   E     = params[0];
   g     = params[1];
   rhow  = params[2];
   rhoi  = params[3];
   nu    = params[4];
   rho   = rhoi/rhow;
   /* */
   D        = E*(h*h*h)/12/(1-nu*nu);
   L        = exp( 0.2*log(D/rhow/om/om) );
   alp_nd   = om*om/g*L;
   h_nd     = h/L;
   zeta_nd  = rho*h_nd;

   /* printf("E=          %f\n",h); */

   /* get wavenumber for ice; */
   varpi = 1/alp_nd-zeta_nd;
   /* [ki,BG2,avc]   = gen_root_ice(varpi,H_nd,guess*L); */
   gen_root_ice(&ki,&BG2,&avc,varpi,H_nd,guess*L);
   /* *kice = ki/L; */
   outputs[1] = ki/L;

   /* get wavenumber for water; */
   varpi = 1/alp_nd;
   Hw_nd = H_nd+zeta_nd;
   /* [kw,BG1] = gen_root_wtr(varpi,Hw_nd,alp_nd); */
   gen_root_wtr(&kw,&BG1,varpi,Hw_nd,alp_nd);
   /* *kwtr  = kw/L; */
   outputs[2] = kw/L;

   /* get intrinsic admittance; */
   /* |R|^2+int_adm*|T|^2=1 */
   /* printf("\noutput: BG1   = %0.9e\n", BG1); */
   /* printf("\noutput: BG2   = %0.9e\n", BG2); */
   /* int_adm  = BG1/BG2; */
   outputs[3] = BG1/BG2;

   /* get viscous attenuation; */
   visc_rp_nd  = visc_rp/rhow/om/L;
   damping_nd  = avc*visc_rp_nd;
   /* *damping    = damping_nd/L; */
   outputs[0] = damping_nd/L;

   /* printf("\L   = %0.9e\n", L); */
   /* printf("\noutput: BG1   = %0.9e\n", outputs[0]); */

   /* printf("outputs[0]  %d\n",outputs[0]); */

   /* get interpolated attenuation coefficient; */
   /* RTparam_fast(ac,modT,argR,argT, */
   /*              alp_nd,h_nd,*int_adm); */

   RTparam_fast(&outputs[4],&outputs[5],&outputs[6],&outputs[7],
                alp_nd,h_nd,outputs[3]);



   /* printf("inp[0]          %d\n",params[0]); */
   /* printf("inp[1]          %d\n",params[1]); */
   /* printf("inp[2]          %d\n",params[2]); */
   /* printf("inp[3]          %d\n",params[3]); */
   /* printf("inp[4]          %d\n",params[4]); */



   /* printf("out[0]          %d\n",outputs[0]); */
   /* printf("out[1]          %d\n",outputs[1]); */
   /* printf("out[2]          %d\n",outputs[2]); */
   /* printf("out[3]          %d\n",outputs[3]); */
   /* printf("out[4]          %d\n",outputs[4]); */
   /* printf("out[5]          %d\n",outputs[5]); */
   /* printf("out[6]          %d\n",outputs[6]); */
   /* printf("out[7]          %d\n",outputs[7]); */

}

/*******************************************************************/
/* %%function [ki,BG2,avc]=gen_root_ice(del,H,guess) */
/* finds the root of the ice dispersion relation nearest to 'guess'. */
/* Newton-Rhapson method solver */
/* of f  = Lam*k*sinh(k*H)-cosh(k*H) = 0, Lam=k^4+del */
int gen_root_ice(double *ki2, double *BG2,double *avc,
                 double del, double H,double guess) {

   double fac,k0,dk,ki;
   double res,Lam,Lampr,denom;

   fac   = 1.0;
   k0    = guess;

   /* Call dispersion relation function; */
   /* dk    = NR_corr_term(k0,del,H,fac); */
   NR_corr_term(&dk,&Lam,&Lampr,k0,del,H,fac);
   ki = k0-dk;

   while(fabs(dk) > EPS) {
     k0 = ki;

     /* Call dispersion relation function; */
     /* dk = NR_corr_term(k0,del,H,fac); */
     NR_corr_term(&dk,&Lam,&Lampr,k0,del,H,fac);
     ki = k0-dk;
   }

   /* Call dispersion relation function; */
   /* [dk,Lam,Lampr] = NR_corr_term(ki,del,H,fac); */
   NR_corr_term(&dk,&Lam,&Lampr,ki,del,H,fac);
   /* */
   denom = H*(Lam*Lam*ki*ki-1)+Lampr;
   res   = -ki/denom;

   /* Outputs; */
   *ki2  = ki;
   *BG2  = Lam*Lam*res;
   *avc  = ki/denom; /* -1i*dk/d(visc_rp) - non-dimensional */
}

/*******************************************************************/
/* %%function [kw,BG1]=gen_root_wtr(del,H,guess) */
/* finds the root of the water dispersion relation nearest to 'guess'. */
/* Newton-Rhapson method solver */
/* of f  = del*k*sinh(k*H)-cosh(k*H) = 0 */
int gen_root_wtr(double *kw2, double *BG1,
                 double del, double H,double guess) {

   double fac,k0,dk,kw;
   double res,Lam,Lampr,denom;

   fac   = 0.0;
   k0    = guess;

   /* Call dispersion relation function; */
   /* dk    = NR_corr_term(k0,del,H,fac); */
   NR_corr_term(&dk,&Lam,&Lampr,k0,del,H,fac);
   kw = k0-dk;

   while(fabs(dk) > EPS) {
      k0       = kw;

      /* Call dispersion relation function; */
      /* dk       = NR_corr_term(k0,del,H,fac); */
      NR_corr_term(&dk,&Lam,&Lampr,k0,del,H,fac);
      kw       = k0-dk;
   }

   /* Call dispersion relation function; */
   /* [dk,Lam,Lampr] = NR_corr_term(kw,del,H,fac); */
   NR_corr_term(&dk,&Lam,&Lampr,kw,del,H,fac);
   denom          = H*(Lam*Lam*kw*kw-1)+Lampr;
   res            = -kw/denom;

   /* printf("(denom,res,Bg2)=(%.9e,%.9e,%.9e)",res,denom,Lam*Lam*res); */

   /* Outputs; */
   *kw2  = kw;
   *BG1  = Lam*Lam*res;
}

/* function [dk,Lam,Lampr] = NR_corr_term(k,del,H,fac) */
/* %% dk=f/f_k, where f has the same zeros as of the dispersion function, */
/* %% is the correction term in the Newton-Rhapson method for finding zeros in f. */
int NR_corr_term(double *dk,double *Lam2,double *Lampr2,
                 double k, double del, double H,double fac) {

   double f,df;
   double Lam,Lampr,x,k4;

   k4    = k*k*k*k;
   Lam   = fac*k4+del;
   Lampr = 5*fac*k4+del;
   x     = 7.5;

   if(fabs(k*H)<=x) {
      f  = Lam*k*sinh(k*H)-cosh(k*H);
      df = Lam*(k*H)*cosh(k*H)+(Lampr-H)*sinh(k*H);
   }
   else {
      f  = Lam*k*tanh(k*H)-1;
      df = Lam*k*H+(Lampr-H)*tanh(k*H);
   }

   /* Outputs; */
   *dk      = f/df;
   *Lam2    = Lam;
   *Lampr2  = Lampr;
}
