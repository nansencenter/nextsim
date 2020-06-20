/* int RTparam_outer(double *,double *,double *,double *, */
/*                  double *,double *,double *,double *, */
/*                  double,double,double,double,double *); */


int RTparam_outer(double outputs[],double h,double om,double visc_rp,double guess,double params[]);

int gen_root_ice(double *, double *,double *,
                 double, double,double);

int gen_root_wtr(double *, double *,
                 double, double,double);

int NR_corr_term(double *,double *,double *,
                 double, double, double,double);
