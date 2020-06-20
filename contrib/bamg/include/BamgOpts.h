/*!\file:  BamgOpts.h
 * \brief place holder for optimization function arguments
 */ 

#ifndef _BAMGOPTS_H_
#define _BAMGOPTS_H_

class BamgOpts{

	public:

		/*Parameters*/
		double  anisomax;
		double  cutoff;
		double  coeff;
		double  errg;
		double  gradation;
		int     Hessiantype;
		double  MaxCornerAngle;
		int     maxnbv;
		double  maxsubdiv;
		int     Metrictype;
		int     nbjacobi;
		int     nbsmooth;
		double  omega;
		double  power;
		bool    random;
		int     verbose;

		/*Flags*/
		int     Crack;
		int     geometricalmetric;
		int     KeepVertices;
		int     splitcorners;

		/*Metric related*/
		double  hmin;
		double  hmax;
		int     hminVerticesSize[2];
		double* hminVertices;
		int     hmaxVerticesSize[2];
		double* hmaxVertices;
		int     hVerticesSize[2];
		double* hVertices;
		int     metricSize[2];
		double* metric;
		int     fieldSize[2];
		double* field;
		int     errSize[2];
		double* err;

		BamgOpts();
		~BamgOpts();

		void Check(void);

};
#endif
