#include <cstdio>
#include <string.h>
#include <cmath>

#include "./shared.h"
//#include "../errorlib/include.h"
#include "./exceptions.h"

#include "GeomEdge.h"
#include "Geometry.h"

using namespace std;

namespace bamg {

	/*Constructor/Destructor*/

	/*Methods*/
	R2 GeomEdge::F(double theta) const{/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/F)*/
		// parametrization of the curve edge

	   R2 A=v[0]->r,B=v[1]->r;
		double ca,cb,cta,ctb;

		//Check that theta is in [0 1]
		_assert_(theta>-1e-12 && theta<1+1e-12);

		if (TgA()){
			if (TgB()){ //Hermite interpolation
				cb =  theta*theta*(3-2*theta);
				ca =  1-cb;
				cta = (1-theta)*(1-theta)*theta;
				ctb = (theta-1)*theta*theta ;
			}
			else {
				double t = theta;
				cb = t*t;
				ca = 1-cb;
				cta= t-cb;
				ctb=0;
			}
		}
		else{
			if (TgB()){
				double t = 1-theta;
				ca = t*t;
				cb = 1-ca;
				ctb= -t+ca;
				cta=0;
			}
			else { // lagrange P1
				ca =(1-theta);
				cb = theta;
				cta=ctb=0;
			}
		}
		return A*ca + B*cb + tg[0]*cta + tg[1]*ctb;
	  }
	/*}}}*/
	int    GeomEdge::Mark()    const  {/*{{{*/
		return type &16;
	}/*}}}*/
	double GeomEdge::R1tg(double theta,R2 & t) const{/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/R1tg)*/
		// 1/R of radius of cuvature

		R2 A=v[0]->r,B=v[1]->r;
		double dca,dcb,dcta,dctb;
		double ddca,ddcb,ddcta,ddctb;
		double tt = theta*theta;

		//check theta
		_assert_(theta>=0 && theta<=1);

		if (TgA()){
			if (TgB()){
				// Tangent A and B provided:
				// interpolation d'hermite
				dcb = 6*theta*(1-theta);
				ddcb = 6*(1-2*theta);
				dca = -dcb;
				ddca = -ddcb;
				dcta =  (3*theta - 4)*theta + 1;
				ddcta=6*theta-4;
				dctb = 3*tt - 2*theta;
				ddctb = 6*theta-2;
			}
			else {
				//Tangent A provided but tangent B not provided
				// 1-t*t, t-t*t, t*t
				double t = theta;
				dcb = 2*t;
				ddcb = 2;
				dca = -dcb;
				ddca = -2;
				dcta = 1-dcb;
				ddcta = -ddcb;
				dctb=0;
				ddctb=0;
			}
		}
		else{
			if (TgB()){
				//Tangent B provided but tangent A not provided
				double t = 1-theta;
				dca = -2*t;
				ddca = 2;
				dcb = -dca;
				ddcb = -2;
				dctb = 1+dca;
				ddctb= ddca;
				dcta =0;
				ddcta =0;
			}
			else {
				//Neither thangent A nor tangent B provided
				// lagrange P1
				t=B-A;
				return 0;
			}
		}
		R2 d  =  A*dca  + B*dcb  + tg[0]* dcta  + tg[1] * dctb;
		R2 dd =  A*ddca + B*ddcb + tg[0]* ddcta + tg[1] * ddctb;
		double d2=(d,d);
		double sd2 = sqrt(d2);
		t=d;
		if(d2>1.0e-20){
			t/=sd2;
			return Abs(Det(d,dd))/(d2*sd2);
		}
		else return 0;
	}
	/*}}}*/
	int    GeomEdge::Required()       {/*{{{*/
		return type &64;
	}/*}}}*/
	void GeomEdge::Set(const GeomEdge & rec,const Geometry & Gh ,Geometry & GhNew){ /*{{{*/
		*this = rec;
		v[0] = GhNew.vertices + Gh.GetId(v[0]);
		v[1] = GhNew.vertices + Gh.GetId(v[1]);
		if (Adj[0]) Adj[0] =  GhNew.edges + Gh.GetId(Adj[0]);
		if (Adj[1]) Adj[1] =  GhNew.edges + Gh.GetId(Adj[1]);
	}
	/*}}}*/
	void   GeomEdge::SetTgA()         { /*{{{*/
		type |=4; /*=>2d digit to 1*/
	}/*}}}*/
	void   GeomEdge::SetTgB()         { /*{{{*/
		type |=8; /*=> 3d digit to 1*/
	}/*}}}*/
	void   GeomEdge::SetMark()        { /*{{{*/
		type |=16;/*=> 4th digiy to 1*/
	}/*}}}*/
	void   GeomEdge::SetUnMark()      { /*{{{*/
		type &= 1007 /* 1023-16 = 000111110111 => 4th digit to 0*/;
	}/*}}}*/
	void   GeomEdge::SetRequired()    { /*{{{*/
		type |= 64;/*=>6th digit to 1*/
	}/*}}}*/
	int    GeomEdge::TgA()     const  {/*{{{*/
		return type &4;
	}/*}}}*/
	int    GeomEdge::TgB()     const  {/*{{{*/
		return type &8;
	}/*}}}*/
}
