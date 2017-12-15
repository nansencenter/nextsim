#include "./include.h"
#include "./shared.h"

/*Constructors/Destructors*/
BamgOpts::BamgOpts(){/*{{{*/

	this->anisomax=0;
	this->cutoff=0;
	this->coeff=0;
	this->errg=0;
	this->gradation=0;
	this->Hessiantype=0;
	this->MaxCornerAngle=0;
	this->maxnbv=0;
	this->maxsubdiv=0;
	this->Metrictype=0;
	this->nbjacobi=0;
	this->nbsmooth=0;
	this->omega=0;
	this->power=0;
	this->random=0;
	this->verbose=0;

	this->geometricalmetric=0;
	this->KeepVertices=0;
	this->splitcorners=0;

	this->hmin=0;
	this->hmax=0;
	this->hminVertices=NULL; this->hminVerticesSize[0]=this->hminVerticesSize[1]=0;
	this->hmaxVertices=NULL; this->hmaxVerticesSize[0]=this->hmaxVerticesSize[1]=0;
	this->hVertices=NULL;    this->hVerticesSize[0]=this->hVerticesSize[1]=0;
	this->metric=NULL;       this->metricSize[0]=this->metricSize[1]=0;
	this->field=NULL;        this->fieldSize[0]=this->fieldSize[1]=0;
	this->err=NULL;          this->errSize[0]=this->errSize[1]=0;

}
/*}}}*/
BamgOpts::~BamgOpts(){/*{{{*/

	delete [] this->hminVertices;
	delete [] this->hmaxVertices;
	delete [] this->hVertices;
	delete [] this->metric;
	delete [] this->field;
	delete [] this->err;

}
/*}}}*/

/*Methods*/
void BamgOpts::Check(void){/*{{{*/

	int i;

	if (this->anisomax<1) _error_("'anisomax' option should be >=1");
	if (this->coeff==0) _error_("'coeff' should be positive");
	if (this->errg<0) _error_("'errg' option should be >0");
	if (this->gradation<1) _error_("'gradation' option should be >=1");
	if (this->Hessiantype!=0  && this->Hessiantype!=1) _error_("'Hessiantype' supported options are 0 and 1");
	if (this->maxnbv<3) _error_("'maxnbv' option should be >3");
	if (this->maxsubdiv<=1) _error_("'maxsubdiv' should be >1");
	if (this->Metrictype!=0   && this->Metrictype!=1 && this->Metrictype!=2) _error_("'Metrictype' supported options are 0, 1 and 2");
	if (this->nbjacobi<=0) _error_("'nbjacobi' option should be >0");
	if (this->nbsmooth<=0) _error_("'nbsmooth' option should be >0");

	if (this->KeepVertices!=0 && this->KeepVertices!=1) _error_("'KeepVertices' supported options are 0 and 1");
	if (this->geometricalmetric!=0  && this->geometricalmetric!=1) _error_("'geometricalmetric' supported options are 0 and 1");

	if (this->hmin<=0) _error_("'hmin' option should be >0");
	if (this->hmax<=0 || this->hmax<this->hmin) _error_("'hmax' option should be between 0 and hmin=" << this->hmin);
	if (this->hminVertices && this->hminVerticesSize[1]!=1) _error_("'hminVertices' should be a column");
	if (this->hmaxVertices && this->hmaxVerticesSize[1]!=1) _error_("'hmaxVertices' should be a column");
	if (this->hVertices && this->hVerticesSize[1]!=1) _error_("'hVertices' should be a column");
	if (this->metric && (this->metricSize[1]!=1 && this->metricSize[1]!=3)) _error_("'metric' should have either 1 (iso) or 3 (aniso) columns.");
	if (this->field){
		if (this->errSize[0]!=1 || this->errSize[1]!=this->fieldSize[1]) _error_("'err' should be of size " << 1 << " x " << this->fieldSize[1]);
		for (i=0;i<this->fieldSize[1];i++) {if (this->err[i]<=0) _error_("'err' option should be >0");};
	}

}
/*}}}*/
