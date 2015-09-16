#ifndef _DATASET_H_
#define _DATASET_H_

#include <vector>
#include "./Object.h"

/* class Object{ */
/* 	public:  */
/* 		virtual       ~Object() {}; */
/* }; */
class DataSet{

	public:
	    std::vector<Object*> objects;
		DataSet();
		DataSet(double* observations_list,double* x,double* y,int n);
		~DataSet();
		int AddObject(Object* obj);
		int Size();
		Object* GetObjectByOffset(int offset);
};

#endif
