#ifndef _DATASET_H_
#define _DATASET_H_

#include <vector>

class Object{
	public: 
		virtual       ~Object() {};
};
class DataSet{
	private:
		std::vector<Object*> objects;
	public:
		DataSet();
		DataSet(double* observations_list,double* x,double* y,int n);
		~DataSet();
		int AddObject(Object* obj);
		int Size();
		Object* GetObjectByOffset(int offset);
};

#endif
