/*includes*/
#include "./DataSet.h"
//#include "../errorlib/include.h"
#include "./include.h"
using namespace std;

/*FUNCTION DataSet::DataSet(){{{*/
DataSet::DataSet(){
	return;
}
/*}}}*/
/*FUNCTION DataSet::~DataSet{{{*/
DataSet::~DataSet(){

	/*  use reverse_iterator for efficiency in matlab memory manager
		 (keeping old code in case it needs to revert back)  */

	//	vector<Object*>::iterator object;
	vector<Object*>::reverse_iterator object;

	//	for ( object=objects.begin() ; object < objects.end(); object++ ){
	//		delete (*object);
	//	}
	for ( object=objects.rbegin() ; object < objects.rend(); object++ ){
		delete (*object);
	}
}
/*}}}*/
/*FUNCTION DataSet::AddObject{{{*/
int  DataSet::AddObject(Object* obj){

	objects.push_back(obj);

	return 1;
}
/*}}}*/
/*FUNCTION DataSet::GetObjectByOffset{{{*/
Object* DataSet::GetObjectByOffset(int offset){

	/*Check index in debugging mode*/
	_assert_(this!=NULL);
	_assert_(offset<this->Size());

	return objects[offset];

}
/*}}}*/
/*FUNCTION DataSet::Size{{{*/
int  DataSet::Size(void){
	_assert_(this!=NULL);
	return objects.size();
}
/*}}}*/
