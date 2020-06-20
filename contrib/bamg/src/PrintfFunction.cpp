/*\file PrintfFunction.c
 *\brief: this function is used by the _printf_ macro, to take into account the
 *fact we may be running on a cluster.
 */

#include <stdarg.h>
#include <cstdio>
#include "./Print.h"
//#include "../Comm/IssmComm.h"
#include "./sharedstring.h"
#include "./MemOps.h"

int PrintfFunctionOnCpu0(const string & message){
	int  my_rank;

	/*recover my_rank:*/
	my_rank=0;//IssmComm::GetRank();

	if(my_rank==0){
		#ifdef _HAVE_ANDROID_JNI_
		__android_log_print(ANDROID_LOG_INFO, "Native",message.c_str());
		#else
		ApiPrintf(message.c_str());
		#endif
	}
	return 1;
}
int PrintfFunctionOnAllCpus(const string & message){

	#ifdef _HAVE_ANDROID_JNI_
	__android_log_print(ANDROID_LOG_INFO, "Native",message.c_str());
	#else
	ApiPrintf(message.c_str());
	#endif

	return 1;
}
