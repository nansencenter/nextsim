/*\file Print.h
 *\brief: print I/O for ISSM
 */

#ifndef _ISSM_PRINT_H_
#define _ISSM_PRINT_H_

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;
/*macros:*/
/* _printf_{{{*/
/* macro to print some string on all cpus */
#define _printf_(StreamArgs)\
  do{std::ostringstream aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy; \
	  aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy << StreamArgs; \
	  PrintfFunctionOnAllCpus(aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy.str());}while(0)
/*}}}*/
/* _printf0_ {{{*/
/* macro to print some string only on cpu 0 */
#define _printf0_(StreamArgs)\
  do{std::ostringstream aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy; \
	  aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy << StreamArgs; \
	  PrintfFunctionOnCpu0(aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy.str());}while(0)
/*}}}*/

/*functions: */
int PrintfFunctionOnCpu0(const string & message);
int PrintfFunctionOnAllCpus(const string & message);

#endif
