/*!\file exceptions.h
 * \brief: two types of exceptions are handled for now. Errors, and
 * warnings. Those exceptions are trapped provided the matlab modules
 * are started using MODULEBOOT, and ended using MODULEEND. These are
 * macros hiding try, catch statements. This header file defines our
 * own exceptions
 */

#ifndef _MY_EXCEPTIONS_H_
#define _MY_EXCEPTIONS_H_

#include <exception>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

/* #ifdef HAVE_CONFIG_H */
/* 	#include <config.h> */
/* #else */
/* #error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!" */
/* #endif */

/*macros: */
/* _assert_ {{{*/
/*Assertion macro: do nothing if macro _ISSM_DEBUG_ undefined*/
#ifdef _ISSM_DEBUG_
#define _assert_(statement)\
  if (!(statement)) _error_("Assertion \""<<#statement<<"\" failed, please report bug to "<<PACKAGE_BUGREPORT)
#else
#define _assert_(ignore)\
  ((void) 0)
#endif
/*}}}*/
/* _error_ {{{*/
/*new Error exception macro*/
#ifdef _INTEL_WIN_
#define _error_(StreamArgs)\
   do{std::ostringstream aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy; \
   aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy << StreamArgs << std::ends; \
   throw ErrorException(aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy.str());}while(0)
#else
#define _error_(StreamArgs)\
	do{std::ostringstream aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy; \
   aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy << StreamArgs << std::ends; \
   throw ErrorException(__FILE__,__func__,__LINE__,aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy.str());}while(0)
#endif
/*}}}*/
/* ExceptionTrapBegin/ExceptionTrapEnd {{{*/

/*The following macros hide the error exception handling in a matlab module. Just put
 * ExceptionTrapBegin(); and ExceptionTrapEnd(); at the beginning and end of a module, and c++ exceptions
 * will be trapped. Really nifty!*/

#define ExceptionTrapBegin(); \
	try{

#define ExceptionTrapEnd(); }\
	catch(ErrorException &exception){\
		exception.Report();\
		return 0;\
	}\
	catch(exception& e) {\
		_printf_("Standard exception: " << e.what() << "\n\n");\
		return 0;\
	}\
	catch(...){\
		_printf_("An unexpected error occurred \n\n");\
		return 0;\
	}
/*}}}*/

/*ISSM exception class: */
class ErrorException: public exception { /*{{{*/

	char* what_str;
	char* function_name;
	char* file_name;
	int   file_line;

	public:
	ErrorException(const string &what_arg); //for windows
	ErrorException(const string &what_file,const string& what_function,int what_line,const string& what_arg);//UNIX
	~ErrorException() throw();
	virtual const char *what() const throw();
	void Report() const;
	const char* MatlabReport() const;
	const char* PythonReport() const;

};
/*}}}*/
#endif
