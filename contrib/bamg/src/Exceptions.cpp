/* \file Exceptions.cpp
 * \brief: implementation of the exceptions.
 */

#include <cstring>
#include <cstdio>
#include "./exceptions.h"
#include "./Print.h"
//#include "./IssmComm.h"
#include "./MemOps.h"

ErrorException::ErrorException(const string & what_arg){/*{{{*/

	int len;
	len           = strlen(what_arg.c_str())+1;
	what_str      = new char[len];
	memcpy(what_str,what_arg.c_str(),len);

	file_name     = NULL;
	function_name = NULL;
	file_line     = 0;

}/*}}}*/
ErrorException::ErrorException(const string& what_file, const string& what_function,int what_line, const string& what_arg){/*{{{*/

	int len;

	len      = strlen(what_arg.c_str())+1;
	what_str = new char[len];
	memcpy(what_str,what_arg.c_str(),len);

	len       = strlen(what_file.c_str())+1;
	file_name = new char[len];
	memcpy(file_name,what_file.c_str(),len);

	len           = strlen(what_function.c_str())+1;
	function_name = new char[len];
	memcpy(function_name,what_function.c_str(),len);

	file_line= what_line;
	/*When error messages are not shown properly, uncomment the following line*/
	//this->Report();

}/*}}}*/
ErrorException::~ErrorException() throw(){/*{{{*/
	delete [] what_str;
	delete [] file_name;
	delete [] function_name;
}/*}}}*/
const char* ErrorException::what() const throw(){/*{{{*/
	return what_str;
}/*}}}*/
void ErrorException::Report() const{/*{{{*/

	/*WINDOWS*/
	if(!function_name || file_line==0){
		_printf_("Error message: " << what());
		return;
	}

	/*recover my_rank and num_procs:*/
	//int my_rank   = IssmComm::GetRank();
	//int num_procs = IssmComm::GetSize();
	int my_rank = 1;
	int num_procs = 1;

	if(num_procs==1){
		_printf_("\n??? Error in ==> " << file_name << ":" << file_line << "\n");
		_printf_(function_name << " error message: " << what() << "\n\n");
	}
	else{
		_printf_("\n[" << my_rank<< "] ??? Error using ==> " << file_name << ":" << file_line << "\n");
		_printf_(  "[" << my_rank << "] " << function_name << " error message: " << what() << "\n\n");
	}

	return;
}/*}}}*/
const char* ErrorException::MatlabReport() const{/*{{{*/

	/*Output*/
	std::ostringstream buffer;
	char *message = NULL;

	/*WINDOWS*/
	if(!function_name || file_line==0){
		buffer << " error message: " << this->what_str;
	}
	else{
		buffer << "\nError in ==> " << this->file_name << ":" << file_line << "\n";
		buffer << this->function_name << " error message: " << this->what_str;
	}

	/*Convert std::ostringstream to std::string and then create char* */
	std::string buffer2 = buffer.str();
	message = xNew<char>(strlen(buffer2.c_str())+1); sprintf(message,"%s",buffer2.c_str());
	return message;
}/*}}}*/
const char* ErrorException::PythonReport() const{/*{{{*/

	/*Output*/
	std::ostringstream buffer;
	char *message = NULL;

	/*WINDOWS*/
	if(!function_name || file_line==0){
		buffer << " error message: " << this->what_str;
	}
	else{
		buffer << "\nError in ==> " << this->file_name << ":" << file_line << "\n";
		buffer << this->function_name << " error message: " << this->what_str;
	}

	/*Convert std::ostringstream to std::string and then create char* */
	std::string buffer2 = buffer.str();
	message = xNew<char>(strlen(buffer2.c_str())+1); sprintf(message,"%s",buffer2.c_str());
	return message;
}/*}}}*/
