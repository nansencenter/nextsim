/*! \file Option.h
 *  \brief: header file for option abstract object
 */

#ifndef _OPTIONOBJECT_H_
#define _OPTIONOBJECT_H_

/*Headers:{{{*/
#include "./shared.h"
#include "./Object.h"
#include "./exceptions.h"
/*}}}*/

class Option: public Object {

	public:

		/*Option constructors, destructors*/
		Option(){};
		~Option(){};

		/*Object virtual functions definitions*/
		virtual void  Echo()= 0;
		virtual void  DeepEcho()= 0;
		virtual void  DeepEcho(char  *indent)=0;
		int           Id(){_error_("Not implemented yet"); };
		int           ObjectEnum(){return OptionEnum;              };
		Object       *copy(){_error_("Not implemented yet"); };

		/*virtual functions: */
		virtual char* Name()=0;
		virtual int   NumEl()=0;
		virtual int   NDims()=0;
		virtual int*  Size()=0;

};
#endif  /* _OPTIONOBJECT_H */
