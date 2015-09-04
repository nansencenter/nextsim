/*!\file:  ApiPrintf: this file only gets compiled and linked for issm.exe. API specific
 * versions of this file (for Python and Matlab) will be used for the wrappers (modules).
 */

#include <stdio.h>

void ApiPrintf(const char* string){
	printf("%s",string);
}
