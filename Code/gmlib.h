#ifndef _GM_LIB_H_
#define _GM_LIB_H_

#include <stdio.h>
#include <stdlib.h>

/*************************
* Error handling for all of my code
*
*************************/

#define GMERR(ENUM) {GM_NoteError(__LINE__,ENUM);goto ERROR;}
#define GMERRH(FUNCNAME,RET) \
ERROR: \
	GM_PrintError(stderr,FUNCNAME,__FILE__); \
	return(RET);

#define GSVA(NAME,N) gsl_vector *NAME = gsl_vector_alloc(N)
#define GSMA(N1,N2) gsl_matrix_alloc(N1,N2)

extern int _gm_errNum;

void GM_NoteError(int lineNumber,int errCode);
void GM_PrintError(FILE *fp,const char* funcName,const char* fileName);

/* Declaration of a function pointer type */
typedef void (*GMFreeFunc_t)(void *p);


int 	GM_Initialize();
int 	GM_Finalize();

void *	GM_Malloc(size_t n);

int	GM_Free(void *p);
int	GM_Freeze(void *pp);
int	GM_FreezeList(char *message,...);
int	GM_FreezeListWithMethod(GMFreeFunc_t myFreeFunc,char *message, ...);

void GM_FreeGSLVector(void *gsl_vector);
void GM_FreeGSLMatrix(void *gsl_matrix);

#endif
