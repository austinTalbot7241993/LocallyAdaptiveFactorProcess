#include "gmlib.h"
#include <stdarg.h>
#include <gsl/gsl_matrix.h>


static int _isInitialized;
static int _LineNumber,_ErrorCode;

void GM_NoteError(int lineNumber,int errCode){
	_LineNumber = lineNumber;
	_ErrorCode = errCode;
}

void GM_PrintError(FILE *fp,const char* funcName,const char* fileName){
	fprintf(fp,"Error %d on line %d in function (%s) in file '%s'\n",\
	    _ErrorCode,_LineNumber,funcName,fileName);\

}

int     GM_Initialize()
{
	/* Call this at the beginning of your app */
	if ( ! _isInitialized ) {
		_isInitialized = 1;
	}
	return( 1 );
}

int     GM_Finalize()
{
	/* Call this at the end of your app */
	if ( _isInitialized ) {
	    _isInitialized = 0;
	}
	return( 1 );
}


void *  GM_Malloc(size_t n)
{
	return( (n > 0 ) ? malloc(n) : NULL );
}

int	GM_Free(void *p)
{
	if ( p == NULL )			goto ERROR;
	free(p);
	return( 0 );
ERROR:
	return( 1 );
}


int	GM_Freeze(void *pp)
{
	/* pp is really a void **p.  If the pointer is nozero it will
	   free it and set it to zero
	*/
    char *pc;
	char **ppc = (char **)pp;
	if ( ppc == NULL )			goto ERROR;
	pc = (char *)(*ppc);
	if ( pc != NULL ) { GM_Free(pc); }
	*ppc = NULL;
        return( 0 );
ERROR:
	return( 1 );
}


static int _freePointerList(int nitems, char **itemList)
{
	/* This takes a list of addresses to pointers and freezes it */
	int i;
	void *pp;
	/* Free them in the order they were given */
	for ( i = 0; i < nitems; i++ ) {
	    if ( (pp = itemList[i]) == NULL )			goto ERROR;
	    if ( GM_Freeze(pp) <= 0 ) {
		goto ERROR;
	    }
	}
	return( 0 );
ERROR:
	return( 1 );
}

int	GM_FreezeList(char *message,...)
{
    /* Give it a NULL-terminated list of ADDRESSES of pointers to stuff.
	This can only free things that were malloced by GM_Malloc.
    */
    void 		*item[128];
    int			maxItems = 127, nitems = 0;
    va_list		ap;
    void		*pv = NULL;
    /* Parse the subroutine argument list until we hit a NULL and gather
       up the list of addresses of pointers to be freed
    */
    va_start(ap,message);
    pv = va_arg(ap,void *);	
    while ( (pv != NULL) && (nitems < maxItems) ) {
	item[nitems++] = pv;
	 pv = va_arg(ap,void *);
    }
    va_end(ap);
    if ( _freePointerList(nitems, (char **) & item[0]) )	goto ERROR;
    return( 0 );
ERROR:
    return( 1 );
}


int     GM_FreezeListWithMethod(GMFreeFunc_t myFreeFunc,
				char *message, ... )
{
    /* Give it a NULL-terminated list of ADDRESSES of pointers to stuff.
        This can only free things that were malloced by GM_Malloc.
		myFreeFunc is presumably the pointer to a free function
		that takes a void *p as argument and returns void.
	For example:
		GM_FreeListWithMethod( &gsl_matrix_free,
					"x",
					& a,
					& b,
					& c,
					& d,
					NULL );
	ALL parameters passed are ADDRESES of variables to be freed.
	All must use one homogeneous function call to free.
    */
    void                *item[128];
    int                 maxItems = 127, nitems = 0;
    va_list             ap;
    GMFreeFunc_t	aFreeFunc;
    void                *pv = NULL;
    char		**ppc;
    int			i;

    /* Parse the subroutine argument list until we hit a NULL and gather
       up the list of addresses of pointers to be freed
    */
    va_start(ap,message);
//    aFreeFunc = va_arg(ap,GMFreeFunc_t);     /* This is the func */
    pv = va_arg(ap,void *);	
    while ( (pv != NULL) && (nitems < maxItems) ) {
        item[nitems++] = pv;
         pv = va_arg(ap,void *);
    }
    va_end(ap);
    if ( myFreeFunc == NULL ){
		fprintf(stderr,"aFreeFunc==NULL\n");
		goto ERROR;
	}
    for (i = 0; i < nitems; i++ ) {
	if ( (ppc = item[i]) != NULL ) {	
		if ( *ppc != NULL ) {
			(*myFreeFunc)( (void *)(*ppc) );
			*ppc = NULL;
		}
	}
    }

    return( 0 );
ERROR:
	fprintf(stderr,"Error in FreezeListWithMethod\n");
    return( 1 );
}

void GM_FreeGSLVector(void *vector){
	gsl_vector_free((gsl_vector *)vector);
}

void GM_FreeGSLMatrix(void *matrix){
	gsl_matrix_free((gsl_matrix *)matrix);
}

