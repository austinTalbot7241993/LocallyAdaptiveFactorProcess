#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include "gmlib.h"
#include "matrix_utils.h"
#include "mdarray_complex.h"





int marray3d_complex_checkSizes(marray3d_complex *self,size_t d1,size_t d2,size_t d3){
	return((self==NULL)||(self->d1!=d1)||(self->d2!=d2)||(self->d3=d3));
}

int marray3d_complex_checkX1orT(marray3d_complex *self,size_t t){
	return((self==NULL)||((self->d1!=1)&&(self->d1!=t)));
}

int marray3d_complex_checkY1orT(marray3d_complex *self,size_t t){
	return((self==NULL)||((self->d2!=1)&&(self->d2!=t)));
}

int marray3d_complex_checkZ1orT(marray3d_complex *self,size_t t){
	return((self==NULL)||((self->d3!=1)&&(self->d3!=t)));
}

int marray3d_complex_checkXY(marray3d_complex *self,size_t d1,size_t d2){
	return((self==NULL)||(self->d1!=d1)||(self->d2!=d2));
}

int marray3d_complex_checkYZ(marray3d_complex *self,size_t d2,size_t d3){
	return((self==NULL)||(self->d3!=d3)||(self->d2!=d2));
}

int marray3d_complex_checkXZ(marray3d_complex *self,size_t d1,size_t d3){
	return((self==NULL)||(self->d1!=d1)||(self->d3!=d3));
}

int marray3d_complex_compareSizes(marray3d_complex *m1,marray3d_complex *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d1!=m2->d1)||(m1->d2!=m2->d2)||(m1->d3!=m2->d3));
}

int marray3d_complex_compareXY(marray3d_complex *m1,marray3d *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d1!=m2->d1)||(m1->d2!=m2->d2));
}

int marray3d_complex_compareYZ(marray3d_complex *m1,marray3d *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d2!=m2->d2)||(m1->d3!=m2->d3));
}

int marray3d_complex_compareXZ(marray3d_complex *m1,marray3d *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d1!=m2->d1)||(m1->d3!=m2->d3));
}

//Allocation and free methods
marray3d_complex * marray3d_alloc(size_t d1,size_t d2,size_t d3){
	int i;
	marray3d_complex *self = (marray3d_complex *)(GM_Malloc(sizeof(marray3d_complex)));
	if(self==NULL)					GMERR(-1);
	self->d1 = d1;
	self->d2 = d2;
	self->d3 = d3;
	self->matrixList = (gsl_matrix **)GM_Malloc(d1*sizeof(gsl_matrix_complex *));
	for(i=0;i<d1;i++){
		self->matrixList[i] = gsl_matrix_complex_alloc(d2,d3);
		if(self->matrixList[i]==NULL)		GMERR(-2);
	}
	return(self);
GMERRH("marray3d_alloc",NULL);
}

int marray3d_complex_free(marray3d_complex *self){
	int i;
	for(i=0;i<self->d1;i++){
		gsl_matrix_complex_free(self->matrixList[i]);
		self->matrixList[i] = NULL;
	}
	GM_Free(self->matrixList);
	self->matrixList = NULL;
	GM_Free(self);
	return(0);
GMERRH("marray3d_complex_free",1);
}

//Initialization methods

gsl_complex marray3d_complex_get(marray3d_complex *self,size_t x,size_t y,size_t z){
	return(temp);
}

int marray3d_complex_set(marray3d *self,size_t x,size_t y,size_t z,gsl_complex val){
	return(0);
GMERRH("marray3d_complex_set",1);
}

int marray3d_complex_set_zero(marray3d_complex *self){
	int i,d1;
	if(self==NULL)						GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_complex_set_zero(self->matrixList[i]);
	}
	return(0);
GMERRH("marray3d_complex_set_zero",1);
}

int marray3d_complex_set_constant(marray3d_complex *self,gsl_complex x){
	int i,d1;
	if(self==NULL)						GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_complex_set_all(self->matrixList[i],x);
	}
	return(0);
GMERRH("marray3d_complex_set_constant",1);
}

int marray3d_complex_set_identity(marray3d_complex *self){
	int i,d1;
	if(self==NULL)						GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_complex_set_identity(self->matrixList[i]);
	}
	return(0);
GMMERH("marray3d_complex_set_identity",1);	
}

int marray3d_complex_set_matrix(marray3d_complex *self,gsl_matrix_complex *m){
	int i,d1;
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((self->d2!=m->size1)||(self->d3!=m->size2))	GMERR(-11);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_complex_memcpy(self->matrixList[i],m);
	}
	return(0);
GMERRH("marray3d_complex_set_matrix",1);
}

//Read and write methods
int marray3d_complex_read(FILE *fp,marray3d_complex *self){
    int d1,d2,d3,i,j,k,bufflen,maxbufflen=1024;
	gsl_complex temp = gsl_complex_rect(0,0);
    int size1=0,size2=0,size3=0;
    double val=0;
    char buff[1024];
    if((fp==NULL)||(self==NULL))                    GMERR(-1);
    d1 = self->d1;
    d2 = self->d2;
    d3 = self->d3;
    if(fgets(buff,1023,fp)==NULL)       GMERR(-11);
    if(buff[0]!='#')                    GMERR(-21);
    sscanf(buff+2,"%d %d %d",&size1,&size2,&size3);
    if(size1!=d1||size2!=d2||size3!=d3) GMERR(-31);
    for(i=0;i<d1;i++){
        for(j=0;j<d2;j++){
            for(k=0;k<d3;k++){
                fscanf(fp,"%lf",&val);
				GSL_SET_REAL(&temp,val);
                fscanf(fp,"%lf",&val);
				GSL_SET_IMAG(&temp,val);
                marray3d_complex_set(self,i,j,k,temp);
            }    
        }    
    }    
    return(0);
GMERRH("marray3d_complex_read",1);
}

int marray3d_complex_write(FILE *fp,marray3d_complex *self){
    int d1,d2,d3,i,j,k;
	gsl_complex temp;
    if((fp==NULL)||(self==NULL))                    GMERR(-1);
    d1 = self->d1;
    d2 = self->d2;
    d3 = self->d3;
    fprintf(fp,"# %d %d %d\n",d1,d2,d3);
    for(i=0;i<d1;i++){
        for(j=0;j<d2;j++){
            for(k=0;k<d3;k++){
				temp = marray3d_complex_get(self,i,j,k);
                fprintf(fp,"%0.6f\n",GSL_REAL(temp));
                fprintf(fp,"%0.6f\n",GSL_IMAG(temp));
            }
        }
    }
    return(0);
GMERRH("marray3d_complex_write",1);
}

int marray3d_complex_print(marray3d_complex *self){
	GMERR(-9999);
	return(0);
GMERRH("marray3d_complex_print",1);
}

int marray3d_complex_fancyWrite(FILE *fp,marray3d_complex *self){
	GMERR(-9999);
	return(0);
GMERRH("marray3d_complex_fancyWrite",1);
}

//Slices
int marray3d_complex_set_X(marray3d_complex *self,gsl_matrix_complex *m,int x){
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((x<0)||(self->d1<=x))						GMERR(-11);
	if((self->d1!=m->size1)||(self->d3!=m->size2))	GMERR(-21);
	gsl_matrix_complex_memcpy(self->matrixList[x],m);
	return(0);
GMERRH("marray3d_complex_set_X",1);
}

int marray3d_complex_set_Y(marray3d_complex *self,gsl_matrix_complex *m,int y){
    int i,j,d1,d3;
    gsl_complex temp;
    if((self==NULL)||(m==NULL))                     GMERR(-1);
    if((y<0)||(self->d2<=y))                        GMERR(-11);
    if((self->d1!=m->size1)||(self->d3!=m->size2))  GMERR(-21);
    d1 = self->d1;
    d3 = self->d3;
    for(i=0;i<d1;i++){
        for(j=0;j<d3;j++){
            temp = gsl_matrix_complex_get(m,i,j);
            gsl_matrix_complex_set(self->matrixList[i],y,j,temp);
        }
    }
    return(0);
GMERRH("marray3d_complex_set_Y",1);
}

int marray3d_complex_set_Z(marray3d_complex *self,gsl_matrix_complex *m,int z){
    int i,j,d1,d2;
    gsl_complex temp;
    if((self==NULL)||(m==NULL))                     GMERR(-1);
    if((z<0)||(self->d3<=z))                        GMERR(-11);
    if((self->d1!=m->size1)||(self->d2!=m->size2))  GMERR(-21);
    d1 = self->d1;
    d2 = self->d2;
    for(i=0;i<d1;i++){
        for(j=0;j<d2;j++){
            temp = gsl_matrix_complex_get(m,i,j);
            gsl_matrix_complex_set(self->matrixList[i],j,z,temp);
        }
    }
    return(0);
GMERRH("marray3d_set_Z",1);
}

int marray3d_complex_get_X(marray3d_complex *self,gsl_matrix_complex *m,int x){
    if((self==NULL)||(m==NULL))                     GMERR(-1);
    if((x<0)||(self->d1<=x))                        GMERR(-11);
    if((self->d2!=m->size1)||(self->d3!=m->size2))  GMERR(-21);
    gsl_matrix_complex_memcpy(m,self->matrixList[x]);
    return(0);
GMERRH("marray3d_complex_get_X",1);
}

int marray3d_complex_get_Y(marray3d_complex *self,gsl_matrix_complex *m,int y){
    int i,j,d1,d3;
    gsl_complex temp;
    if((self==NULL)||(m==NULL))                     GMERR(-1);
    if((y<0)||(self->d2<=y))                        GMERR(-11);
    if((self->d1!=m->size1)||(self->d3!=m->size2))  GMERR(-21);
    d1 = self->d1;
    d3 = self->d3;
    for(i=0;i<d1;i++){
        for(j=0;j<d3;j++){
            temp = gsl_matrix_get(self->matrixList[i],y,j);
            gsl_matrix_set(m,i,j,temp);
        }
    }

    return(0);
GMERRH("marray3d_complex_get_Y",1);
}

int marray3d_complex_get_Z(marray3d_complex *self,gsl_matrix_complex *m,int z){
    int i,j,d1,d2;
    gsl_complex temp;
    if((self==NULL)||(m==NULL))                     GMERR(-1);
    if((z<0)||(self->d3<=z))                        GMERR(-11);
    if((self->d1!=m->size1)||(self->d2!=m->size2))  GMERR(-21);
    d1 = self->d1;
    d2 = self->d2;
    for(i=0;i<d1;i++){
        for(j=0;j<d2;j++){
            temp = gsl_matrix_get(self->matrixList[i],j,z);
            gsl_matrix_complex_set(m,i,j,temp);
        }
    }

    return(0);
GMERRH("marray3d_complex_get_Z",1);
}

int marray3d_complex_set_pencilX(marray3d_complex *self,gsl_vector_complex *v,int y,int z){
	if((self==NULL)||(v==NULL))				GMERR(-1);
	return(0);
GMERRH("",1);
}

int marray3d_complex_set_pencilY(marray3d_complex *self,gsl_vector_complex *v,int x,int z){
	if((self==NULL)||(v==NULL))				GMERR(-1);
	return(0);
GMERRH("",1);
}

int marray3d_complex_set_pencilZ(marray3d_complex *self,gsl_vector_complex *v,int x,int y){
	if((self==NULL)||(v==NULL))				GMERR(-1);
	return(0);
GMERRH("",1);
}

int marray3d_complex_get_pencilX(marray3d_complex *self,gsl_vector_complex *v,int y,int z){
	if((self==NULL)||(v==NULL))				GMERR(-1);
	return(0);
GMERRH("",1);
}

int marray3d_complex_get_pencilY(marray3d_complex *self,gsl_vector *v,int x,int z){
	if((self==NULL)||(v==NULL))				GMERR(-1);
	return(0);
GMERRH("",1);
}

int marray3d_complex_get_pencilZ(marray3d_complex *self,gsl_vector *v,int x,int y){
	if((self==NULL)||(v==NULL))				GMERR(-1);
	return(0);
GMERRH("",1);
}


int marray3d_complex_get_matrixSliceX(marray3d_complex *self,gsl_matrix_complex *m,int x,int y1,int y2,gsl_vector_complex *work);
int marray3d_complex_set_matrixSliceX(marray3d_complex *self,gsl_matrix_complex *m,int x,int y1,int y2,gsl_vector_complex *work);












