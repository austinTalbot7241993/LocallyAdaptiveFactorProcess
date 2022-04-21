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
#include "mdarray.h"

/*Code by Austin Bryan Talbot
 *        a.talbot@verizon.net
 *
 *Description:
 *
 *Completed:
 *
 *Revision History:
 *
*/

/************************************
*************************************
**								   **
** Methods for 3 dimensional array **
**								   **
*************************************
************************************/


/************************
* Dimension comparisons *
************************/

//works
int marray3d_checkSizes(marray3d *self,size_t d1,size_t d2,size_t d3){
	return((self==NULL)||(self->d1!=d1)||(self->d2!=d2)||(self->d3!=d3));
}

int marray3d_checkX1orT(marray3d *self,size_t t){
	return((self==NULL)||((self->d1!=1)&&(self->d1!=t)));
}

int marray3d_checkY1orT(marray3d *self,size_t t){
	return((self==NULL)||((self->d2!=1)&&(self->d2!=t)));
}

int marray3d_checkZ1orT(marray3d *self,size_t t){
	return((self==NULL)||((self->d3!=1)&&(self->d3!=t)));
}

int marray3d_checkXY(marray3d *self,size_t d1,size_t d2){
	return((self==NULL)||(self->d1!=d1)||(self->d2!=d2));
}

int marray3d_checkXZ(marray3d *self,size_t d1,size_t d3){
	return((self==NULL)||(self->d1!=d1)||(self->d3!=d3));
}

int marray3d_checkYZ(marray3d *self,size_t d2,size_t d3){
	return((self==NULL)||(self->d2!=d2)||(self->d3!=d3));
}

//works
int marray3d_compareSizes(marray3d *m1,marray3d *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d1!=m2->d1)||(m1->d2!=m2->d2)||(m1->d3!=m2->d3));
}

int marray3d_compareXY(marray3d *m1,marray3d *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d1!=m2->d1)||(m1->d2!=m2->d2));
}

int marray3d_compareYZ(marray3d *m1,marray3d *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d2!=m2->d2)||(m1->d3!=m2->d3));
}

int marray3d_compareXZ(marray3d *m1,marray3d *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d1!=m2->d1)||(m1->d3!=m2->d3));
}


/******************************
* Allocation and free methods *
******************************/

//works
marray3d * marray3d_alloc(size_t d1,size_t d2,size_t d3){
	int i;
	marray3d *self = (marray3d *)(GM_Malloc(sizeof(marray3d)));
	if(self==NULL)					GMERR(-1);
	self->d1 = d1;
	self->d2 = d2;
	self->d3 = d3;
	self->matrixList = (gsl_matrix **)GM_Malloc(d1*sizeof(gsl_matrix *));
	for(i=0;i<d1;i++){
		self->matrixList[i] = gsl_matrix_alloc(d2,d3);
		if(self->matrixList[i]==NULL)	GMERR(-2);
	}
	return(self);
GMERRH("marray3d_alloc",NULL);
}

//works
int marray3d_free(marray3d *self){
	int i;
	for(i=0;i<self->d1;i++){
		gsl_matrix_free(self->matrixList[i]);
		self->matrixList[i] = NULL;
	}
	GM_Free(self->matrixList);
	self->matrixList = NULL;
	GM_Free(self);
	return(0);
GMERRH("marray3d_free",1);
}




/*************************
* Initialization methods *
*************************/

//works
int marray3d_set_zero(marray3d *self){
	int i,d1;
	if(self==NULL)						GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_set_zero(self->matrixList[i]);
	}
	return(0);
GMERRH("marray3d_set_zero",1);
}

//works
int marray3d_set_constant(marray3d *self,double x){
	int i,d1;
	if(self==NULL)						GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_set_all(self->matrixList[i],x);
	}
	return(0);
GMERRH("marray3d_set_constant",1);
}

//works
int marray3d_set_identity(marray3d *self){
	int i,d1;
	if(self==NULL)						GMERR(-1);
	if(self->d2!=self->d3)				GMERR(-11);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_set_identity(self->matrixList[i]);
	}
	return(0);
GMERRH("marray3d_set_identity",1);
}

//works
int marray3d_set_matrix(marray3d *self,gsl_matrix *m){
	int i,d1;
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((self->d2!=m->size1)||(self->d3!=m->size2))	GMERR(-11);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_memcpy(self->matrixList[i],m);
	}
	return(0);
GMERRH("marray3d_set_matrix",1);
}

/*************************
* Read and write methods *
*************************/

int marray3d_read(FILE *fp,marray3d *self){
	int d1,d2,d3,i,j,k,bufflen,maxbufflen=1024;
	int size1=0,size2=0,size3=0;
	double val=0;
	char buff[1024];
	if((fp==NULL)||(self==NULL))					GMERR(-1);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	if(fgets(buff,1023,fp)==NULL)		GMERR(-11);
	if(buff[0]!='#')					GMERR(-21);
	sscanf(buff+2,"%d %d %d",&size1,&size2,&size3);
	if(size1!=d1||size2!=d2||size3!=d3)	GMERR(-31);
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				fscanf(fp,"%lf",&val);
				marray3d_set(self,i,j,k,val);
			}
		}
	}
	return(0);
GMERRH("marray3d_read",1);
}

int marray3d_write(FILE *fp,marray3d *self){
	int d1,d2,d3,i,j,k;
	if((fp==NULL)||(self==NULL))					GMERR(-1);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	fprintf(fp,"# %d %d %d\n",d1,d2,d3);
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				fprintf(fp,"%0.6f\n",marray3d_get(self,i,j,k));
			}
		}
	}

	return(0);
GMERRH("marray3d_write",1);
}

int marray3d_fancyWrite(FILE *fp,marray3d *self){
	int i,j,k;
	int d1,d2,d3;
	if((fp==NULL)||(self==NULL))					GMERR(-1);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				fprintf(fp,"%0.15f ",marray3d_get(self,i,j,k));
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	return(0);
GMERRH("marray3d_fancyWrite",1);
}

//works
int marray3d_print(marray3d *self){
	int d1,i;
	if(self==NULL)									GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		printGSLMatrix(self->matrixList[i]);
		printf("\n");
	}
	return(0);
GMERRH("marray3d_print",1);
}

/******************************************************
* Random initialization with a common random variable *
******************************************************/

int marray3d_set_randN(marray3d *self,gsl_rng *r,double mu,double sigma){
	int d1,d2,d3,i,j,k;
	double temp;
	if(self==NULL)									GMERR(-1);
	if(sigma<=0)									GMERR(-11);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				temp = gsl_ran_gaussian(r,sigma) + mu;
				gsl_matrix_set(self->matrixList[i],j,k,temp);
			}
		}
	}
	return(0);
GMERRH("marray3d_set_randN",1);
}

int marray3d_set_randG(marray3d *self,gsl_rng *r,double rate,double scale){
	int d1,d2,d3,i,j,k;
	double temp;
	if(self==NULL)									GMERR(-1);
	if((rate<=0)||(scale<=0))						GMERR(-11);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				temp = gsl_ran_gamma(r,rate,scale);
				gsl_matrix_set(self->matrixList[i],j,k,temp);
			}
		}
	}
	return(0);
GMERRH("marray3d_set_randG",1);
}

int marray3d_set_randIG(marray3d *self,gsl_rng *r,double rate,double scale){
	int d1,d2,d3,i,j,k;
	double temp;
	if(self==NULL)									GMERR(-1);
	if((rate<=0)||(scale<=0))						GMERR(-11);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				temp = 1.0/gsl_ran_gamma(r,rate,1.0/scale);
				gsl_matrix_set(self->matrixList[i],j,k,temp);
			}
		}
	}
	return(0);
GMERRH("marray3d_set_randIG",1);
}

int marray3d_set_randU(marray3d *self,gsl_rng *r,double a,double b){
	int d1,d2,d3,i,j,k;
	double temp;
	if(self==NULL)									GMERR(-1);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				temp = gsl_ran_flat(r,a,b);
				gsl_matrix_set(self->matrixList[i],j,k,temp);
			}
		}
	}
	return(0);
GMERRH("marray3d_set_randU",1);
}


//Random initializaiton with a function
//int marray3d_set_function(marray3d *self,fcnptr)

/**********************
* Setters and getters *
**********************/

//Individual points

//works
double marray3d_get(marray3d *self,size_t x,size_t y,size_t z){
	return(gsl_matrix_get(self->matrixList[x],y,z));
}

//works
int marray3d_set(marray3d *self,size_t x,size_t y,size_t z,double val){
	if(self==NULL)									GMERR(-1);
	gsl_matrix_set(self->matrixList[x],y,z,val);
	return(0);
GMERRH("marray3d_set",1);
}


//Slices
int marray3d_set_X(marray3d *self,gsl_matrix *m,int x){
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((x<0)||(self->d1<=x))						GMERR(-11);
	if((self->d2!=m->size1)||(self->d3!=m->size2))	GMERR(-21);
	gsl_matrix_memcpy(self->matrixList[x],m);
	return(0);
GMERRH("marray3d_set_X",1);
}

int marray3d_set_Y(marray3d *self,gsl_matrix *m,int y){
	int i,j,d1,d3;
	double temp;
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((y<0)||(self->d2<=y))						GMERR(-11);
	if((self->d1!=m->size1)||(self->d3!=m->size2))	GMERR(-21);
	d1 = self->d1;
	d3 = self->d3;
	for(i=0;i<d1;i++){
		for(j=0;j<d3;j++){
			temp = gsl_matrix_get(m,i,j);
			gsl_matrix_set(self->matrixList[i],y,j,temp);
		}
	}

	return(0);
GMERRH("marray3d_set_Y",1);
}

int marray3d_set_Z(marray3d *self,gsl_matrix *m,int z){
	int i,j,d1,d2;
	double temp;
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((z<0)||(self->d3<=z))						GMERR(-11);
	if((self->d1!=m->size1)||(self->d2!=m->size2))	GMERR(-21);
	d1 = self->d1;
	d2 = self->d2;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			temp = gsl_matrix_get(m,i,j);
			gsl_matrix_set(self->matrixList[i],j,z,temp);
		}
	}

	return(0);
GMERRH("marray3d_set_Z",1);
}

int marray3d_get_X(marray3d *self,gsl_matrix *m,int x){
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((x<0)||(self->d1<=x))						GMERR(-11);
	if((self->d2!=m->size1)||(self->d3!=m->size2))	GMERR(-21);
	gsl_matrix_memcpy(m,self->matrixList[x]);
	return(0);
GMERRH("marray3d_get_X",1);
}

int marray3d_get_Y(marray3d *self,gsl_matrix *m,int y){
	int i,j,d1,d3;
	double temp;
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((y<0)||(self->d2<=y))						GMERR(-11);
	if((self->d1!=m->size1)||(self->d3!=m->size2))	GMERR(-21);
	d1 = self->d1;
	d3 = self->d3;
	for(i=0;i<d1;i++){
		for(j=0;j<d3;j++){
			temp = gsl_matrix_get(self->matrixList[i],y,j);
			gsl_matrix_set(m,i,j,temp);
		}
	}

	return(0);
GMERRH("marray3d_get_Y",1);
}

int marray3d_get_Z(marray3d *self,gsl_matrix *m,int z){
	int i,j,d1,d2;
	double temp;
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((z<0)||(self->d3<=z))						GMERR(-11);
	if((self->d1!=m->size1)||(self->d2!=m->size2))	GMERR(-21);
	d1 = self->d1;
	d2 = self->d2;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			temp = gsl_matrix_get(self->matrixList[i],j,z);
			gsl_matrix_set(m,i,j,temp);
		}
	}

	return(0);
GMERRH("marray3d_get_Z",1);
}


int marray3d_get_pencilX(marray3d *self,gsl_vector *v,int y,int z){
	int i,dim1;
	GMERR(-1111);
	if((self==NULL)||(v==NULL))						GMERR(-1);
	if(v->size!=self->d1)							GMERR(-11);
	dim1=self->d1;
	for(i=0;i<dim1;i++){
		gsl_vector_set(v,i,gsl_matrix_get(self->matrixList[i],y,z));
	}

	return(0);
GMERRH("marray3d_get_pencilX",1);
}

int marray3d_get_pencilY(marray3d *self,gsl_vector *v,int x,int z){
	int i,dim2;
	GMERR(-1111);
	if((self==NULL)||(v==NULL))						GMERR(-1);
	if(v->size!=self->d2)							GMERR(-11);
	dim2=self->d2;
	for(i=0;i<dim2;i++){
		gsl_vector_set(v,i,gsl_matrix_get(self->matrixList[x],i,z));
	}

	return(0);
GMERRH("marray3d_get_pencilY",1);
}

int marray3d_get_pencilZ(marray3d *self,gsl_vector *v,int x,int y){
	int i,dim3;
	GMERR(-1111);
	if((self==NULL)||(v==NULL))						GMERR(-1);
	if(v->size!=self->d3)							GMERR(-11);
	dim3=self->d3;
	for(i=0;i<dim3;i++){
		gsl_vector_set(v,i,gsl_matrix_get(self->matrixList[x],y,i));
	}
	
	return(0);
GMERRH("marray3d_get_pencilZ",1);
}

int marray3d_set_pencilX(marray3d *self,gsl_vector *v,int y,int z){
	int i,dim1;
	GMERR(-1111);
	if((self==NULL)||(v==NULL))						GMERR(-1);
	if(v->size!=self->d1)							GMERR(-11);
	dim1=self->d1;
	for(i=0;i<dim1;i++){
		gsl_matrix_set(self->matrixList[i],y,z,gsl_vector_get(v,i));
	}

	return(0);
GMERRH("marray3d_set_pencilX",1);
}

int marray3d_set_pencilY(marray3d *self,gsl_vector *v,int x,int z){
	int i,dim2;
	GMERR(-1111);
	if((self==NULL)||(v==NULL))						GMERR(-1);
	if(v->size!=self->d2)							GMERR(-11);
	dim2=self->d2;
	for(i=0;i<dim2;i++){
		gsl_matrix_set(self->matrixList[x],i,z,gsl_vector_get(v,i));
	}

	return(0);
GMERRH("marray3d_set_pencilY",1);
}

int marray3d_set_pencilZ(marray3d *self,gsl_vector *v,int x,int y){
	int dim3,i;
	if((self==NULL)||(v==NULL))						GMERR(-1);
	if(v->size!=self->d3)							GMERR(-11);
	dim3=self->d3;
	gsl_matrix_set_row(self->matrixList[x],y,v);

	return(0);
GMERRH("marray3d_set_pencilZ",1);
}

int marray3d_get_matrixSliceX(marray3d *self,gsl_matrix *m,int x,int y1,int y2,gsl_vector *work){
	int t,n;
	if((self==NULL)||(m==NULL))											GMERR(-1);
	if((x<0)||(x>=self->d1))											GMERR(-11);
	if((y1<0)||(y2>self->d3))											GMERR(-21);
	if((m->size1!=self->d2)||(m->size2!=(y2-y1)))						GMERR(-31);
	if((work->size!=self->d2))											GMERR(-41);
	n = y2-y1;
	for(t=0;t<n;t++){
		if(gsl_matrix_get_col(work,self->matrixList[x],t+y1))			GMERR(-51);
		if(gsl_matrix_set_col(m,t,work))								GMERR(-61);
	}
	return(0);
GMERRH("marray3d_get_matrixSliceX",1);
}

int marray3d_set_matrixSliceX(marray3d *self,gsl_matrix *m,int x,int y1,int y2,gsl_vector *work){
	int t,n;
	if((self==NULL)||(m==NULL))						GMERR(-1);
	if((x<0)||(x>=self->d1))						GMERR(-11);
	if((y1<0)||(y2>self->d3))						GMERR(-21);
	if((m->size1!=self->d2)||(m->size2!=(y2-y1)))	GMERR(-31);
	if((work->size!=self->d2))											GMERR(-41);
	n = y2-y1;
	for(t=0;t<n;t++){
		if(gsl_matrix_get_col(work,m,t))								GMERR(-61);
		if(gsl_matrix_set_col(self->matrixList[x],t+y1,work))			GMERR(-51);
	}
	return(0);
GMERRH("marray3d_set_matrixSliceX",1);
}
/*****************
* Linear algebra *
*****************/

int marray3d_LU(marray3d *m,marray3d *result,gsl_permutation *p,int *s){
	if((m==NULL)||(result==NULL)||(p==NULL))	GMERR(-1);
	if(p->size!=m->d2)							GMERR(-11);
	if(m->d3!=m->d2)							GMERR(-11);
	return(0);
GMERRH("marray3d_LU",1);
}

int marray3d_cholesky(marray3d *m,marray3d *result){
	return(0);
GMERRH("marray3d_cholesky",1);
}

int marray3d_inv(marray3d *m,marray3d *result){
	return(0);
GMERRH("marray3d_inv",1);
}

int marray3d_eigenValues(marray3d *m,gsl_matrix *evalues){
	return(0);
GMERRH("marray3d_eigenValues",1);
}

int marray3d_eigenVectors(marray3d *m,gsl_matrix *evalues,marray3d *evecs){
	return(0);
GMERRH("marray3d_eigenVectors",1);
}

/**************************
* Mathematical Operations *
**************************/

//works
int marray3d_add(marray3d *m1,marray3d *m2){
	int i,d1;
	if((m1==NULL)||(m2==NULL))						GMERR(-1);
	if(marray3d_compareSizes(m1,m2))				GMERR(-11);
	d1 = m1->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_add(m1->matrixList[i],m2->matrixList[i]);
	}
	return(0);
GMERRH("marray3d_add",1);
}

//works
int marray3d_sub(marray3d *m1,marray3d *m2){
	int i,d1;
	if((m1==NULL)||(m2==NULL))						GMERR(-1);
	if(marray3d_compareSizes(m1,m2))				GMERR(-11);
	d1 = m1->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_sub(m1->matrixList[i],m2->matrixList[i]);
	}
	
	return(0);
GMERRH("marray3d_sub",1);
}

//works
int marray3d_scale(marray3d *m1,double scalar){
	int i,d1;
	if(m1==NULL)									GMERR(-1);
	d1 = m1->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_scale(m1->matrixList[i],scalar);
	}
	return(0);
GMERRH("marray3d_scale",1);
}

//works
int marray3d_cpy(marray3d *m1,marray3d *m2){
	int i,d1;
	if((m1==NULL)||(m2==NULL))						GMERR(-1);
	if(marray3d_compareSizes(m1,m2))				GMERR(-11);
	d1 = m1->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_memcpy(m1->matrixList[i],m2->matrixList[i]);
	}

	return(0);
GMERRH("marray3d_cpy",1);
}

//works
int marray3d_pointwiseMult(marray3d *m1,marray3d *m2){
	int i,d1;
	if((m1==NULL)||(m2==NULL))						GMERR(-1);
	if(marray3d_compareSizes(m1,m2))				GMERR(-11);
	d1 = m1->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_mul_elements(m1->matrixList[i],m2->matrixList[i]);
	}
	return(0);
GMERRH("marray3d_pointwiseMult",1);
}

//works
int marray3d_pointwiseDivide(marray3d *m1,marray3d *m2){
	int i,d1;
	if((m1==NULL)||(m2==NULL))						GMERR(-1);
	if(marray3d_compareSizes(m1,m2))				GMERR(-11);
	d1 = m1->d1;
	for(i=0;i<d1;i++){
		gsl_matrix_div_elements(m1->matrixList[i],m2->matrixList[i]);
	}
	return(0);
GMERRH("marray3d_pointwiseDivide",1);
}

int marray3d_mulX(marray3d *m1,marray3d *m2,marray3d *m3){
	int d1,i;
	if((m1==NULL)||(m2==NULL)||(m3==NULL))			GMERR(-1);
	d1 = m1->d1;
	if((m2->d1!=d1)||(m3->d1!=d1))					GMERR(-11);
	if(m1->d3!=m2->d2)								GMERR(-21);
	if(m1->d2!=m3->d2)								GMERR(-31);
	if(m2->d3!=m3->d3)								GMERR(-41);
	for(i=0;i<d1;i++){
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1->matrixList[i],m2->matrixList[i],0.0,m3->matrixList[i]);
	}
	
	return(0);
GMERRH("marray3d_mulX",1);
}

int marray3d_mulY(marray3d *m1,marray3d *m2,marray3d *m3){
	if((m1==NULL)||(m2==NULL)||(m3==NULL))			GMERR(-1);
	printf("Should use marray3d_mulX\n");
	GMERR(-11);
	return(0);
GMERRH("marray3d_mulY",1);
}

int marray3d_mulZ(marray3d *m1,marray3d *m2,marray3d *m3){
	if((m1==NULL)||(m2==NULL)||(m3==NULL))			GMERR(-1);
	printf("Should use marray3d_mulX\n");
	GMERR(-11);
	return(0);
GMERRH("marray3d_mulZ",1);
}

int marray3d_mulMatY(marray3d *m1,gsl_matrix *m2,marray3d *m3){
	if((m1==NULL)||(m2==NULL)||(m3==NULL))			GMERR(-1);
	printf("Should use marray3d_mulMatX\n");
	GMERR(-11);
	return(0);
GMERRH("marray3d_mulMatY",1);
}

int marray3d_mulMatX(marray3d *m1,gsl_matrix *m2,marray3d *m3){
	int i,d1;
	if((m1==NULL)||(m2==NULL)||(m3==NULL))			GMERR(-1);
	if(m1->d1!=m3->d1)								GMERR(-11);
	if(m1->d3!=m2->size1)							GMERR(-21);
	if(m1->d2!=m3->d2)								GMERR(-31);
	if(m2->size2!=m3->d3)							GMERR(-41);
	d1 = m1->d1;
	for(i=0;i<d1;i++){
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1->matrixList[i],m2,0.0,m3->matrixList[i]);
	}

	return(0);
GMERRH("marray3d_mulMatX",1);
}

int marray3d_mulMatZ(marray3d *m1,gsl_matrix *m2,marray3d *m3){
	if((m1==NULL)||(m2==NULL)||(m3==NULL))			GMERR(-1);
	printf("Should use marray3d_mulMatX\n");
	return(0);
GMERRH("marray3d_mulMatZ",1);
}

/************************************
*************************************
**								   **
** Methods for 4 dimensional array **
**								   **
*************************************
************************************/

/************************
* Dimension comparisons *
************************/
int marray4d_checkSizes(marray4d *self,size_t d1,size_t d2,size_t d3,size_t d4){
	return((self==NULL)||(self->d1!=d1)||(self->d2!=d2)||(self->d3!=d3)||(self->d4!=d4));
}

int marray4d_checkYZW(marray4d *self,size_t d2,size_t d3,size_t d4){
	return((self==NULL)||(self->d2!=d2)||(self->d3!=d3)||(self->d4!=d4));
}

int marray4d_checkX1orT(marray4d *self,size_t t){
	return((self==NULL)||((self->d1!=t)&&(self->d1!=1)));
}

int marray4d_checkY1orT(marray4d *self,size_t t){
	return((self==NULL)||((self->d2!=t)&&(self->d2!=1)));
}

int marray4d_checkZ1orT(marray4d *self,size_t t){
	return((self==NULL)||((self->d3!=t)&&(self->d3!=1)));
}

int marray4d_checkW1orT(marray4d *self,size_t t){
	return((self==NULL)||((self->d4!=t)&&(self->d4!=1)));
}

int marray4d_compareSizes(marray4d *m1,marray4d *m2){
	return((m1==NULL)||(m2==NULL)||(m1->d1!=m2->d1)||(m1->d2!=m2->d2)||(m1->d3!=m2->d3)||(m1->d4!=m2->d4));
}

/******************************
* Allocation and free methods *
******************************/
marray4d * marray4d_alloc(size_t d1,size_t d2,size_t d3,size_t d4){
	int i;
	marray4d *self = (marray4d *)(GM_Malloc(sizeof(marray4d)));
	if(self==NULL)				GMERR(-1);
	self->d1 = d1;
	self->d2 = d2;
	self->d3 = d3;
	self->d4 = d4;
	self->marrayList = (marray3d **)GM_Malloc(d1*sizeof(marray3d));
	for(i=0;i<d1;i++){
		self->marrayList[i] = marray3d_alloc(d2,d3,d4);
		if(self->marrayList[i]==NULL)	GMERR(-2);
	}
	return(self);
GMERRH("marray4d_alloc",NULL);
}

int marray4d_free(marray4d *self){
	int i;
	for(i=0;i<self->d1;i++){
		marray3d_free(self->marrayList[i]);
		self->marrayList[i] = NULL;
	}
	GM_Free(self->marrayList);
	self->marrayList = NULL;
	GM_Free(self);
	return(0);
GMERRH("marray4d_free",1);
}

/*************************
* Initialization methods *
*************************/

int marray4d_set_zero(marray4d *self){
	int i,d1;
	if(self==NULL)						GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		if(marray3d_set_zero(self->marrayList[i])) GMERR(-11);
	}
	return(0);
GMERRH("marray4d_set_zero",1);
}

int marray4d_set_constant(marray4d *self,double x){ 
	int i,d1;
	if(self==NULL)						GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		if(marray3d_set_constant(self->marrayList[i],x))GMERR(-11);
	}
	return(0);
GMERRH("marray4d_set_constant",1);
}

int marray4d_set_identity(marray4d *self){
	int i,d1;
	if(self==NULL)						GMERR(-1);
	if(self->d3!=self->d4)				GMERR(-11);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		if(marray3d_set_identity(self->marrayList[i]))	GMERR(-21);
	}
	return(0);
GMERRH("marray4d_set_identity",1);
}

int marray4d_set_marray3d(marray4d *self,marray3d *m){
	int i,d1;
	if((self==NULL)||(m==NULL))			GMERR(-1);
	if((self->d2!=m->d1)||(self->d3!=m->d2)||(self->d4!=m->d3))GMERR(-11);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		if(marray3d_cpy(self->marrayList[i],m))	GMERR(-21);
	}
	return(0);
GMERRH("marray4d_set_marray3d",1);
}

/*************************
* Read and write methods *
*************************/

int marray4d_read(FILE *fp,marray4d *self){
	int d1,d2,d3,d4,l,i,j,k,bufflen,maxbufflen=1024;
	int size1=0,size2=0,size3=0,size4=0;
	double val=0;
	char buff[1024];
	if((fp==NULL)||(self==NULL))					GMERR(-1);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	d4 = self->d4;
	if(fgets(buff,1023,fp)==NULL)		GMERR(-11);
	if(buff[0]!='#')					GMERR(-21);
	sscanf(buff+2,"%d %d %d %d",&size1,&size2,&size3,&size4);
	if((size1!=d1)||(size2!=d2)||(size3!=d3)||(size4!=d4))	GMERR(-31);
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				for(l=0;l<d4;l++){
					fscanf(fp,"%lf",&val);
					marray4d_set(self,i,j,k,l,val);
				}
			}
		}
	}
	return(0);
GMERRH("marray4d_read",1);
}

int marray4d_write(FILE *fp,marray4d *self){
	int d1,d2,d3,d4,i,j,k,l;
	if((fp==NULL)||(self==NULL))					GMERR(-1);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	d4 = self->d4;
	fprintf(fp,"# %d %d %d %d\n",d1,d2,d3,d4);
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				for(l=0;l<d4;l++){
					fprintf(fp,"%0.15f\n",marray4d_get(self,i,j,k,l));
				}
			}
		}
	}
	return(0);
GMERRH("marray4d_write",1);
}

int marray4d_print(marray4d *self){
	int d1,i;
	if(self==NULL)		GMERR(-1);
	d1 = self->d1;
	for(i=0;i<d1;i++){
		if(marray3d_print(self->marrayList[i]))	GMERR(-11);
	}
	return(0);
GMERRH("marray4d_print",1);
}

int marray4d_fancyWrite(FILE *fp,marray4d *self){
	int i,j,k,l;
	int d1,d2,d3,d4;
	if((fp==NULL)||(self==NULL))					GMERR(-1);
	d1 = self->d1;
	d2 = self->d2;
	d3 = self->d3;
	d4 = self->d4;
	for(i=0;i<d1;i++){
		for(j=0;j<d2;j++){
			for(k=0;k<d3;k++){
				for(l=0;l<d4;l++){
					fprintf(fp,"%0.15f ",marray4d_get(self,i,j,k,l));
				}
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n\n");
	}
	return(0);
GMERRH("marray4d_fancyWrite",1);
}

/******************************************************
* Random initialization with a common random variable *
******************************************************/

int marray4d_set_randN(marray4d *self,double mu,double sigma2){
	return(0);
GMERRH("marray4d_set_randN",1);
}

int marray4d_set_randG(marray4d *self,double rate,double scale){
	return(0);
GMERRH("marray4d_set_randG",1);
}

int marray4d_set_randIG(marray4d *self,double rate,double scale){
	return(0);
GMERRH("marray4d_set_randIG",1);
}

int marray4d_set_randU(marray4d *self,double a,double b){ 
	return(0);
GMERRH("marray4d_set_randU",1);
}


/**********************
* Setters and getters *
**********************/

//Individual points

double marray4d_get(marray4d *self,size_t x,size_t y,size_t z,size_t w){
	return(marray3d_get(self->marrayList[x],y,z,w));
}

int marray4d_set(marray4d *self,size_t x,size_t y,size_t z,size_t w,double val){
	if(self==NULL)									GMERR(-1);
	if(marray3d_set(self->marrayList[x],y,z,w,val))	GMERR(-11);
	return(0);
GMERRH("marray4d_set",1);
}

//Slices

//Set
int marray4d_set_ZW(marray4d *self,gsl_matrix *m,int x,int y){
	if((self==NULL)||(m==NULL))			GMERR(-1);
	if((self->d3!=m->size1)||(self->d4!=m->size2))	GMERR(-11);
	if(marray3d_set_X(self->marrayList[x],m,y))	GMERR(-21);
	return(0);
GMERRH("marray4d_set_YZ",1);
}

int marray4d_set_X(marray4d *self,marray3d *m,int x){
	if((self==NULL)||(m==NULL))			GMERR(-1);
	if((x<0||self->d1<x))				GMERR(-11);
	if(marray3d_cpy(self->marrayList[x],m))			GMERR(-21);
	return(0);
GMERRH("marray4d_set_X",1);
}

//get
int marray4d_get_ZW(marray4d *self,gsl_matrix *m,int x,int y){
	if((self==NULL)||(m==NULL))			GMERR(-1);
	if((self->d3!=m->size1)||(self->d4!=m->size2))	GMERR(-11);
	if(marray3d_get_X(self->marrayList[x],m,y))	GMERR(-21);
	return(0);
GMERRH("marray4d_get_YZ",1);
}

int marray4d_get_X(marray4d *self,marray3d *m,int x){
	if((self==NULL)||(m==NULL))			GMERR(-1);
	if((x<0||self->d1<x))				GMERR(-11);
	if(marray3d_cpy(m,self->marrayList[x]))			GMERR(-21);
	return(0);
GMERRH("marray4d_set_X",1);
}

/*****************
* Hybrid methods *
*****************/

int marray34_compareYZW(marray4d *m4,marray3d *m3){
	return((m4==NULL)||(m3==NULL)||(m4->d2!=m3->d1)||(m4->d3!=m3->d2)||(m4->d4!=m3->d3));	
}



