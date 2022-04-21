#include <math.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include "matrix_utils.h"
#include "gmlib.h"

/* These are my customized allocation methods. These are 
 * combined with the macro given in matrix_utils.h. This
 * lets me test my code using my own allocation methods
 * so memory leaks are easier to track down. Once those 
 * are fixed I can switch the macro being used in 
 * matrix_utils.h and use the original allocation method
 * which I assume is slightly faster. Thus when I compile
 * the final version I can use the fastest methods. */

gsl_vector * D_gsl_vector_alloc(int n){
	gsl_vector *result = gsl_vector_alloc(n);
	return(result);
 }

gsl_matrix * D_gsl_matrix_alloc(int n1,int n2){
	gsl_matrix *result = gsl_matrix_alloc(n1,n2);
	return(result);
}

int matrixCheckSize(gsl_matrix *m,int x,int y){
	return((m==NULL)||(m->size1!=x)||(m->size2!=y));
}

int nL(){
	printf("\n");
	return(0);
}

/*
* Slicing methods meant to recreate python
*/

int gsl_vector_slice(gsl_vector *v,int start,int end,gsl_vector *out){
	int n=end-start;
	int i;
	if(vectorCheckSize(out,n))					GMERR(-11);
	for(i=0;i<n;i++){
		gsl_vector_set(out,i,gsl_vector_get(v,i+start));
	}
	return(0);
GMERRH("gsl_vector_slice",1);
}

int gsl_matrix_slice(gsl_matrix *m,int x1,int x2,int y1,int y2,gsl_matrix *out){
	int n1,n2,i,j;
	if(x2==-1)	n1 = m->size1-x1;
	else		n1 = x2-x1;
	if(y2==-1)	n2 = m->size2-y2;
	else		n2 = y2-y1;
	if(matrixCheckSize(out,n1,n2))				GMERR(-11);
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			gsl_matrix_set(out,i,j,gsl_matrix_get(m,i+x1,j+y1));	
		}
	}

	return(0);
GMERRH("gsl_matrix_slice",1);
}

int reverseVector(gsl_vector *v,gsl_vector *rv){
	int t,n;
	double temp;
	if((v==NULL)||(rv==NULL)||(v->size!=rv->size))	GMERR(-11);
	n=v->size;

	for(t=0;t<n;t++){
		temp = gsl_vector_get(v,t);
		gsl_vector_set(rv,n-t-1,temp);
	}
	return(0);
GMERRH("gsl_vector_reverse",1);
}

int mpmTd2(gsl_matrix *m){
	int i,j,n;
	double temp;
	if(m==NULL||(m->size1!=m->size2))				GMERR(-11);
	n = m->size1;

	for(i=0;i<n;i++){
		for(j=i;j<n;j++){
			temp = gsl_matrix_get(m,i,j)+gsl_matrix_get(m,j,i);
			gsl_matrix_set(m,i,j,temp);
			gsl_matrix_set(m,j,i,temp);
		}
	}
	gsl_matrix_scale(m,0.5);
	return(0);
GMERRH("mpmTd2",1);
}

/*
int matrixCheckSize(gsl_matrix *m,int x,int y){
	int sX = m->size1;
	int sY = m->size2;
	if(sX!=x)									GMERR(-11);
	if(sY!=y)									GMERR(-21);
	return(0);
GMERRH("matrixCheckSize",1);
}*/

int vectorCheckSize(gsl_vector *v,int x){
	return((v==NULL)||(v->size!=x));
}


int cmatrixIsZero(const gsl_matrix *m){
	int n1=m->size1;
	int n2=m->size2;
	int i,j;
	double temp;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			temp = gsl_matrix_get(m,i,j);
			if(temp!=0.0)			goto EXIT;
		}
	}

	return(1);
EXIT:
	return(0);
}

int cvectorIsZero(const gsl_vector *v){
	int n = v->size;
	int i;
	double temp;
	for(i=0;i<n;i++){
		temp = gsl_vector_get(v,i);
		if(temp!=0.0)				goto EXIT;
	}
	return(1);
EXIT:
	return(0);
}
int matrixIsZero(gsl_matrix *m){
	int n1 = m->size1;
	int n2 = m->size2;
	int i,j;
	double temp;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			temp = gsl_matrix_get(m,i,j);
			if(temp!=0.0)			goto EXIT;
		}
	}

	return(1);
EXIT:
	return(0);
}
int vectorIsZero(gsl_vector *v){
	int n = v->size;
	int i;
	double temp;
	for(i=0;i<n;i++){
		temp = gsl_vector_get(v,i);
		if(temp!=0.0)				goto EXIT;
	}

	return(1);
EXIT:
	return(0);
}


/*This sets a diagonal matrix of the size of vector *v*/

//Works
int printGSLVector(gsl_vector *v){
	int i;
	int n1;
	if(v==NULL)									GMERR(-11);
	n1 = v->size;
	for(i=0;i<n1;i++){
		printf("%0.15f\n",gsl_vector_get(v,i));
	}
	return(0);
GMERRH("printGSLVector",1);
}

int printGSLVectorR(gsl_vector *v){
	int i;
	int n1;
	if(v==NULL)									GMERR(-11);
	n1 = v->size;
	for(i=0;i<n1;i++){
		printf("%0.2f\n",gsl_vector_get(v,i));
	}
	return(0);
GMERRH("printGSLVectorR",1);
}

//Works
int printGSLVectorT(gsl_vector *v){
	int i;
	int n1;
	if(v==NULL)									GMERR(-11);
	n1 = v->size;
	for(i=0;i<n1;i++){
		printf("%0.15f ",gsl_vector_get(v,i));
	}
	printf("\n");
	return(0);
GMERRH("printGSLVectorT",1);
}

int printGSLVectorRT(gsl_vector *v){
	int i;
	int n1;
	if(v==NULL)									GMERR(-11);
	n1 = v->size;
	for(i=0;i<n1;i++){
		printf("%0.2f ",gsl_vector_get(v,i));
	}
	printf("\n");
	return(0);
GMERRH("printGSLVectorRT",1);
}

//Works
int printGSLMatrix(gsl_matrix *m){
	int i,j;
	int n1,n2;
	if(m==NULL)									GMERR(-11);
	n1 = m->size1;
	n2 = m->size2;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			printf("%0.15f ",gsl_matrix_get(m,i,j));
		}
		printf("\n");
	}
	return(0);
GMERRH("printGSLMatrix",1);
}

int printGSLMatrixR(gsl_matrix *m){
	int i,j;
	int n1,n2;
	if(m==NULL)									GMERR(-11);
	n1 = m->size1;
	n2 = m->size2;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			printf("%0.2f ",gsl_matrix_get(m,i,j));
		}
		printf("\n");
	}
	return(0);
GMERRH("printGSLMatrixR",1);
}

int printGSLVectorSize(gsl_vector *v){
	int s;
	if(v==NULL)				GMERR(-1);
	s = v->size;
	printf("%d\n",s);
	return(0);
GMERRH("printGSLVectorSize",1);
}

int printGSLMatrixSize(gsl_matrix *m){
	int m1,m2;
	if(m==NULL)				GMERR(-1);
	m1 = m->size1;
	m2 = m->size2;
	printf("%d %d\n",m1,m2);
	return(0);
GMERRH("printGSLMatrixSize",1);
}

int printGSLVectorComplex(gsl_vector_complex *v){
	return(0);
GMERRH("printGSLVectorComplex",1);
}

int printGSLVectorComplexT(gsl_vector_complex *v){
	return(0);
GMERRH("printGSLVectorComplexT",1);
}

int printGSLMatrixComplex(gsl_matrix_complex *m){
	return(0);
GMERRH("printGSLMatrixComplex",1);
}


//Works
gsl_matrix * gsl_matrix_set_diagonal(gsl_vector *v){
	int N = v->size;
	int i;
	gsl_matrix *m = gsl_matrix_alloc(N,N);

	if(v==NULL)									GMERR(-11);
	gsl_matrix_set_zero(m);
	for(i=0;i<N;i++){
		gsl_matrix_set(m,i,i,gsl_vector_get(v,i));
	}
	return(m);
GMERRH("gsl_amtrix_set_diagonal",NULL);
}

int positiveVector(gsl_vector *v){
	int i;
	int N = v->size;
	for(i=0;i<N;i++){
		if(gsl_vector_get(v,i)<=0)	goto MYRETURN;
	}
	return(0);
MYRETURN:	return(1);
}

int positiveMatrix(gsl_matrix *m){
	int i,j;
	int N1=m->size1;int N2=m->size2;
	for(i=0;i<N1;i++){
		for(j=0;j<N2;j++){
			if(gsl_matrix_get(m,i,j)<=0)	goto MYRETURN;
		}
	}
	return(0);
MYRETURN:	return(1);
}

//Returns a one if not positive definite. I did
//this because it is mainly used to see if a 
//covariance matrix is valid while error checking

//Returns 1 if not symmetric
int symmetricMatrix(gsl_matrix *m){
	int n = m->size1;
	int i,j;
	if(m->size1!=m->size2)					goto MYRETURN;
	for(i=0;i<n;i++){
		for(j=i;j<n;j++){
			if(gsl_matrix_get(m,i,j)!=gsl_matrix_get(m,j,i))	goto MYRETURN;
		}
	}
	return(0);
MYRETURN:	return(1);
}

int positiveDefinite(gsl_matrix *m){
	int n = m->size1;
	gsl_vector *eval = gsl_vector_alloc(n);
	gsl_matrix *m_cpy = gsl_matrix_alloc(n,n);
	gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(n);
	if(m->size1!=m->size2)		GMERR(-11);
	gsl_matrix_memcpy(m_cpy,m);
	//ensures real eigenvalues
	gsl_eigen_symm(m_cpy,eval,w);
	if(gsl_vector_ispos(eval)==0)	GMERR(-31);
	//Free the temporary variables
	gsl_eigen_symm_free(w);
	gsl_matrix_free(m_cpy);
	gsl_vector_free(eval);
	return(0);
GMERRH("positiveDefinite",1);
}
//Works
gsl_matrix * gsl_matrix_set_blockDiagonal2(gsl_matrix *m1,gsl_matrix *m2){
	int n11 = m1->size1;int n21 = m2->size1;
	int n12 = m1->size2;int n22 = m2->size2;
	int N1 = n11+n21;
	int N2 = n12+n22;
	int i,j;
	gsl_matrix *m;
	if(m1==NULL)									GMERR(-11);
	if(m2==NULL)									GMERR(-21);
	m = gsl_matrix_alloc(N1,N2);
	gsl_matrix_set_zero(m);
	for(i=0;i<n11;i++){
		for(j=0;j<n12;j++){
			gsl_matrix_set(m,i,j,gsl_matrix_get(m1,i,j));
		}
	}
	for(i=0;i<n21;i++){
		for(j=0;j<n22;j++){
			gsl_matrix_set(m,i+n11,j+n12,gsl_matrix_get(m2,i,j));
		}
	}
	return(m);
GMERRH("gsl_matrix_set_blockDiagonal2",NULL);
}

//Should work
gsl_matrix * gsl_matrix_set_blockDiagonal3(gsl_matrix *m1,gsl_matrix *m2,gsl_matrix *m3){
	int n11 = m1->size1;int n21 = m2->size1;int n31 = m3->size1;
	int n12 = m1->size2;int n22 = m2->size2;int n32 = m3->size2;
	int N1 = n11+n21+n31;
	int N2 = n12+n22+n32;
	int i,j;
	gsl_matrix *m;
	if(m1==NULL)									GMERR(-11);
	if(m2==NULL)									GMERR(-21);
	if(m3==NULL)									GMERR(-31);
	m = gsl_matrix_alloc(N1,N2);
	gsl_matrix_set_zero(m);
	for(i=0;i<n11;i++){
		for(j=0;j<n12;j++){
			gsl_matrix_set(m,i,j,gsl_matrix_get(m1,i,j));
		}
	}
	for(i=0;i<n21;i++){
		for(j=0;j<n22;j++){
			gsl_matrix_set(m,i+n11,j+n12,gsl_matrix_get(m2,i,j));
		}
	}
	for(i=0;i<n31;i++){
		for(j=0;j<n32;j++){
			gsl_matrix_set(m,i+n11+n21,j+n12+n22,gsl_matrix_get(m2,i,j));
		}
	}
	return(m);
GMERRH("gsl_matrix_set_blockDiagonal3",NULL);
}

//Works
gsl_vector * gsl_vector_set_ones(int N){
	gsl_vector *v;
	int i;
	if(N<1)											GMERR(-11);
	v = gsl_vector_alloc(N);
	gsl_vector_set_zero(v);
	for(i=0;i<N;i++){
		gsl_vector_set(v,i,1.0);
	}
	return(v);
GMERRH("gsl_vector_set_ones",NULL);
}

//Works
gsl_matrix * gsl_matrix_set_ones(int N1,int N2){
	gsl_matrix *m;
	int i,j;
	if(N1<1||N2<1)									GMERR(-11);
	m = gsl_matrix_alloc(N1,N2);
	gsl_matrix_set_zero(m);
	for(i=0;i<N1;i++){
		for(j=0;j<N2;j++){
			gsl_matrix_set(m,i,j,1.0);
		}
	}
	return(m);
GMERRH("gsl_matrix_set_ones",NULL);
}

int randnV(gsl_vector *v,gsl_rng *rand,double mean,double std){
	int i,n;
	double temp;
	n = v->size;
	for(i=0;i<n;i++){
		temp = gsl_ran_gaussian(rand,std) + mean;
		gsl_vector_set(v,i,temp);
	}
	return(0);
GMERRH("randnV",1);
}

int randnM(gsl_matrix *m,gsl_rng *rand,double mean,double std){
	int i,n1,n2,j;
	double temp;
	n1 = m->size1;
	n2 = m->size2;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			temp = gsl_ran_gaussian(rand,std) + mean;
			gsl_matrix_set(m,i,j,temp);
		}
	}
	return(0);
GMERRH("randnM",1);
}

int randgV(gsl_vector *v,gsl_rng *rand,double rate,double scale){
	int i,n;
	double temp;
	n = v->size;
	for(i=0;i<n;i++){
		temp = gsl_ran_gamma(rand,rate,scale);
		gsl_vector_set(v,i,temp);
	}
	return(0);
GMERRH("randgV",1);
}

int randgM(gsl_matrix *m,gsl_rng *rand,double rate,double scale){
	int i,j,n1,n2;
	double temp;
	n1 = m->size1;
	n2 = m->size2;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			temp = gsl_ran_gamma(rand,rate,scale);
			gsl_matrix_set(m,i,j,temp);
		}
	}
	return(0);
GMERRH("randgM",1);
}

int randigV(gsl_vector *v,gsl_rng *rand,double rate,double scale){
	int i;
	double temp;
	int n = v->size;
	for(i=0;i<n;i++){
		temp = 1.0/gsl_ran_gamma(rand,rate,1.0/scale);
		gsl_vector_set(v,i,temp);
	}
	return(0);
GMERRH("randigV",1);
}

int randigM(gsl_matrix *m,gsl_rng *rand,double rate,double scale){
	int i,j;
	double temp;
	int n1 = m->size1;
	int n2 = m->size2;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			temp = 1.0/gsl_ran_gamma(rand,rate,1.0/scale);
			gsl_matrix_set(m,i,j,temp);
		}
	}
	return(0);
GMERRH("randigM",1);
}

int randDiagN(gsl_matrix *m,gsl_rng *rand,double mean,double std){
	int n1 = m->size1;int n2 = m->size2;
	int N,i;
	double temp;
	if(n1<n2){
		N = n1;
	}else{
		N = n2;
	}
	for(i=0;i<N;i++){
		temp = gsl_ran_gaussian(rand,std) + mean;
		gsl_matrix_set(m,i,i,temp);
	}
	return(0);
GMERRH("randDiagN",1);
}

int randDiagG(gsl_matrix *m,gsl_rng *rand,double rate,double scale){
	int n1 = m->size1;int n2 = m->size2;
	int N,i;
	double temp;
	if(n1<n2){
		N = n1;
	}else{
		N = n2;
	}
	for(i=0;i<N;i++){
		temp = gsl_ran_gamma(rand,rate,scale);
		gsl_matrix_set(m,i,i,temp);
	}
	return(0);
GMERRH("randDiagG",1);
}

int randDiagIG(gsl_matrix *m,gsl_rng *rand,double rate,double scale){
	int n1 = m->size1;int n2 = m->size2;
	int N,i;
	double temp;
	if(n1<n2){
		N = n1;
	}else{
		N = n2;
	}
	for(i=0;i<N;i++){
		temp = 1.0/gsl_ran_gamma(rand,rate,1.0/scale);
		gsl_matrix_set(m,i,i,temp);
	}
	return(0);
GMERRH("randDiagIG",1);
}

/*Here are some common operations I want with vectors that for
some reason aren't included*/

int kron_mm(gsl_matrix *m1,gsl_matrix *m2,gsl_matrix *result){
	double x,y;
	int i,j,k,l;
	int s11,s12,s21,s22;
	if((m1->size1*m2->size1)!=result->size1)		GMERR(-11);
	if((m1->size2*m2->size2)!=result->size2)		GMERR(-21);
	s11 = m1->size1;
	s12 = m1->size2;
	s21 = m2->size1;
	s22 = m2->size2;
	for(i=0;i<s11;i++){
		for(j=0;j<s12;j++){
			x = gsl_matrix_get(m1,i,j);
			for(k=0;k<s21;k++){
				for(l=0;l<s22;l++){
					y = gsl_matrix_get(m2,k,l);
					gsl_matrix_set(result,i*s11+k,j*s12+l,x*y);
				}
			}
		}
	}
	return(0);	
GMERRH("kron_mm",1);
}

int hankel(gsl_vector *c,gsl_vector *r,gsl_matrix *out){
	int n1,n2,i,j,k;
	double temp;
	if((c==NULL)||(out==NULL))	GMERR(-1);
	gsl_matrix_set_zero(out);
	if(r!=NULL){if(gsl_vector_get(c,c->size-1)!=gsl_vector_get(r,0))
			printf("Warning: (hankel) c[end]!=r[0]\n");}
	if(r==NULL){
		n1 = c->size;
		n2 = c->size;
		if((out->size1!=n1)||(out->size2!=n2))	GMERR(-111);
		for(i=0;i<n1;i++){
			for(j=0;j<n1-i;j++){
				temp = gsl_vector_get(c,i+j);
				gsl_matrix_set(out,i,j,temp);
			}
		}

	}else{
		n1 = c->size;
		n2 = r->size;
		if((out->size1!=n1)||(out->size2!=n2))	GMERR(-211);
		if(n1>n2){
			for(j=0;j<n2;j++){
				temp = gsl_vector_get(r,n2-1-j);
				for(i=0;i<j+1;i++){
					gsl_matrix_set(out,n1-i-1,n2-j-1+i,temp);
				}
			}
			for(i=0;i<n1;i++){
				temp = gsl_vector_get(c,i);
				for(j=0;j<n2;j++){
					if(i-j>=0){
						gsl_matrix_set(out,i-j,j,temp);
					}
				}
			}

		}else{
			for(j=0;j<n2;j++){
				temp = gsl_vector_get(r,n2-1-j);
				for(i=0;i<j+1;i++){
					if((n1-i-1>=0)&&(n2-j-1+i>=0)){
						gsl_matrix_set(out,n1-i-1,n2-j-1+i,temp);
					}
				}
			}

			for(i=0;i<n1;i++){
				temp = gsl_vector_get(c,i);
				for(j=0;j<n2;j++){
					if(i-j>=0){
						gsl_matrix_set(out,i-j,j,temp);
					}
				}
			}
		}
	}
	 
	return(0);
GMERRH("hankel",1);
}

int meanVector(gsl_vector *v,double *mean){
	double d = 0;
	int i;
	int N = v->size;
	for(i=0;i<N;i++){
		d = d + gsl_vector_get(v,i);
	}
	d = d/N;
	*mean = d;

	return(0);
GMERRH("meanVector",1);
}

int meanMatrix(gsl_matrix *m,int axis,gsl_vector *v){
	double d = 0;
	int i,j;
	int n1=m->size1;
	int n2=m->size2;
	if(axis){
		for(i=0;i<n1;i++){
			d = 0;
			for(j=0;j<n2;j++){
				d = d + gsl_matrix_get(m,i,j);
			}
			d = d/n2;
			gsl_vector_set(v,i,d);
		}
	}else{
		for(j=0;j<n2;j++){
			d = 0;
			for(i=0;i<n1;i++){
				d = d + gsl_matrix_get(m,i,j);
			}
			d = d/n1;
			gsl_vector_set(v,j,d);
		}
	}
	
	return(0);
GMERRH("meanMatrix",1);
}


int demeanVector(gsl_vector *v,gsl_vector *v1){
	int i,N;
	double avg,temp;
	if((v==NULL)||(v1==NULL))					GMERR(-1);
	if(v->size!=v1->size)						GMERR(-11);
	N = v->size;
	for(i=0;i<N;i++){
		avg = avg + gsl_vector_get(v,i);
	}
	avg = avg/N;
	for(i=0;i<N;i++){
		temp = gsl_vector_get(v,i);
		temp = temp-avg;
		gsl_vector_set(v1,i,temp);
	}
	return(0);
GMERRH("demeanVector",1);
}

int demeanMatrix(gsl_matrix *m,int axis,gsl_matrix *m1){
	int N1,N2;
	int i,j;
	double avg,temp;
	if((m==NULL)||(m1==NULL))					GMERR(-1);
	N1 = m->size1;
	N2 = m->size2;
	if(matrixCheckSize(m1,N1,N2))				GMERR(-11);
	if(axis==1){
		for(i=0;i<N1;i++){
			avg = 0;
			for(j=0;j<N2;j++){
				avg += gsl_matrix_get(m,i,j);
			}
			avg = avg/N2;
			for(j=0;j<N2;j++){
				temp = gsl_matrix_get(m,i,j);
				temp = temp-avg;
				gsl_matrix_set(m,i,j,temp);
			}
		}
	}else if(axis==2){
		for(j=0;j<N2;j++){
			avg = 0;
			for(i=0;i<N1;i++){
				avg += gsl_matrix_get(m,i,j);
			}
			avg = avg/N1;
			for(i=0;i<N2;i++){
				temp = gsl_matrix_get(m,i,j);
				temp = temp-avg;
				gsl_matrix_set(m,i,j,temp);
			}
		}
	}else{
		GMERR(-999);
	}
	
	
	return(0);
GMERRH("demeanMatrix",1);
}




int kron_mv(gsl_matrix *m1,gsl_vector *v,gsl_matrix *result){

	return(0);	
GMERRH("kron_mm",1);
}


int kron_vm(gsl_vector *v,gsl_matrix *m,gsl_matrix *result){
	int i,j,k;
	int v1,m1,m2;
	double vi;
	double mij;
	if(result->size2!=(v->size*m->size2))		GMERR(-11);
	if(result->size1!=m->size1)					GMERR(-21);
	v1 = v->size;
	m1 = m->size1;
	m2 = m->size2;
	for(i=0;i<v1;i++){
		vi = gsl_vector_get(v,i);
		for(j=0;j<m1;j++){
			for(k=0;k<m2;k++){
				mij = gsl_matrix_get(m,j,k);
				gsl_matrix_set(result,j,i*v1+k,vi*mij);
			}
		}
	}
	
	return(0);	
GMERRH("kron_vm",1);
}


int kron_vv(gsl_vector *v1,gsl_vector *v2,gsl_matrix *result){
	
	return(0);	
GMERRH("kron_mm",1);
}

int sum_col(gsl_matrix *m,gsl_vector *v){
	int t,i;
	double sum;
	int m1,m2;
	if(m->size2!=v->size)		GMERR(-11);
	m1 = m->size1;
	m2 = m->size2;
	for(i=0;i<m2;i++){
		sum = 0;
		for(t=0;t<m1;t++){
			sum += gsl_matrix_get(m,t,i);
		}
		gsl_vector_set(v,i,sum);
	}
	return(0);
GMERRH("sum_col",1);
}

int sum_absCol(gsl_matrix *m,gsl_vector *v){
	int t,i;
	double sum;
	int m1,m2;
	if(m->size2!=v->size)		GMERR(-11);
	m1 = m->size1;
	m2 = m->size2;
	for(i=0;i<m2;i++){
		sum = 0;
		for(t=0;t<m1;t++){
			sum += fabs(gsl_matrix_get(m,t,i));
		}
		gsl_vector_set(v,i,sum);
	}
	return(0);
GMERRH("sum_absCol",1);
}

int sum_abs2Col(gsl_matrix *m,gsl_vector *v){
	int t,i;
	double sum,s;
	int m1,m2;
	if(m->size2!=v->size)		GMERR(-11);
	m1 = m->size1;
	m2 = m->size2;
	for(i=0;i<m2;i++){
		sum = 0;
		for(t=0;t<m1;t++){
			s = gsl_matrix_get(m,t,i);
			sum += s*s;
		}
		gsl_vector_set(v,i,sum);
	}
	return(0);
GMERRH("sum_absCol",1);
}

int sum_row(gsl_matrix *m,gsl_vector *v){
	int t,i;
	double sum;
	int m1,m2;
	if(m->size1!=v->size)		GMERR(-11);
	m1 = m->size1;
	m2 = m->size2;
	for(i=0;i<m1;i++){
		sum = 0;
		for(t=0;t<m2;t++){
			sum += gsl_matrix_get(m,i,t);
		}
		gsl_vector_set(v,i,sum);
	}
	return(0);
GMERRH("sum_col",1);
}

int sum_absRow(gsl_matrix *m,gsl_vector *v){
	int t,i;
	double sum;
	int m1,m2;
	if(m->size1!=v->size)		GMERR(-11);
	m1 = m->size1;
	m2 = m->size2;
	for(i=0;i<m1;i++){
		sum = 0;
		for(t=0;t<m2;t++){
			sum += fabs(gsl_matrix_get(m,i,t));
		}
		gsl_vector_set(v,i,sum);
	}
	return(0);
GMERRH("sum_absRow",1);
}

int sum_abs2Row(gsl_matrix *m,gsl_vector *v){
	int t,i;
	double sum,s;
	int m1,m2;
	if(m->size1!=v->size)		GMERR(-11);
	m1 = m->size1;
	m2 = m->size2;
	for(i=0;i<m1;i++){
		sum = 0;
		for(t=0;t<m2;t++){
			s = gsl_matrix_get(m,i,t);
			sum += s*s;
		}
		gsl_vector_set(v,i,sum);
	}
	return(0);
GMERRH("sum_abs2Row",1);
}

int cumProd(gsl_vector *v,gsl_vector *result){
	double prod=1;
	int t,Nt;
	if(v->size!=result->size)			GMERR(-11);
	Nt = v->size;
	for(t=0;t<Nt;t++){
		prod = prod*gsl_vector_get(v,t);
		gsl_vector_set(result,t,prod);
	}
	return(0);
GMERRH("cumProd",1);
}

int cumSum(gsl_vector *v,gsl_vector *result){
	double sum = 0;
	int t,Nt;
	if(v->size!=result->size)			GMERR(-11);
	Nt = v->size;
	for(t=0;t<Nt;t++){
		sum = sum + gsl_vector_get(v,t);
		gsl_vector_set(result,t,sum);
	}
	return(0);
GMERRH("cumSum",1);
}

int reshape(gsl_vector *v,int s1,int s2,gsl_matrix *m){
	int i,j,count;
	if(matrixCheckSize(m,s1,s2))		GMERR(-11);
	count = 0;
	for(j=0;j<s2;j++){
		for(i=0;i<s1;i++){
			gsl_matrix_set(m,i,j,gsl_vector_get(v,count));
			count = count + 1;
		}
	}
	
	
	return(0);
GMERRH("reshape",1);
}


//Works
double vecProd(gsl_vector *v){
	double prod = 1;
	int i,N;
	N = v->size;
	for(i=0;i<N;i++){
		prod = prod*gsl_vector_get(v,i);
	}
	return(prod);
}

//Works
double vecSum(gsl_vector *v){
	double sum = 0;
	int i,N;
	N = v->size;
	for(i=0;i<N;i++){
		sum = sum + gsl_vector_get(v,i);
	}
	return(sum);
}

//Works
int getDiag(gsl_matrix *m,gsl_vector *v){
	int n1 = m->size1;
	int n2 = m->size2;
	int n = GSL_MIN(n1,n2);
	int i;
	v = gsl_vector_alloc(n);
	gsl_vector_set_zero(v);
	for(i=0;i<n;i++){
		gsl_vector_set(v,i,gsl_matrix_get(m,i,i));
	}
	return(0);
GMERRH("getDiag",1);
}

//Works
int getOffDiag(gsl_matrix *m,gsl_vector *v,int offset){
	/*For now we assume that this is a square matrix*/
	int n = m->size1;
	int i,l;
	if(abs(offset)>=n)									GMERR(-11);
	if(m->size1!=m->size2)								GMERR(-21);
	l = n-abs(offset);
	v = gsl_vector_alloc(l);
	for(i=0;i<l;i++){
		if(offset>0){
			gsl_vector_set(v,i,gsl_matrix_get(m,i,i+offset));
		}else{
			gsl_vector_set(v,i,gsl_matrix_get(m,i-offset,i));
		}
	}
	return(0);
GMERRH("getOffDiag",1);
}

//Works
int setDiag(gsl_matrix *m,gsl_vector *v){
	/*For now we assume that this is a square matrix*/
	int i;
	int n = v->size;
	if((m->size1!=n)||(m->size2!=n))					GMERR(-11);
	for(i=0;i<n;i++){
		gsl_matrix_set(m,i,i,gsl_vector_get(v,i));
	}
	return(0);
GMERRH("setDiag",1);
}

//Works
int setOffDiag(gsl_matrix *m,gsl_vector *v,int offset){
	int i;
	int n =v->size;
	if((m->size1!=(n+abs(offset)))||(m->size2!=(n+abs(offset))))	GMERR(-11);
	for(i=0;i<n;i++){
		if(offset>0){
			gsl_matrix_set(m,i,i+offset,gsl_vector_get(v,i));
		}else{
			gsl_matrix_set(m,i-offset,i,gsl_vector_get(v,i));
		}
	}
	return(0);
GMERRH("setOffDiag",1);
}

int setDiagConst(gsl_matrix *m,double d){
	int i,n1,n2,minimum;
	if(m==NULL)							GMERR(-1);
	n1 = m->size1;
	n2 = m->size2;
	if(n1>n2){
		minimum = n2;
	}else{
		minimum = n1;
	}

	for(i=0;i<minimum;i++){
		gsl_matrix_set(m,i,i,d);
	}
	
	return(0);
GMERRH("setDiagConst",1);
}

int setOffDiagConst(gsl_matrix *m,double d,int offset){
	int i,n1,n2,minimum;
	if(m==NULL)							GMERR(-1);
	n1 = m->size1;
	n2 = m->size2;
	if(n1>n2){
		minimum = n2;
	}else{
		minimum = n1;
	}

	if(abs(offset)>minimum)				GMERR(-11);
	for(i=0;i<minimum;i++){
		if(offset>0){
			gsl_matrix_set(m,i,i+offset,d);
		}else{
			gsl_matrix_set(m,i-offset,i,d);
		}
	}
	return(0);
GMERRH("setOffDiagConst",1);
}

int sqrtVector(gsl_vector *v){
	int n1,i;
	if(v==NULL)							GMERR(-1);
	n1 = v->size;
	for(i=0;i<n1;i++){
		gsl_vector_set(v,i,sqrt(gsl_vector_get(v,i)));
	}
	return(0);
GMERRH("sqrtVector",1);
}

int sqrtMatrix(gsl_matrix *m){
	int n1,n2,i,j;
	if(m==NULL)							GMERR(-1);
	n1 = m->size1;
	n2 = m->size2;
	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			gsl_matrix_set(m,i,j,sqrt(gsl_matrix_get(m,i,j)));
		}
	}
	return(0);
GMERRH("sqrtMatrix",1);
}

//Works
double trace(gsl_matrix *m){
	double trace = 0;
	int n = GSL_MIN(m->size1,m->size2);
	int i;
	for(i=0;i<n;i++){
		trace = trace + gsl_matrix_get(m,i,i);
	}
	return(trace);
}


int diff(gsl_vector *v,gsl_vector *d){
	int i,n;
	if((v->size-1)!=d->size)			GMERR(-1);
	n = v->size-1;
	for(i=0;i<n;i++){
		gsl_vector_set(d,i,gsl_vector_get(v,i+1)-gsl_vector_get(v,i));
	}
	
	return(0);
GMERRH("diff",1);
}

int vector_diff2(gsl_vector *v1,gsl_vector *v2,double *result){
	int i,n;
	double temp,sum;
	if((v1==NULL)||(v2==NULL)||(v1->size!=v2->size))	GMERR(-1);
	n = v1->size;
	sum = 0;
	for(i=0;i<n;i++){
		temp = gsl_vector_get(v1,i)-gsl_vector_get(v2,i);
		temp = temp*temp;
		sum += temp;
	}
	*result = sum;
	return(0);
GMERRH("vector_diff2",1);
}


/* Here we want do matrix vector and vector vector multiplication
using standard notation of python. All that is required is specify
types of data when calling the function. This doesnt assume symmetry*/

//Works
int dot_vv(gsl_vector *v1,gsl_vector *v2,double* result){
	if(v1->size!=v2->size)							GMERR(-11);
	gsl_blas_ddot(v1,v2,result);
	return(0);
GMERRH("dot_vv",1);
}

//Works
int dot_mv(gsl_matrix *m,gsl_vector *v,gsl_vector *result){
	if(m->size2!=v->size)							GMERR(-11);
	if(m->size1!=result->size)						GMERR(-21);

	if(gsl_blas_dgemv(CblasNoTrans,1.0,m,v,0.0,result))	GMERR(-31);
	return(0);
GMERRH("dot_mv",1);
}

//Works
int dot_vm(gsl_vector *v,gsl_matrix *m,gsl_vector *result){
	gsl_matrix *mT = gsl_matrix_alloc(m->size2,m->size1);
	if(v->size!=m->size1)							GMERR(-11);
	if(m->size2!=result->size)						GMERR(-21);
	if(gsl_matrix_transpose_memcpy(mT,m))			GMERR(-31);

	if(gsl_blas_dgemv(CblasNoTrans,1.0,mT,v,0.0,result))	GMERR(-41);
	gsl_matrix_free(mT);
	return(0);
GMERRH("dot_vm",1);
}

//works
int dot_mm(gsl_matrix *m1,gsl_matrix *m2,gsl_matrix *result){
	if(m1->size2!=m2->size1)						GMERR(-11);
	if(m1->size1!=result->size1)					GMERR(-21);
	if(m2->size2!=result->size2)					GMERR(-31);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,m2,0.0,result))	GMERR(-41);
	return(0);
GMERRH("dot_mm",1);
}

//Works
int dot_amaT(gsl_matrix *a,gsl_matrix *m,gsl_matrix *result){
	gsl_matrix *aT = gsl_matrix_alloc(a->size2,a->size1);
	gsl_matrix *r1 = gsl_matrix_alloc(a->size1,m->size2);
	if(a->size2!=m->size1)							GMERR(-11);
	if(m->size1!=m->size2)							GMERR(-21);
	if(result->size1!=result->size2)				GMERR(-31);
	if(result->size1!=a->size1)						GMERR(-41);
	if(gsl_matrix_transpose_memcpy(aT,a))			GMERR(-46);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,a,m,0.0,r1))		GMERR(-51);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,r1,aT,0.0,result))	GMERR(-61);
	gsl_matrix_free(aT);							
	gsl_matrix_free(r1);
	return(0);
GMERRH("dot_amaT",1);
}

//Works
int dot_aTma(gsl_matrix *a,gsl_matrix *m,gsl_matrix *result){
	gsl_matrix *aT = gsl_matrix_alloc(a->size2,a->size1);
	gsl_matrix *r1 = gsl_matrix_alloc(a->size2,m->size1);
	if(a->size1!=m->size1)							GMERR(-11);
	if(m->size1!=m->size2)							GMERR(-21);
	if(result->size1!=result->size2)				GMERR(-31);
	if(a->size2!=result->size1)						GMERR(-41);
	if(gsl_matrix_transpose_memcpy(aT,a))			GMERR(-51);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,aT,m,0.0,r1))		GMERR(-61);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,r1,a,0.0,result))	GMERR(-71);
	gsl_matrix_free(aT);
	gsl_matrix_free(r1);
	return(0);
GMERRH("dot_aTma",1);
}

//Works
int dot_vTmv(gsl_vector *v,gsl_matrix *m,double *result){
	gsl_vector *r1 = gsl_vector_alloc(v->size);
	printGSLVectorT(v);
	nL();
	printGSLMatrix(m);
	if(v->size!=m->size1)							GMERR(-11);
	if(m->size1!=m->size2)							GMERR(-21);
	if(gsl_blas_dgemv(CblasNoTrans,1.0,m,v,0.0,r1)) GMERR(-31);
	if(gsl_blas_ddot(r1,v,result))					GMERR(-41);

	gsl_vector_free(r1);
	return(0);
GMERRH("dot_vTmv",1);
}

//Works
int dot_mTm(gsl_matrix *m,gsl_matrix *result){
	gsl_matrix *mT = gsl_matrix_alloc(m->size2,m->size1);
	if((m->size2!=result->size1)||(m->size2!=result->size2)) GMERR(-11);
	if(gsl_matrix_transpose_memcpy(mT,m))			GMERR(-21);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,mT,m,0.0,result))	GMERR(-31);
	gsl_matrix_free(mT);
	return(0);
GMERRH("dot_mTm",1);
}

//Works
int dot_mmT(gsl_matrix *m,gsl_matrix *result){
	gsl_matrix *mT = gsl_matrix_alloc(m->size2,m->size1);
	if(m->size1!=result->size1)						GMERR(-11);
	if(result->size1!=result->size2)				GMERR(-21);
	if(gsl_matrix_transpose_memcpy(mT,m))			GMERR(-31);
	if(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m,mT,0.0,result))	GMERR(-41);
	gsl_matrix_free(mT);
	return(0);
GMERRH("dot_mmT",1);
}

//Works
int dot_vvT(gsl_vector *v,gsl_matrix *result){
	int N = v->size;
	gsl_matrix *m = gsl_matrix_alloc(N,1);
	if(gsl_matrix_set_col(m,0,v))					GMERR(-11);
	if(dot_mmT(m,result))							GMERR(-21);
	gsl_matrix_free(m);
	return(0);
GMERRH("dot_vvT",1);
}

/*here are the norms*/

//Works
int vnorm2(gsl_vector *v,double *result){
	*result = gsl_blas_dnrm2(v);
	return(0);
GMERRH("vnorm2",1);
}

//Works
int vnorm1(gsl_vector *v,double *result){
	int i,N;
	double norm = 0;
	N = v->size;
	for(i=0;i<N;i++){
		norm = norm + fabs(gsl_vector_get(v,i));
	}
	*result = norm;
	return(0);
GMERRH("vnorm1",1);
}

//Works
int vnorm00(gsl_vector *v,double *result){
	int i,N;
	double j = 0;
	double norm = 0;
	N = v->size;
	for(i=0;i<N;i++){
		j = fabs(gsl_vector_get(v,i));
		if(j>norm)	norm=j;
	}
	*result = norm;
	return(0);
GMERRH("vnorm00",1);
}

//Does not work
int vnormp(gsl_vector *v,double p,double *result){
	int i,N;
	if(p<1)											GMERR(-11);
	double norm=0;
	N = v->size;
	for(i=0;i<N;i++){
		norm = norm + pow((gsl_vector_get(v,i)),p);
	}
	return(0);
GMERRH("vnormp",1);
}

int mnorm1(gsl_matrix *m,double *result){
	gsl_vector *col = gsl_vector_alloc(m->size1);
	double norm = 0;
	double colnorm = 0;
	int i;
	int N = m->size2;
	for(i=0;i<N;i++){
		gsl_matrix_get_col(col,m,i);
		if(vnorm1(col,&colnorm))					GMERR(-11);
		if(colnorm>norm)	norm=colnorm;
	}
	*result = norm;
	gsl_vector_free(col);
	return(0);
GMERRH("mnorm1",1);
}

int mnormF(gsl_matrix *m,double *result){
	gsl_matrix *mT =gsl_matrix_alloc(m->size2,m->size1);
	if(gTranspose(m,mT))							GMERR(-10);
	int N = m->size1;
	double tr = 0;
	gsl_matrix *mmt = gsl_matrix_alloc(N,N);
	if(dot_mm(mT,m,mmt))							GMERR(-11);
	tr = trace(mmt);
	tr = sqrt(tr);
	*result = tr;
	return(0);
GMERRH("mnormF",1);
}

int mnorm00(gsl_matrix *m,double *result){
	gsl_vector *row = gsl_vector_alloc(m->size2);
	double norm = 0;
	double rownorm = 0;
	int i;
	int N = m->size1;
	for(i=0;i<N;i++){
		gsl_matrix_get_row(row,m,i);
		if(vnorm1(row,&rownorm))					GMERR(-11);
		if(rownorm>norm)	norm=rownorm;
	}
	*result = norm;
	gsl_vector_free(row);
	return(0);
GMERRH("mnorm00",1);
}


/*Get matrix inverses. Use LU for nonsymmetric matrices
 *and chol for symmetric ones */
/*
int invLU(gsl_matrix *m,gsl_matrix *result){
	//if(gsl_linalg_LU_det())
	return(0);
GMERRH("invLU",1);
}

int invChol(gsl_matrix *m,gsl_matrix *result){

	return(0);
GMERRH("invChol",1);
}
*/
int det_LU(gsl_matrix *A,double *result){
	gsl_permutation *p = gsl_permutation_alloc(A->size1);
	gsl_matrix *tmpA = gsl_matrix_alloc(A->size1,A->size2);
	int signum;
	if(A->size1!=A->size2)							GMERR(-11);
	if(gsl_matrix_memcpy(tmpA,A))						GMERR(-21);
	if(gsl_linalg_LU_decomp(tmpA,p,&signum))		GMERR(-31);
	*result = gsl_linalg_LU_det(tmpA,signum);
	gsl_permutation_free(p);
	gsl_matrix_free(tmpA);
	return(0);
GMERRH("determinant",1);
}


int conditionNumber(gsl_matrix *m,double *result){

	return(0);
GMERRH("conditionNumber",1);
}

int get_symmetric_eig(gsl_matrix *m,gsl_vector *v){
	if(m->size1!=m->size2)							GMERR(-11);
	return(0);
GMERRH("get_symmetric_eig",1);
}


int get_asymmetric_eig(gsl_matrix *m,gsl_vector_complex *v){
	if(m->size1!=m->size2)							GMERR(-11);
	return(0);
GMERRH("get_asymmetric_eig",1);
}


int get_max_eig_rad(gsl_matrix *m,double *result){
	
	if(m->size1!=m->size2)							GMERR(-11);

	return(0);
GMERRH("get_max_eig_rad",1);
}

int get_min_eig_rad(gsl_matrix *m,double *result){
	if(m->size1!=m->size2)							GMERR(-11);

	return(0);
GMERRH("get_min_eig_rad",1);
}

int inv_LU(gsl_matrix *m,gsl_matrix *minv){
	int N = m->size1;
	gsl_permutation *p = gsl_permutation_alloc(N);
	gsl_matrix *mcpy = gsl_matrix_alloc(m->size1,m->size2);
	int s;
	if(m->size1!=m->size2)							GMERR(-11);
	if(minv->size1!=minv->size2)					GMERR(-21);
	if(minv->size1!=m->size1)						GMERR(-31);
	if(gsl_matrix_memcpy(mcpy,m))					GMERR(-41);
	if(gsl_linalg_LU_decomp(mcpy,p,&s))				GMERR(-51);
	if(gsl_linalg_LU_invert(mcpy,p,minv))			GMERR(-61);
	gsl_matrix_free(mcpy);
	gsl_permutation_free(p);
	return(0);
GMERRH("inv_LU",1);
}

int solve_LU(gsl_matrix *m,gsl_vector *b,gsl_vector *result){
	int N = m->size1;
	gsl_permutation *p = gsl_permutation_alloc(N);
	gsl_matrix *mcpy;
	int s;
	if(m->size1!=m->size2)							GMERR(-11);
	if(b->size!=m->size1)							GMERR(-21);
	if(result->size!=b->size)						GMERR(-31);
	if(gsl_matrix_memcpy(mcpy,m))					GMERR(-41);
	if(gsl_linalg_LU_decomp(mcpy,p,&s))				GMERR(-51);
	if(gsl_linalg_LU_solve(mcpy,p,b,result))		GMERR(-61);
	gsl_matrix_free(mcpy);
	gsl_permutation_free(p);
	return(0);
GMERRH("solveLU",1);
}

int inv_chol(gsl_matrix *m,gsl_matrix *minv){
	
	return(0);
GMERRH("inv_chol",1);
}

int vector_outer(gsl_vector *v,gsl_matrix *m){
	int i,j,n;
	double temp1,temp2;
	if((v==NULL)||(m==NULL))						GMERR(-1);
	n = v->size;
	if((m->size1!=n)||(m->size2!=n))				GMERR(-11);
	for(i=0;i<n;i++){
		temp1 = gsl_vector_get(v,i);
		for(j=0;j<n;j++){
			temp2 = gsl_vector_get(v,j);
			gsl_matrix_set(m,i,j,temp1*temp2);
		}
	}
	
	return(0);
GMERRH("vector_outer",1);
}


//Works
int gTranspose(gsl_matrix *m,gsl_matrix *mT){
	if(m->size1!=mT->size2)							GMERR(-11);
	if(m->size2!=mT->size1)							GMERR(-21);
	gsl_matrix_transpose_memcpy(mT,m);	
	return(0);
GMERRH("gTranspose",1);
}

int gHermetian(gsl_matrix_complex *c,gsl_matrix_complex *cH){
	if(c->size1!=cH->size2)							GMERR(-11);
	if(c->size2!=cH->size1)							GMERR(-21);
	return(0);
GMERRH("gHermetian",1);
}

//////////////////////////////////////////////////////////
int tgetSEigenvalues(gsl_matrix *m,gsl_vector *v){
	gsl_matrix *mcpy;
	gsl_eigen_symm_workspace *w;
	int n;
	if((m==NULL)||(v==NULL))						GMERR(-1);
	mcpy = gsl_matrix_alloc(m->size1,m->size2);
	n = m->size1;
	if(v->size!=n)									GMERR(-11);
	if(m->size1!=m->size2)							GMERR(-21);
	w = gsl_eigen_symm_alloc(n);
	gsl_matrix_memcpy(mcpy,m);
	if(gsl_eigen_symm(mcpy,v,w))					GMERR(-31);
	gsl_matrix_free(mcpy);
	gsl_eigen_symm_free(w);

	return(0);
GMERRH("tgetSEigenvalues",1);
}

int tgetSEigenvectors(gsl_matrix *m,gsl_vector *eigenvalues,gsl_matrix *eigenvectors){
	gsl_matrix *mcpy;
	gsl_eigen_symmv_workspace *w;
	int n;
	if((m==NULL)||(eigenvalues==NULL)||(eigenvectors==NULL))GMERR(-1);
	mcpy = gsl_matrix_alloc(m->size1,m->size2);
	n = m->size1;
	if((eigenvectors->size1!=m->size1))				GMERR(-21);
	if((eigenvectors->size2!=m->size2))				GMERR(-31);
	if(eigenvalues->size!=n)						GMERR(-41);
	
	gsl_eigen_symmv(mcpy,eigenvalues,eigenvectors,w);
	gsl_eigen_symmv_sort(eigenvalues,eigenvectors,GSL_EIGEN_SORT_ABS_ASC);

	gsl_eigen_symmv_free(w);
	gsl_matrix_free(mcpy);
	
	return(0);
GMERRH("tgetSEigenvectors",1);
}

int getSEigenvalues(gsl_matrix *m,gsl_matrix *mcpy,gsl_eigen_symm_workspace *w,gsl_vector *eval){
	int n;
	if((m==NULL)||(mcpy==NULL)||(w==NULL)||(eval==NULL))	GMERR(-1);
	n = m->size1;
	if((m->size1!=m->size2))						GMERR(-11);
	if((mcpy->size1!=mcpy->size2))					GMERR(-21);
	if((m->size1!=mcpy->size1))						GMERR(-31);
	if(eval->size!=n)								GMERR(-41);
	if(gsl_eigen_symm(mcpy,eval,w))					GMERR(-51);


	return(0);
GMERRH("getSEigenvalues",1);
}

int getSEigenvectors(gsl_matrix *m,gsl_matrix *mcpy,gsl_eigen_symmv_workspace *w,gsl_vector *eval,gsl_matrix *evec){
	if((m==NULL)||(mcpy==NULL)||(w==NULL)||(eval==NULL)||(evec==NULL))	GMERR(-1);
	if(gsl_eigen_symmv(mcpy,eval,evec,w))			GMERR(-11);
	if(gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC))	GMERR(-21);
	return(0);
GMERRH("getSEigenvectors",1);
}


int tgetEigenvalues(gsl_matrix *m,gsl_vector_complex *v){
	int n;
	gsl_eigen_nonsymm_workspace *w;
	gsl_matrix *mcpy = gsl_matrix_alloc(m->size1,m->size2);
	if((m==NULL)||(v==NULL))						GMERR(-1);
	n = m->size1;
	if((m->size2!=n)||(v->size!=n))					GMERR(-11);
	w = gsl_eigen_nonsymm_alloc(n);
	gsl_matrix_memcpy(mcpy,m);
	gsl_eigen_nonsymm(mcpy,v,w);

	gsl_matrix_free(mcpy);
	gsl_eigen_nonsymm_free(w);
	
	return(0);
GMERRH("tgetEigenvalues",1);
}

int tgetEigenvectors(gsl_matrix *m,gsl_vector_complex *eigenvalues,gsl_matrix_complex *eigenvectors){
	int n;
	gsl_eigen_nonsymmv_workspace *w;
	gsl_matrix *mcpy = gsl_matrix_alloc(m->size1,m->size2);
	if((m==NULL)||(eigenvalues==NULL)||(eigenvectors==NULL))	GMERR(-1);
	n = m->size1;
	if((m->size2!=n)||(eigenvalues->size!=n))					GMERR(-11);
	if((eigenvectors->size1!=n)||(eigenvectors->size2!=n))		GMERR(-21);
	w = gsl_eigen_nonsymmv_alloc(n);
	gsl_matrix_memcpy(mcpy,m);
	gsl_eigen_nonsymmv(mcpy,eigenvalues,eigenvectors,w);
	gsl_eigen_nonsymmv_sort(eigenvalues,eigenvectors,GSL_EIGEN_SORT_ABS_ASC);

	gsl_matrix_free(mcpy);
	gsl_eigen_nonsymmv_free(w);

	return(0);
GMERRH("tgetEigenvectors",1);
}

int getEigenvalues(gsl_matrix *m,gsl_matrix *mcpy,gsl_eigen_nonsymm_workspace *w,gsl_vector_complex *eval){
	if((m==NULL)||(mcpy==NULL)||(w==NULL)||(eval==NULL))		GMERR(-1);
	gsl_matrix_memcpy(mcpy,m);
	gsl_eigen_nonsymm(mcpy,eval,w);
	return(0);
GMERRH("getEigenvalues",1);
}

int getEigenvectors(gsl_matrix *m,gsl_matrix *mcpy,gsl_eigen_nonsymmv_workspace *w,gsl_vector_complex *eval,gsl_matrix_complex *evec){
	if((m==NULL)||(mcpy==NULL)||(w==NULL)||(eval==NULL)||(evec==NULL))	GMERR(-1);
	gsl_matrix_memcpy(mcpy,m);
	gsl_eigen_nonsymmv(mcpy,eval,evec,w);
	gsl_eigen_nonsymmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);
	return(0);
GMERRH("getEigenvectors",1);
}

