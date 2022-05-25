#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "gmlib.h"
#include "matrix_utils.h"
#include "mdarray.h"


/*
Observation matrix for state space model. In this case, we assume that the
Np components of the observation vector each track the state of the same 
index (more specifically, its U component). All other terms are ignored.
Np is the dimension of the observation vector.
Ns is the number of underlying states.
*/

//Works
int NGP_Z(int Np,int Ns,marray3d *Z){
	int minInd,i;
	//if(marrayCheckYZ(Z,Np,3*Ns))					GMERR(-11);
	if(marray3d_checkYZ(Z,Np,3*Ns))					GMERR(-11);
	//if((Z->size1!=Np)||(Z->size2!=(3*Ns)))				GMERR(-11);
	//marray_set_zero(Z);
	marray3d_set_zero(Z);
	if(Np>3*Ns){
		minInd = 3*Ns;
	}else{
		minInd = Np;
	}
	for(i=0;i<minInd;i++){
		marray3d_set(Z,0,i,i,1.0);
	}
	return(0);
GMERRH("NGP2_Z",1);
}

/*
State transition matrix. Defined in (5) of Zhu and Dunson. Order of entries
is all state variables, followed by all derivatives, then all latents.
Ns is the dimension of the underlying state space.
*/

//Works
int NGP_G(int Ns,double delta,int approx,gsl_matrix *G){
	int i;
	if((G->size1!=(3*Ns))||(G->size2!=(3*Ns)))			GMERR(-11);
	gsl_matrix_set_identity(G);
	for(i=0;i<2*Ns;i++){
		gsl_matrix_set(G,i,i+Ns,delta);
	}
	if(approx){
		for(i=0;i<Ns;i++){
			gsl_matrix_set(G,i,i+2*Ns,delta*delta/2);
		}
	}
	return(0);
GMERRH("NGP_G2",1);
}

/*
Covariance of state disturbances. Defined in (5) of Zhu and Dunson. Order of
entries is all state variables, followed by all derivatives, then all latents.
*/

//Works
int NGP_W(int Ns,double delta,gsl_vector *varU,gsl_vector *varA,int approx,gsl_matrix *W){
	int indx;
	double del2=delta*delta;
	double del3=delta*delta*delta;
	double del4=delta*delta*delta*delta;
	double del5=delta*delta*delta*delta*delta;
	double varUi,varAi;
	int count=0;
	if((varU->size!=Ns)&&(varU->size!=1))					GMERR(-1);
	if((varA->size!=Ns)&&(varA->size!=1))					GMERR(-2);
	if(approx){
		if((W->size1!=(2*Ns))||(W->size2!=(2*Ns)))			GMERR(-11);
		gsl_matrix_set_zero(W);
		for(indx=0;indx<Ns;indx++){
			gsl_matrix_set(W,indx,indx,delta*gsl_vector_get(varU,indx));
		}
		for(indx=Ns;indx<(2*Ns);indx++){
			gsl_matrix_set(W,indx,indx,delta*gsl_vector_get(varA,count));
			count = count + 1;
		}
	}else{
		if((W->size1!=(3*Ns))||(W->size2!=(3*Ns)))			GMERR(-21);
		gsl_matrix_set_zero(W);
		for(indx=0;indx<Ns;indx++){
			varUi = gsl_vector_get(varU,indx);
			varAi = gsl_vector_get(varA,indx);
			gsl_matrix_set(W,indx,indx,(del3/3)*varUi+(del5/20)*varAi);
			gsl_matrix_set(W,indx,Ns+indx,(del2/2)*varUi+(del4/8)*varAi);
			gsl_matrix_set(W,indx,2*Ns+indx,(del3/6)*varAi);
			gsl_matrix_set(W,Ns+indx,indx,(del2/2)*varUi+(del4/8)*varAi);
			gsl_matrix_set(W,Ns+indx,Ns+indx,delta*varUi+(del3/3)*varAi);
			gsl_matrix_set(W,Ns+indx,2*Ns+indx,(del2/2)*varAi);
			gsl_matrix_set(W,2*Ns+indx,indx,(del3/6)*varAi);
			gsl_matrix_set(W,2*Ns+indx,Ns+indx,(del2/2)*varAi);
			gsl_matrix_set(W,2*Ns+indx,2*Ns+indx,delta*varAi);
		}
	}
	return(0);
GMERRH("getW2",1);
}

/*
Take nested GP parameters and return matrices suitable for feeding into
state space model.
Np is the dimension of the observation vector
Ns is the number of underlying states
δ is a vector of differences between observation times. If observations are 
uniformly spaced, δ can be a scalar or a length-1 vector.
*/


static int _checkSizes(int Np,int Ns,gsl_vector *delta,double *sigEps,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_vector *a0,gsl_matrix *P0){
	int Nt = delta->size;
	int Nm = 3*Ns;
	int Nr;
	if(approx){
		Nr=2*Ns;
	}else{
		Nr=3*Ns;
	}
	if(marray3d_checkYZ(H,Np,Np))								GMERR(-11);
	if(marray3d_checkSizes(T,Nt,Nm,Nm))							GMERR(-21);
	if(marray3d_checkSizes(Q,Nt,Nr,Nr))							GMERR(-31);
	if(marray3d_checkYZ(R,Nm,Nr))								GMERR(-41);
	if(a0->size!=Nm)											GMERR(-51);
	if((P0->size1!=Nm)||(P0->size2!=Nm))						GMERR(-61);
	if(marray3d_checkYZ(Z,Np,Nm))								GMERR(-71);
	
	return(0);
GMERRH("_checkSizes",1);
}

int NGP_assemble_matrices(int Np,int Ns,gsl_vector *delta,double *sigEps,gsl_vector *sigU,gsl_vector *sigA,double *sigMu,double *sigAlph,int approx,marray3d *Z,marray3d *H,marray3d *T,marray3d *R,marray3d *Q,gsl_vector *a0,gsl_matrix *P0){
	//calculate dimensions
	//Note: for an observation vector of size Np, we introduce 3 latent
	//variables per observation and 2 or 3 noise variables (depending
	//on approximation)
	int Nt = delta->size;
	int Nr,t;
	int Nm = 3*Ns;
	gsl_matrix *RR,*TT,*QQ,*HH;
	gsl_vector *sigU2 = gsl_vector_alloc(sigU->size);
	gsl_vector *sigA2 = gsl_vector_alloc(sigA->size);
	if(approx)
		Nr = 2*Ns;
	else
		Nr = 3*Ns;
	if(_checkSizes(Np,Ns,delta,sigEps,sigU,sigA,sigMu,sigAlph,approx,Z,H,T,R,Q,a0,P0))		GMERR(-11);
	RR = gsl_matrix_alloc(Nm,Nr);
	HH = gsl_matrix_alloc(Np,Np);
	TT = gsl_matrix_alloc(Nm,Nm);
	QQ = gsl_matrix_alloc(Nr,Nr);
	sigU2 = gsl_vector_alloc(sigU->size);
	sigA2 = gsl_vector_alloc(sigA->size);
	gsl_vector_memcpy(sigU2,sigU);
	gsl_vector_mul(sigU2,sigU);
	gsl_vector_memcpy(sigA2,sigA);
	gsl_vector_mul(sigA2,sigA);

	gsl_matrix_set_identity(HH);
	gsl_matrix_scale(HH,*sigEps);
	gsl_matrix_scale(HH,*sigEps);
	marray3d_set_X(H,HH,0);

	if(approx){
		//R = sparse([2:3:Nm ; 3:3:Nm], [1:2:Nr ; 2:2:Nr], 1., Nm, Nr)
		gsl_matrix_set_zero(RR);
		for(t=0;t<Ns;t++){
			gsl_matrix_set(RR,3*t+1,2*t,1.0);
			gsl_matrix_set(RR,3*t+2,2*t+1,1.0);
		}
		if(marray3d_set_X(R,RR,0))							GMERR(-31);

	}else{
		//R = speye(Float64, Nm)
		gsl_matrix_set_identity(RR);
		if(marray3d_set_X(R,RR,0))							GMERR(-31);
	}
	
	for(t=0;t<Nt;t++){
		//T[:,:,t] = G(Ns,delta[t],approx=approx)
		if(NGP_G(Ns,gsl_vector_get(delta,t),approx,TT))		GMERR(-41);
		if(marray3d_set_X(T,TT,t))							GMERR(-51);
		//Q[:,:,t] = W(Ns,delta[t],sigU.^2,sigA.^2,approx =approx)
		if(NGP_W(Ns,gsl_vector_get(delta,t),sigU2,sigA2,approx,QQ))		GMERR(-61);
		if(marray3d_set_X(Q,QQ,t))							GMERR(-71);

	}

	//a0 = zeros(Nm)
	gsl_vector_set_zero(a0);
	gsl_matrix_set_zero(P0);

	//P0 = diagm(repeat([σμ, σμ, σα], inner=[Ns]))
	for(t=0;t<Ns;t++){
		gsl_matrix_set(P0,3*t,3*t,*sigMu);
		gsl_matrix_set(P0,3*t+1,3*t+1,*sigMu);
		gsl_matrix_set(P0,3*t+2,3*t+2,*sigAlph);
	}

	//Z=Z(Np,Ns)
	if(NGP_Z(Np,Ns,Z))										GMERR(-71);

	gsl_vector_free(sigA2);
	gsl_vector_free(sigU2);

	return(0);
GMERRH("NGP2_assemble_matrices",1);
}

