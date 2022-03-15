
/*
 (rho v, u) = ( v , div(C : grad(u)) ) + (v,F)
 
 where (v,F) = (v, Mp delta_1) + (v, -div(-Ms delta_2 ))
 
 (rho v, u) = -(grad(v) , C : grad(u)) + (v, sigma.n)_ds
              + (v, Mp delta_1) + (grad(v), Ms delta_2 ) - (v, (Ms delta_2 ).n)_ds

            = -(grad(v) , C : grad(u))
              +(v, Mp delta_1 ) + (grad(v), Ms delta_2 )
              + (v, sigma.n)_ds
              - (v, (-Ms delta_2 ).n)_ds  ==> which will be dropped as we assume the source does not intersect with the boundary
 
 
 
*/

#include <petsc.h>
#include <petsctime.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>

#include "rupture.h"

typedef enum { TENS2D_XX=0, TENS2D_YY=1, TENS2D_XY=2 } VoigtTensor2d;

typedef struct _p_SpecFECtx *SpecFECtx;

typedef struct {
  PetscInt region;
  PetscReal lambda,mu;
  PetscReal rho;
} QPntIsotropicElastic;

struct _p_SpecFECtx {
  PetscMPIInt rank,size;
  PetscInt basisorder;
  PetscInt mx,my,mz;
  PetscInt mx_g,my_g,mz_g,nx_g,ny_g,nz_g;
  //PetscReal dx,dy,dz;
  PetscInt dim;
  PetscInt dofs;
  DM dm;
  PetscInt npe,npe_1d,ne,ne_g;
  PetscInt *element;
  PetscReal *xi1d,*w1d,*w;
  PetscReal *elbuf_coor,*elbuf_field,*elbuf_field2;
  PetscInt  *elbuf_dofs;
  PetscInt nqp;
  QPntIsotropicElastic *cell_data;
  PetscReal **dN_dxi,**dN_deta;
  PetscReal **dN_dx,**dN_dy;
  //PetscInt  source_implementation; /* DR */
  
  /* DR */
  PetscReal delta;         /* fault thickness */
  PetscReal *elbuf_field3; /* additional element buffer (will hold velocity) */
  DRVar     *dr_qp_data;   /* stores data like slip, slip-rate */
  PetscReal mu_s,mu_d,D_c; /* linear slip weakening parameters */
};


typedef struct {
  PetscReal xi[2];
  PetscInt  nbasis;
  PetscInt  *element_indices;
  PetscReal *element_values;
  PetscReal *buffer;
} PointwiseContext;




/**
 * Function to calculate weighting for the traction
*/ 
PetscErrorCode PetscTanHWeighting(PetscReal *Result, PetscReal ValueTrial, PetscReal CritValue,  PetscReal phi, PetscReal Amplitude, PetscReal Offset)
{
  PetscReal weight;

  weight = 0.5 * PetscTanhReal((PetscAbsReal(phi)-Offset) * Amplitude)  + 0.5;

  Result[0] =  CritValue * (1.0 - weight)  + ValueTrial * weight ;
  PetscFunctionReturn(0);
}

/**
 * Function to calculate the new KV timestep following Galvez (2014) eq 27.
 * dT_KV = (sqrt(1+(eta*eta)/(dT*dT))-eta/dT)
*/ 
void GetStableTimeStep(double dT, double eta, double * dT_KV)
{
  dT_KV[0] = (sqrt(1 + (eta/dT) * (eta/dT)) - eta / dT)*dT;
} 

/*
 warp for dr mesh
 
 get ymax
 
 plot (exp(4*x)-1)/exp(4),x
 
 s = y / ymax
 s' = (exp(4*s)-1)/exp(4)
 
*/
PetscErrorCode warp_y_exp(SpecFECtx c,PetscReal factor)
{
  PetscInt i,N;
  PetscReal ymax = -1.0e32,s[2],sp[2];
  Vec coor;
  PetscScalar *_coor;
  PetscErrorCode ierr;
  
  DMGetCoordinates(c->dm,&coor);
  VecGetSize(coor,&N);
  N = N / 2;
  VecGetArray(coor,&_coor);
  for (i=0; i<N; i++) {
    ymax = PetscMax(ymax,_coor[2*i+1]);
  }

  for (i=0; i<N; i++) {
    s[0] = _coor[2*i+0];
    s[1] = _coor[2*i+1];

    // normalize to 1
    s[1] = s[1] / ymax;

    sp[0] = s[0];
    sp[1] = s[1];
    if (s[1] >= 0.0) {
      sp[1] = (PetscExpReal(factor * s[1]) - 1.0)/PetscExpReal(factor);
    } else {
      PetscReal _s = PetscAbsReal(s[1]);
      
      sp[1] = -(PetscExpReal(factor * _s) - 1.0)/PetscExpReal(factor);
    }

    sp[1] *= ymax;
   
    _coor[2*i+0] = sp[0];
    _coor[2*i+1] = sp[1];
  }
  
  VecRestoreArray(coor,&_coor);
  
  PetscFunctionReturn(0);
}


/* N = polynomial order */
PetscErrorCode CreateGLLCoordsWeights(PetscInt N,PetscInt *_npoints,PetscReal **_xi,PetscReal **_w)
{
  PetscInt N1;
  PetscReal *xold,*x,*w,*P;
  PetscReal eps,res;
  PetscInt i,j,k;
  PetscErrorCode ierr;
  
  
  // Truncation + 1
  N1 = N + 1;
  
  ierr = PetscMalloc(sizeof(PetscReal)*N1,&xold);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*N1,&x);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*N1,&w);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*N1*N1,&P);CHKERRQ(ierr);
  
  // Use the Chebyshev-Gauss-Lobatto nodes as the first guess
  for (i=0; i<N1; i++) {
    x[i]=PetscCosReal(PETSC_PI*i/(PetscReal)N);
  }
  
  // The Legendre Vandermonde Matrix
  for (i=0; i<N1; i++) {
    for (j=0; j<N1; j++) {
      P[i+j*N1] = 0.0;
    }
  }
  
  // Compute P_(N) using the recursion relation
  // Compute its first and second derivatives and
  // update x using the Newton-Raphson method.
  for (i=0; i<N1; i++) {
    xold[i]=2.0;
  }
  
  res = 1.0;
  eps = 1.0e-12;
  while (res > eps) {
    
    //xold=x;
    for (i=0; i<N1; i++) {
      xold[i] = x[i];
    }
    
    //P(:,1)=1;    P(:,2)=x;
    for (i=0; i<N1; i++) {
      for (j=0; j<N1; j++) {
        P[i+0*N1] = 1.0;
        P[i+1*N1] = x[i];
      }
    }
    
    //for k=2:N
    //    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    //end
    for (i=0; i<N1; i++) {
      for (k=1; k<N; k++) {
        P[i+(k+1)*N1] = ( (2.0*(k+1)-1.0)*x[i] * P[i+k*N1] - (k+1.0-1.0) * P[i+(k-1)*N1] ) / (PetscReal)(k+1.0);
      }
    }
    
    //x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
    for (i=0; i<N1; i++) {
      x[i] = xold[i] - (x[i] * P[i+(N1-1)*N1] - P[i+(N-1)*N1]) / ( N1 * P[i+(N1-1)*N1] );
    }
    
    res = 0.0;
    for (i=0; i<N1; i++) {
      res += (x[i] - xold[i])*(x[i] - xold[i]);
    }
    res = PetscSqrtReal(res);
  }
  
  // w=2./(N*N1*P(:,N1).^2);
  for (i=0; i<N1; i++) {
    PetscReal pp = P[i+(N1-1)*N1];
    w[i] = 2.0 / (N*N1*pp*pp);
  }
  
  if (_xi) {
    /* flip order so they are ordered from -1 to 1 */
    for (i=0; i<N1/2; i++) {
      PetscReal tmp;
      
      tmp = x[i];
      x[i] = x[N1-1-i];
      x[N1-1-i] = tmp;
    }
    *_xi = x;
  } else {
    ierr = PetscFree(x);CHKERRQ(ierr);
  }
  
  if (_npoints) {
    *_npoints = N1;
  }
  if (_w) {
    *_w = w;
  } else {
    ierr = PetscFree(w);CHKERRQ(ierr);
  }
  ierr = PetscFree(xold);CHKERRQ(ierr);
  ierr = PetscFree(P);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode MatComputeConditionNumber(Mat A,PetscReal *cond)
{
  PetscReal *realpt,*complexpt,*nrmeigs;
  PetscInt rank,i;
  KSP kspV;
  PC pc;
  Vec x,y;
  PetscErrorCode ierr;
  
  ierr = MatCreateVecs(A,&y,&x);CHKERRQ(ierr);
  ierr = VecSet(y,1.0);CHKERRQ(ierr);
  
  ierr = KSPCreate(PETSC_COMM_SELF,&kspV);CHKERRQ(ierr);
  ierr = KSPSetOperators(kspV,A,A);CHKERRQ(ierr);
  ierr = KSPSetType(kspV,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(kspV,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
  ierr = KSPSolve(kspV,y,x);CHKERRQ(ierr);
  
  ierr = MatGetSize(A,&rank,0);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*rank,&realpt);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*rank,&complexpt);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*rank,&nrmeigs);CHKERRQ(ierr);
  
  ierr = KSPComputeEigenvaluesExplicitly(kspV,rank,realpt,complexpt);CHKERRQ(ierr);
  for (i=0; i<rank; i++) {
    nrmeigs[i] = PetscSqrtReal( realpt[i]*realpt[i] + complexpt[i]*complexpt[i]);
  }
  ierr = PetscSortReal(rank,nrmeigs);CHKERRQ(ierr);
  
  *cond = nrmeigs[rank-1]/nrmeigs[0];
  
  ierr = PetscFree(nrmeigs);CHKERRQ(ierr);
  ierr = PetscFree(realpt);CHKERRQ(ierr);
  ierr = PetscFree(complexpt);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = KSPDestroy(&kspV);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode TabulateBasis1d_CLEGENDRE(PetscInt npoints,PetscReal xi[],PetscInt order,PetscInt *_nbasis,PetscReal ***_Ni)
{
  PetscErrorCode ierr;
  PetscReal **Ni,*xilocal,**basis_coeff, *monomials;
  PetscInt i,j,k,p;
  PetscInt nbasis,cnt;
  Mat A;
  Vec x,y;
  KSP ksp;
  PC pc;
  
  
  ierr = CreateGLLCoordsWeights(order,&nbasis,&xilocal,NULL);CHKERRQ(ierr);
  
  ierr = PetscMalloc(sizeof(PetscReal)*nbasis,&monomials);CHKERRQ(ierr);
  
  ierr = PetscMalloc(sizeof(PetscReal*)*npoints,&Ni);CHKERRQ(ierr);
  for (i=0; i<npoints; i++) {
    ierr = PetscMalloc(sizeof(PetscReal)*nbasis,&Ni[i]);CHKERRQ(ierr);
  }
  
  ierr = PetscMalloc(sizeof(PetscReal*)*nbasis,&basis_coeff);CHKERRQ(ierr);
  for (i=0; i<nbasis; i++) {
    ierr = PetscMalloc(sizeof(PetscReal)*nbasis,&basis_coeff[i]);CHKERRQ(ierr);
  }
  
  /* generate all the basis coefficients */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,nbasis,nbasis,NULL,&A);CHKERRQ(ierr);
  for (k=0; k<nbasis; k++) {
    PetscReal xil,Aij;
    
    xil  = xilocal[k];
    
    cnt = 0;
    for (i=0; i<nbasis; i++) {
      Aij = PetscPowReal(xil,(PetscReal)i);
      ierr = MatSetValue(A,k,cnt,Aij,INSERT_VALUES);CHKERRQ(ierr);
      cnt++;
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  {
    PetscReal cond;
    PetscBool compute_vandermonde_condition = PETSC_FALSE;
    
    ierr = PetscOptionsGetBool(NULL,NULL,"-compute_vandermonde_condition",&compute_vandermonde_condition,NULL);CHKERRQ(ierr);
    if (compute_vandermonde_condition) {
      
      PetscPrintf(PETSC_COMM_WORLD,"Computing condition number of Vandermonde matrix\n");
      ierr = MatComputeConditionNumber(A,&cond);CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"cond(V) = %1.6e \n",cond);
    }
  }
  
  ierr = MatCreateVecs(A,&x,&y);CHKERRQ(ierr);
  
  ierr = KSPCreate(PETSC_COMM_SELF,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(ksp,"basis_");CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  
  for (k=0; k<nbasis; k++) {
    const PetscScalar *LA_x;
    
    ierr = VecZeroEntries(y);CHKERRQ(ierr);
    ierr = VecSetValue(y,k,1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(y);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y);CHKERRQ(ierr);
    
    ierr = KSPSolve(ksp,y,x);CHKERRQ(ierr);
    
    ierr = VecGetArrayRead(x,&LA_x);CHKERRQ(ierr);
    for (i=0; i<nbasis; i++) {
      basis_coeff[k][i] = LA_x[i];
    }
    ierr = VecRestoreArrayRead(x,&LA_x);CHKERRQ(ierr);
  }
  
  /* evaluate basis at each xi[] */
  for (p=0; p<npoints; p++) {
    
    /* generate all monomials for point, p */
    cnt = 0;
    for (i=0; i<nbasis; i++) {
      monomials[cnt] = PetscPowReal((PetscReal)xi[p],(PetscReal)i);
      cnt++;
    }
    
    for (i=0; i<nbasis; i++) {
      Ni[p][i] = 0.0;
      
      for (j=0; j<nbasis; j++) {
        Ni[p][i] += basis_coeff[i][j] * monomials[j];
      }
      if (PetscAbsReal(Ni[p][i]) < 1.0e-12) {
        Ni[p][i] = 0.0;
      }
    }
    
    /*
     printf("p = %d (xi = %+1.4e) N = [",p,xi[p]);
     for (i=0; i<nbasis; i++) {
     printf(" %+1.4e ",Ni[p][i]);
     }
     printf("]\n");
     */
  }
  
  //for (p=0; p<npoints; p++) {
  //  ierr = PetscFree(Ni[p]);CHKERRQ(ierr);
  //}
  //ierr = PetscFree(Ni);CHKERRQ(ierr);
  *_Ni = Ni;
  *_nbasis = nbasis;
  
  ierr = PetscFree(monomials);CHKERRQ(ierr);
  for (i=0; i<nbasis; i++) {
    ierr = PetscFree(basis_coeff[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(basis_coeff);CHKERRQ(ierr);
  ierr = PetscFree(xilocal);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode TabulateBasisDerivatives1d_CLEGENDRE(PetscInt npoints,PetscReal xi[],PetscInt order,PetscInt *_nbasis,PetscReal ***_GNix)
{
  PetscErrorCode ierr;
  PetscReal **GNix,*xilocal,**basis_coeff, *monomials;
  PetscInt i,j,k,p;
  PetscInt nbasis,cnt;
  Mat A;
  Vec x,y;
  KSP ksp;
  PC pc;
  
  
  ierr = CreateGLLCoordsWeights(order,&nbasis,&xilocal,NULL);CHKERRQ(ierr);
  
  ierr = PetscMalloc(sizeof(PetscReal)*nbasis,&monomials);CHKERRQ(ierr);
  
  ierr = PetscMalloc(sizeof(PetscReal*)*npoints,&GNix);CHKERRQ(ierr);
  for (i=0; i<npoints; i++) {
    ierr = PetscMalloc(sizeof(PetscReal)*nbasis,&GNix[i]);CHKERRQ(ierr);
  }
  
  ierr = PetscMalloc(sizeof(PetscReal*)*nbasis,&basis_coeff);CHKERRQ(ierr);
  for (i=0; i<nbasis; i++) {
    ierr = PetscMalloc(sizeof(PetscReal)*nbasis,&basis_coeff[i]);CHKERRQ(ierr);
  }
  
  /* generate all the basis coefficients */
  ierr = MatCreateSeqDense(PETSC_COMM_SELF,nbasis,nbasis,NULL,&A);CHKERRQ(ierr);
  for (k=0; k<nbasis; k++) {
    PetscReal xil,Aij;
    
    xil  = xilocal[k];
    
    cnt = 0;
    for (i=0; i<nbasis; i++) {
      Aij = PetscPowReal(xil,(PetscReal)i);
      ierr = MatSetValue(A,k,cnt,Aij,INSERT_VALUES);CHKERRQ(ierr);
      cnt++;
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  {
    PetscReal cond;
    PetscBool compute_vandermonde_condition = PETSC_FALSE;
    
    ierr = PetscOptionsGetBool(NULL,NULL,"-compute_vandermonde_condition",&compute_vandermonde_condition,NULL);CHKERRQ(ierr);
    if (compute_vandermonde_condition) {
      
      PetscPrintf(PETSC_COMM_WORLD,"Computing condition number of Vandermonde matrix\n");
      ierr = MatComputeConditionNumber(A,&cond);CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"cond(V) = %1.6e \n",cond);
    }
  }
  
  ierr = MatCreateVecs(A,&x,&y);CHKERRQ(ierr);
  
  ierr = KSPCreate(PETSC_COMM_SELF,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(ksp,"basis_");CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  
  for (k=0; k<nbasis; k++) {
    const PetscScalar *LA_x;
    
    ierr = VecZeroEntries(y);CHKERRQ(ierr);
    ierr = VecSetValue(y,k,1.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(y);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y);CHKERRQ(ierr);
    
    ierr = KSPSolve(ksp,y,x);CHKERRQ(ierr);
    
    ierr = VecGetArrayRead(x,&LA_x);CHKERRQ(ierr);
    for (i=0; i<nbasis; i++) {
      basis_coeff[k][i] = LA_x[i];
    }
    ierr = VecRestoreArrayRead(x,&LA_x);CHKERRQ(ierr);
  }
  
  /* evaluate basis at each xi[] */
  for (p=0; p<npoints; p++) {
    
    /* generate all monomials for point, p */
    cnt = 0;
    for (i=0; i<nbasis; i++) {
      PetscReal dm_dx;
      
      if (i == 0) {
        dm_dx = 0.0;
        } else {
          dm_dx = ((PetscReal)i)*PetscPowReal((PetscReal)xi[p],(PetscReal)(i-1));
        }
      
      monomials[cnt] = dm_dx;
      cnt++;
    }
    
    for (i=0; i<nbasis; i++) {
      GNix[p][i] = 0.0;
      
      for (j=0; j<nbasis; j++) {
        GNix[p][i] += basis_coeff[i][j] * monomials[j];
      }
      if (PetscAbsReal(GNix[p][i]) < 1.0e-12) {
        GNix[p][i] = 0.0;
      }
    }
    
    /*
     printf("p = %d (xi = %+1.4e) dN_dx = [",p,xi[p]);
     for (i=0; i<nbasis; i++) {
     printf(" %+1.4e ",GNix[p][i]);
     }
     printf("]\n");
     */
  }
  
  ierr = PetscFree(monomials);CHKERRQ(ierr);
  for (i=0; i<nbasis; i++) {
    ierr = PetscFree(basis_coeff[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(basis_coeff);CHKERRQ(ierr);
  ierr = PetscFree(xilocal);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  
  *_nbasis = nbasis;
  *_GNix = GNix;
  
  PetscFunctionReturn(0);
}

PetscErrorCode TabulateBasisDerivativesTensorProduct2d(PetscInt order,PetscReal ***_dN_dxi,PetscReal ***_dN_deta)
{
  PetscErrorCode ierr;
  PetscReal *xiq,**dphi_xi,**dN_dxi,**dN_deta;
  PetscInt qpoint,k,i,j,qi,qj,nqp,nbasis;
  
  ierr = CreateGLLCoordsWeights(order,&nqp,&xiq,NULL);CHKERRQ(ierr);
  ierr = TabulateBasisDerivatives1d_CLEGENDRE(nqp,xiq,order,&nbasis,&dphi_xi);CHKERRQ(ierr);
  
  ierr = PetscMalloc(sizeof(PetscReal*)*nqp*nqp,&dN_dxi);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal*)*nqp*nqp,&dN_deta);CHKERRQ(ierr);
  for (i=0; i<nqp*nqp; i++) {
    ierr = PetscMalloc(sizeof(PetscReal)*nbasis*nbasis,&dN_dxi[i]);CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal)*nbasis*nbasis,&dN_deta[i]);CHKERRQ(ierr);
  }
  
  qpoint = 0;
  for (qj=0; qj<nqp; qj++) {
    for (qi=0; qi<nqp; qi++) {
      
      k = 0;
      for (j=0; j<nbasis; j++) {
        for (i=0; i<nbasis; i++) {
          PetscReal phi_xi,phi_eta;
          
          phi_xi = 0.0;
          if (qi == i) phi_xi = 1.0;
          
          phi_eta = 0.0;
          if (qj == j) phi_eta = 1.0;
          
          dN_dxi[qpoint][k]  = dphi_xi[qi][i] * phi_eta;
          dN_deta[qpoint][k] = phi_xi * dphi_xi[qj][j];
          
          k++;
        }}
      qpoint++;
    }}
  
  /* viewer */
  /*
   for (k=0; k<nqp*nqp; k++) {
   printf("qp[%d]: dNdxi  = [ ",k);
   for (j=0; j<nbasis*nbasis; j++) {
   printf(" %+1.4e ",dN_dxi[k][j]);
   } printf("]\n");
   
   printf("qp[%d]: dNdeta = [ ",k);
   for (j=0; j<nbasis*nbasis; j++) {
   printf(" %+1.4e ",dN_deta[k][j]);
   } printf("]\n");
   }
   */
  
  /* free up mempry */
  ierr = PetscFree(xiq);CHKERRQ(ierr);
  for (k=0; k<nqp; k++) {
    ierr = PetscFree(dphi_xi[k]);CHKERRQ(ierr);
  }
  ierr = PetscFree(dphi_xi);CHKERRQ(ierr);
  
  if (_dN_dxi) { *_dN_dxi = dN_dxi; }
  else {
    for (k=0; k<nqp*nqp; k++) {
      ierr = PetscFree(dN_dxi[k]);CHKERRQ(ierr);
    }
    ierr = PetscFree(dN_dxi);CHKERRQ(ierr);
  }
  
  if (_dN_deta) { *_dN_deta = dN_deta; }
  else {
    for (k=0; k<nqp*nqp; k++) {
      ierr = PetscFree(dN_deta[k]);CHKERRQ(ierr);
    }
    ierr = PetscFree(dN_deta);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode TabulateBasisDerivativesAtPointTensorProduct2d(PetscReal xiq[],PetscInt order,PetscReal ***_dN_dxi,PetscReal ***_dN_deta)
{
  PetscErrorCode ierr;
  PetscReal **dphi_xi,**Ni_xi,**dphi_eta,**Ni_eta,**dN_dxi,**dN_deta;
  PetscInt qpoint,k,i,j,q,nqp,nbasis;
  
  nqp = 1;
  ierr = TabulateBasisDerivatives1d_CLEGENDRE(nqp,&xiq[0],order,&nbasis,&dphi_xi);CHKERRQ(ierr);
  ierr = TabulateBasis1d_CLEGENDRE(nqp,&xiq[0],order,&nbasis,&Ni_xi);CHKERRQ(ierr);
  
  ierr = TabulateBasisDerivatives1d_CLEGENDRE(nqp,&xiq[1],order,&nbasis,&dphi_eta);CHKERRQ(ierr);
  ierr = TabulateBasis1d_CLEGENDRE(nqp,&xiq[1],order,&nbasis,&Ni_eta);CHKERRQ(ierr);
  
  ierr = PetscMalloc(sizeof(PetscReal*)*nqp,&dN_dxi);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal*)*nqp,&dN_deta);CHKERRQ(ierr);
  for (i=0; i<nqp; i++) {
    ierr = PetscMalloc(sizeof(PetscReal)*nbasis*nbasis,&dN_dxi[i]);CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal)*nbasis*nbasis,&dN_deta[i]);CHKERRQ(ierr);
  }
  
  qpoint = 0;
  for (q=0; q<nqp; q++) {
    
    k = 0;
    for (j=0; j<nbasis; j++) {
      for (i=0; i<nbasis; i++) {
        PetscReal phi_xi,phi_eta;
        
        phi_xi = Ni_xi[q][i];
        phi_eta = Ni_eta[q][j];
        
        dN_dxi[qpoint][k]  = dphi_xi[q][i] * phi_eta;
        dN_deta[qpoint][k] = phi_xi * dphi_eta[q][j];
        k++;
      }}
    qpoint++;
  }
  
  /* viewer */
  /*
   for (k=0; k<nqp; k++) {
   printf("qp[%d]: dNdxi  = [ ",k);
   for (j=0; j<nbasis*nbasis; j++) {
   printf(" %+1.4e ",dN_dxi[k][j]);
   } printf("]\n");
   
   printf("qp[%d]: dNdeta = [ ",k);
   for (j=0; j<nbasis*nbasis; j++) {
   printf(" %+1.4e ",dN_deta[k][j]);
   } printf("]\n");
   }
   */
  
  /* free up mempry */
  for (k=0; k<nqp; k++) {
    ierr = PetscFree(dphi_xi[k]);CHKERRQ(ierr);
    ierr = PetscFree(dphi_eta[k]);CHKERRQ(ierr);
  }
  ierr = PetscFree(dphi_xi);CHKERRQ(ierr);
  ierr = PetscFree(dphi_eta);CHKERRQ(ierr);
  
  for (k=0; k<nqp; k++) {
    ierr = PetscFree(Ni_xi[k]);CHKERRQ(ierr);
    ierr = PetscFree(Ni_eta[k]);CHKERRQ(ierr);
  }
  ierr = PetscFree(Ni_xi);CHKERRQ(ierr);
  ierr = PetscFree(Ni_eta);CHKERRQ(ierr);
  
  if (_dN_dxi) { *_dN_dxi = dN_dxi; }
  else {
    for (k=0; k<nqp*nqp; k++) {
      ierr = PetscFree(dN_dxi[k]);CHKERRQ(ierr);
    }
    ierr = PetscFree(dN_dxi);CHKERRQ(ierr);
  }
  
  if (_dN_deta) { *_dN_deta = dN_deta; }
  else {
    for (k=0; k<nqp*nqp; k++) {
      ierr = PetscFree(dN_deta[k]);CHKERRQ(ierr);
    }
    ierr = PetscFree(dN_deta);CHKERRQ(ierr);
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxCreate(SpecFECtx *c)
{
  SpecFECtx ctx;
  PetscErrorCode ierr;
  
  ierr = PetscMalloc(sizeof(struct _p_SpecFECtx),&ctx);CHKERRQ(ierr);
  ierr = PetscMemzero(ctx,sizeof(struct _p_SpecFECtx));CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&ctx->rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&ctx->size);CHKERRQ(ierr);
  *c = ctx;
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxCreateENMap2d_SEQ(SpecFECtx c)
{
  PetscErrorCode ierr;
  PetscInt ni0,nj0,i,j,ei,ej,ecnt,*emap,nid;
  
  ierr = PetscMalloc(sizeof(PetscInt)*c->ne*c->npe,&c->element);CHKERRQ(ierr);
  ierr = PetscMemzero(c->element,sizeof(PetscInt)*c->ne*c->npe);CHKERRQ(ierr);
  
  ecnt = 0;
  for (ej=0; ej<c->my; ej++) {
    nj0 = ej*(c->npe_1d-1);
    
    for (ei=0; ei<c->mx; ei++) {
      ni0 = ei*(c->npe_1d-1);
      
      emap = &c->element[c->npe*ecnt];
      
      for (j=0; j<c->npe_1d; j++) {
        for (i=0; i<c->npe_1d; i++) {
          
          nid = (ni0 + i) + (nj0 + j) * c->nx_g;
          emap[i+j*c->npe_1d] = nid;
        }
      }
      
      ecnt++;
    }
  }
  
  PetscFunctionReturn(0);
}

/* Creates domain over [0,1]^d - scale later */
PetscErrorCode SpecFECtxCreateMeshCoords2d_SEQ(SpecFECtx c)
{
  PetscErrorCode ierr;
  Vec coor;
  DM cdm;
  DMDACoor2d **LA_coor2d;
  PetscInt ei,ej,i,j,ni0,nj0;
  PetscReal dx,dy,x0,y0;
  
  ierr = DMDASetUniformCoordinates(c->dm,0.0,1.0,0.0,1.0,0,0);CHKERRQ(ierr);
  ierr = DMGetCoordinates(c->dm,&coor);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(c->dm,&cdm);CHKERRQ(ierr);
  
  dx = 1.0/((PetscReal)c->mx_g);
  dy = 1.0/((PetscReal)c->my_g);
  ierr = DMDAVecGetArray(cdm,coor,&LA_coor2d);CHKERRQ(ierr);
  for (ej=0; ej<c->my; ej++) {
    
    for (ei=0; ei<c->mx; ei++) {
      x0 = 0.0 + ei*dx;
      y0 = 0.0 + ej*dy;
      
      ni0 = ei*(c->npe_1d-1);
      nj0 = ej*(c->npe_1d-1);
      
      for (j=0; j<c->npe_1d; j++) {
        for (i=0; i<c->npe_1d; i++) {
          LA_coor2d[nj0+j][ni0+i].x = 0.5*(c->xi1d[i]+1.0)*dx + x0;
          LA_coor2d[nj0+j][ni0+i].y = 0.5*(c->xi1d[j]+1.0)*dy + y0;
          
          //if ((ej==0) && (j==0)) {
          //  printf("[e %d,i %d] xc %+1.4e\n",ei,i,LA_coor2d[nj0+j][ni0+i].x*4.0e3-2.0e3);
          //}
          
        }
      }
    }
  }
  ierr = DMDAVecRestoreArray(cdm,coor,&LA_coor2d);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxScaleMeshCoords(SpecFECtx c,PetscReal scale[],PetscReal shift[])
{
  PetscErrorCode ierr;
  Vec coor,lcoor;
  DM cdm;
  
  ierr = DMGetCoordinates(c->dm,&coor);CHKERRQ(ierr);
  
  if (scale) {
    if (c->dim >= 1) ierr = VecStrideScale(coor,0,scale[0]);CHKERRQ(ierr);
    if (c->dim >= 2) ierr = VecStrideScale(coor,1,scale[1]);CHKERRQ(ierr);
    if (c->dim == 3) ierr = VecStrideScale(coor,2,scale[2]);CHKERRQ(ierr);
  }
  if (shift) {
    Vec ss;
    
    ierr = VecDuplicate(coor,&ss);CHKERRQ(ierr);
    
    if (c->dim >= 1) {
      ierr = VecZeroEntries(ss);CHKERRQ(ierr);
      ierr = VecStrideSet(ss,0,shift[0]);CHKERRQ(ierr);
      ierr = VecAXPY(coor,1.0,ss);CHKERRQ(ierr);
    }
    if (c->dim >= 2) {
      ierr = VecZeroEntries(ss);CHKERRQ(ierr);
      ierr = VecStrideSet(ss,1,shift[1]);CHKERRQ(ierr);
      ierr = VecAXPY(coor,1.0,ss);CHKERRQ(ierr);
    }
    if (c->dim >= 3) {
      ierr = VecZeroEntries(ss);CHKERRQ(ierr);
      ierr = VecStrideSet(ss,2,shift[2]);CHKERRQ(ierr);
      ierr = VecAXPY(coor,1.0,ss);CHKERRQ(ierr);
    }
    ierr = VecDestroy(&ss);CHKERRQ(ierr);
  }
  
  ierr = DMGetCoordinateDM(c->dm,&cdm);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(c->dm,&lcoor);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(c->dm,coor,INSERT_VALUES,lcoor);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,coor,INSERT_VALUES,lcoor);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxCreateMesh_SEQ(SpecFECtx c,PetscInt dim,PetscInt mx,PetscInt my,PetscInt mz,PetscInt basisorder,PetscInt ndofs)
{
  PetscErrorCode ierr;
  PetscInt stencil_width,i,j;
  
  c->dim = dim;
  c->mx = mx;
  c->my = my;
  c->mz = mz;
  c->mx_g = mx;
  c->my_g = my;
  c->mz_g = mz;
  c->basisorder = basisorder;
  c->dofs = ndofs;
  
  c->nx_g = basisorder*mx + 1;
  c->ny_g = basisorder*my + 1;
  c->nz_g = basisorder*mz + 1;
  
  ierr = CreateGLLCoordsWeights(basisorder,&c->npe_1d,&c->xi1d,&c->w1d);CHKERRQ(ierr);
  
  stencil_width = 1;
  switch (dim) {
    case 2:
    c->npe = c->npe_1d * c->npe_1d;
    c->ne = mx * my;
    c->ne_g = mx * my;
    
    ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
                        c->nx_g,c->ny_g,PETSC_DECIDE,PETSC_DECIDE,ndofs,stencil_width,NULL,NULL,&c->dm);CHKERRQ(ierr);
    ierr = DMSetUp(c->dm);CHKERRQ(ierr);
    ierr = SpecFECtxCreateENMap2d_SEQ(c);CHKERRQ(ierr);
    ierr = SpecFECtxCreateMeshCoords2d_SEQ(c);CHKERRQ(ierr);
    
    /* tensor product for weights */
    ierr = PetscMalloc(sizeof(PetscReal)*c->npe,&c->w);CHKERRQ(ierr);
    for (j=0; j<c->npe_1d; j++) {
      for (i=0; i<c->npe_1d; i++) {
        c->w[i+j*c->npe_1d] = c->w1d[i] * c->w1d[j];
      }
    }
    
    ierr = TabulateBasisDerivativesTensorProduct2d(basisorder,&c->dN_dxi,&c->dN_deta);CHKERRQ(ierr);
    ierr = TabulateBasisDerivativesTensorProduct2d(basisorder,&c->dN_dx,&c->dN_dy);CHKERRQ(ierr);
    
    break;
  }
  
  c->nqp = c->npe;
  
  ierr = PetscMalloc(sizeof(QPntIsotropicElastic)*c->ne,&c->cell_data);CHKERRQ(ierr);
  ierr = PetscMemzero(c->cell_data,sizeof(QPntIsotropicElastic)*c->ne);CHKERRQ(ierr);

  ierr = PetscMalloc1(c->ne * c->nqp,&c->dr_qp_data);CHKERRQ(ierr);
  ierr = PetscMemzero(c->dr_qp_data,sizeof(DRVar)*c->ne*c->nqp);CHKERRQ(ierr);
  
  ierr = PetscMalloc(sizeof(PetscReal)*c->npe*c->dim,&c->elbuf_coor);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*c->npe*c->dofs,&c->elbuf_field);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*c->npe*c->dofs,&c->elbuf_field2);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*c->npe*c->dofs,&c->elbuf_field3);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscInt)*c->npe*c->dofs,&c->elbuf_dofs);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*
 Degree 4 has 5 basis in each direction
 |           |
 0--1--2--3--4
*/
PetscErrorCode SpecFECtxGetCornerBasis_MPI(SpecFECtx c,PetscInt *si,PetscInt *si_g,PetscInt *sj,PetscInt *sj_g)
{
  PetscInt gi,gj,m,n,k;
  PetscErrorCode ierr;

  ierr = DMDAGetGhostCorners(c->dm,&gi,&gj,NULL,&m,&n,NULL);CHKERRQ(ierr);
  /*printf("rank %d: gi,gj %d %d  npe %d\n",c->rank,gi,gj,c->npe_1d);*/
  for (k=0; k<m; k++) {
    if (((gi+k) % (c->npe_1d-1)) == 0) {
      *si = k;
      *si_g = gi+k;
      break;
    }
  }
  for (k=0; k<n; k++) {
    if (((gj+k) % (c->npe_1d-1)) == 0) {
      *sj = k;
      *sj_g = gj+k;
      break;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxGetLocalBoundingBox(SpecFECtx c,PetscReal gmin[],PetscReal gmax[])
{
  PetscErrorCode ierr;
  PetscInt si[]={0,0},si_g[]={0,0},m,n,ii,jj;
  const PetscReal *LA_coor;
  Vec coor;
  
  ierr = SpecFECtxGetCornerBasis_MPI(c,&si[0],&si_g[0],&si[1],&si_g[1]);CHKERRQ(ierr);
  ierr = DMDAGetGhostCorners(c->dm,NULL,NULL,NULL,&m,&n,NULL);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  ii = si[0];
  jj = si[1];
  gmin[0] = LA_coor[2*(ii + jj*m)+0];
  gmin[1] = LA_coor[2*(ii + jj*m)+1];
  ii = si[0] + c->mx * c->basisorder;
  jj = si[1] + c->my * c->basisorder;
  gmax[0] = LA_coor[2*(ii + jj*m)+0];
  gmax[1] = LA_coor[2*(ii + jj*m)+1];
  ierr = VecRestoreArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxCreateENMap2d_MPI(SpecFECtx c)
{
  PetscErrorCode ierr;
  PetscInt ni0,nj0,i,j,ei,ej,ecnt,*emap,nid;
  PetscInt si,si_g,sj,sj_g,nx_local;
  
  
  ierr = PetscMalloc(sizeof(PetscInt)*c->ne*c->npe,&c->element);CHKERRQ(ierr);
  ierr = PetscMemzero(c->element,sizeof(PetscInt)*c->ne*c->npe);CHKERRQ(ierr);
  
  ierr = SpecFECtxGetCornerBasis_MPI(c,&si,&si_g,&sj,&sj_g);CHKERRQ(ierr);
  ierr = DMDAGetGhostCorners(c->dm,NULL,NULL,NULL,&nx_local,NULL,NULL);CHKERRQ(ierr);
  /*printf("rank %d : %d %d x %d %d\n",c->rank,si,si_g,sj,sj_g);*/
  
  ecnt = 0;
  for (ej=0; ej<c->my; ej++) {
    nj0 = sj + ej*(c->npe_1d-1);
    
    for (ei=0; ei<c->mx; ei++) {
      ni0 = si + ei*(c->npe_1d-1);
      
      emap = &c->element[c->npe*ecnt];
      
      for (j=0; j<c->npe_1d; j++) {
        for (i=0; i<c->npe_1d; i++) {
          
          nid = (ni0 + i) + (nj0 + j) * nx_local;
          emap[i+j*c->npe_1d] = nid;
          //if (c->rank == 0) {
          //  printf("e %d : %d [max %d]\n",ecnt,nid,c->ne*c->npe);
          //}
        }
      }
      
      ecnt++;
    }
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxCreateMeshCoords2d_MPI(SpecFECtx c)
{
  PetscErrorCode ierr;
  Vec coor,gcoor;
  DM cdm;
  PetscInt ei,ej,i,j,ni0,nj0,si,si_g,sj,sj_g,gi,gj,m,n;
  PetscReal dx,dy,x0,y0;
  PetscReal *LA_coor;
  
  ierr = DMDASetUniformCoordinates(c->dm,0.0,1.0,0.0,1.0,0,0);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(c->dm,&cdm);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);

  ierr = DMGetCoordinates(c->dm,&gcoor);CHKERRQ(ierr);
  ierr = VecZeroEntries(gcoor);CHKERRQ(ierr);
  
  ierr = SpecFECtxGetCornerBasis_MPI(c,&si,&si_g,&sj,&sj_g);CHKERRQ(ierr);
  ierr = DMDAGetGhostCorners(c->dm,&gi,&gj,NULL,&m,&n,NULL);CHKERRQ(ierr);

  dx = 1.0/((PetscReal)c->mx_g);
  dy = 1.0/((PetscReal)c->my_g);
  ierr = VecGetArray(coor,&LA_coor);CHKERRQ(ierr);
  for (ej=0; ej<c->my; ej++) {
    
    for (ei=0; ei<c->mx; ei++) {
      if ( si >= m*n) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Out of range-si");
      if ( sj >= m*n) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Out of range-sj");

      x0 = LA_coor[2*(si + sj*m)+0] + ei*dx;
      y0 = LA_coor[2*(si + sj*m)+1] + ej*dy;

      ni0 = si + ei*(c->npe_1d-1);
      nj0 = sj + ej*(c->npe_1d-1);
      
      //printf("rank %d : (%d,%d) -> %d %d  %+1.4e %+1.4e\n",c->rank,ei,ej,ni0,nj0,x0,y0);
      
      for (j=0; j<c->npe_1d; j++) {
        for (i=0; i<c->npe_1d; i++) {
          if ( (ni0+i)+(nj0+j)*m >= m*n) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Local index out of range");

          LA_coor[2*((ni0+i) + (nj0+j)*m)+0] = x0 + 0.5*(c->xi1d[i]+1.0)*dx;
          LA_coor[2*((ni0+i) + (nj0+j)*m)+1] = y0 + 0.5*(c->xi1d[j]+1.0)*dy;
          
          ierr = VecSetValueLocal(gcoor,2*((ni0+i) + (nj0+j)*m)+0,x0 + 0.5*(c->xi1d[i]+1.0)*dx,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecSetValueLocal(gcoor,2*((ni0+i) + (nj0+j)*m)+1,y0 + 0.5*(c->xi1d[j]+1.0)*dy,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = VecRestoreArray(coor,&LA_coor);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(gcoor);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(gcoor);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxCreateMesh_MPI(SpecFECtx c,PetscInt dim,PetscInt mx,PetscInt my,PetscInt mz,PetscInt basisorder,PetscInt ndofs)
{
  PetscErrorCode ierr;
  PetscInt stencil_width,i,j;
  DM dm_ref;
  PetscInt ranks[3];
  PetscInt r,*lx,*ly;
  const PetscInt *lx_ref,*ly_ref;
  DMDALocalInfo info;
  
  c->dim = dim;
  c->mx_g = mx;
  c->my_g = my;
  c->mz_g = mz;
  c->basisorder = basisorder;
  c->dofs = ndofs;
  
  c->nx_g = basisorder*mx + 1;
  c->ny_g = basisorder*my + 1;
  c->nz_g = basisorder*mz + 1;
  
  ierr = CreateGLLCoordsWeights(basisorder,&c->npe_1d,&c->xi1d,&c->w1d);CHKERRQ(ierr);
  
  stencil_width = 1;
  switch (dim) {
    case 2:
    c->npe = c->npe_1d * c->npe_1d;
    c->ne_g = mx * my;

    ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
                        c->mx_g,c->my_g,PETSC_DECIDE,PETSC_DECIDE,1,0,NULL,NULL,&dm_ref);CHKERRQ(ierr);
    ierr = DMSetUp(dm_ref);CHKERRQ(ierr);
    /*ierr = DMView(dm_ref,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */

    ierr = DMDAGetInfo(dm_ref,NULL,NULL,NULL,NULL,&ranks[0],&ranks[1],NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetOwnershipRanges(dm_ref,&lx_ref,&ly_ref,NULL);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dm_ref,&info);CHKERRQ(ierr);
    
    c->mx = info.xm;
    c->my = info.ym;
    c->ne = c->mx * c->my;
    
    ierr = PetscMalloc1(ranks[0],&lx);CHKERRQ(ierr);
    ierr = PetscMalloc1(ranks[1],&ly);CHKERRQ(ierr);
    for (r=0; r<ranks[0]; r++) {
      lx[r] = lx_ref[r] * (c->npe_1d - 1);
    }
    lx[ranks[0]-1]++;

    /*for (r=0; r<ranks[0]; r++)  PetscPrintf(PETSC_COMM_WORLD,"npoints-i[%D] %D \n",r,lx[r]);*/
    
    for (r=0; r<ranks[1]; r++) {
      ly[r] = ly_ref[r] * (c->npe_1d - 1);
    }
    ly[ranks[1]-1]++;
    
    /*for (r=0; r<ranks[1]; r++)  PetscPrintf(PETSC_COMM_WORLD,"npoints-j[%D] %D \n",r,ly[r]);*/
    
    ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
                        c->nx_g,c->ny_g,ranks[0],ranks[1],ndofs,stencil_width,lx,ly,&c->dm);CHKERRQ(ierr);
    ierr = DMSetUp(c->dm);CHKERRQ(ierr);
    /*ierr = DMView(c->dm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */

    ierr = SpecFECtxCreateENMap2d_MPI(c);CHKERRQ(ierr);
    ierr = SpecFECtxCreateMeshCoords2d_MPI(c);CHKERRQ(ierr);
    
    /* tensor product for weights */
    ierr = PetscMalloc(sizeof(PetscReal)*c->npe,&c->w);CHKERRQ(ierr);
    for (j=0; j<c->npe_1d; j++) {
      for (i=0; i<c->npe_1d; i++) {
        c->w[i+j*c->npe_1d] = c->w1d[i] * c->w1d[j];
      }
    }
    
    ierr = TabulateBasisDerivativesTensorProduct2d(basisorder,&c->dN_dxi,&c->dN_deta);CHKERRQ(ierr);
    ierr = TabulateBasisDerivativesTensorProduct2d(basisorder,&c->dN_dx,&c->dN_dy);CHKERRQ(ierr);
    
    ierr = PetscFree(lx);CHKERRQ(ierr);
    ierr = PetscFree(ly);CHKERRQ(ierr);
    ierr = DMDestroy(&dm_ref);CHKERRQ(ierr);
    break;
  }
  
  c->nqp = c->npe;
  
  ierr = PetscMalloc(sizeof(QPntIsotropicElastic)*c->ne,&c->cell_data);CHKERRQ(ierr);
  ierr = PetscMemzero(c->cell_data,sizeof(QPntIsotropicElastic)*c->ne);CHKERRQ(ierr);
  
  ierr = PetscMalloc(sizeof(PetscReal)*c->npe*c->dim,&c->elbuf_coor);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*c->npe*c->dofs,&c->elbuf_field);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*c->npe*c->dofs,&c->elbuf_field2);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*c->npe*c->dofs,&c->elbuf_field3);CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscInt)*c->npe*c->dofs,&c->elbuf_dofs);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxCreateMesh(SpecFECtx c,PetscInt dim,PetscInt mx,PetscInt my,PetscInt mz,PetscInt basisorder,PetscInt ndofs)
{
  PetscMPIInt size;
  PetscErrorCode ierr;
  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size == 1) {
    ierr = SpecFECtxCreateMesh_SEQ(c,dim,mx,my,mz,basisorder,ndofs);CHKERRQ(ierr);
  } else {
    ierr = SpecFECtxCreateMesh_MPI(c,dim,mx,my,mz,basisorder,ndofs);CHKERRQ(ierr);
  }
  ierr = DMDASetFieldName(c->dm,0,"_x");CHKERRQ(ierr);
  ierr = DMDASetFieldName(c->dm,1,"_y");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxSetConstantMaterialProperties(SpecFECtx c,PetscReal lambda,PetscReal mu,PetscReal rho)
{
  PetscInt q;
  
  for (q=0; q<c->ne; q++) {
    c->cell_data[q].lambda = lambda;
    c->cell_data[q].mu     = mu;
    c->cell_data[q].rho    = rho;
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxSetConstantMaterialProperties_Velocity(SpecFECtx c,PetscReal Vp,PetscReal Vs,PetscReal rho)
{
  PetscErrorCode ierr;
  PetscReal mu,lambda;
  
  mu = Vs * Vs * rho;
  lambda = Vp * Vp * rho - 2.0 * mu;
  ierr = SpecFECtxSetConstantMaterialProperties(c,lambda,mu,rho);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"  [material]     Vp = %1.8e\n",Vp);
  PetscPrintf(PETSC_COMM_WORLD,"  [material]     Vs = %1.8e\n",Vs);
  PetscPrintf(PETSC_COMM_WORLD,"  [material] lambda = %1.8e\n",lambda);
  PetscPrintf(PETSC_COMM_WORLD,"  [material]     mu = %1.8e\n",mu);
  PetscPrintf(PETSC_COMM_WORLD,"  [material]    rho = %1.8e\n",rho);
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxSetPerturbedMaterialProperties_Velocity(SpecFECtx c,PetscReal Vp0,PetscReal delta_Vp,PetscReal Vs0,PetscReal delta_Vs,PetscReal rho0,PetscReal delta_rho)
{
  PetscErrorCode ierr;
  Vec Vp,Vs,rho;
  PetscRandom r;
  const PetscReal *LA_Vp,*LA_Vs,*LA_rho;
  PetscInt e;
  
  ierr = VecCreate(PETSC_COMM_WORLD,&Vp);CHKERRQ(ierr);
  ierr = VecSetSizes(Vp,c->ne,c->ne_g);CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vp);CHKERRQ(ierr);
  ierr = VecDuplicate(Vp,&Vs);CHKERRQ(ierr);
  ierr = VecDuplicate(Vp,&rho);CHKERRQ(ierr);

  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
  ierr = PetscRandomSetType(r,PETSCRAND48);CHKERRQ(ierr);

  ierr = PetscRandomSetInterval(r,Vp0-delta_Vp,Vp0+delta_Vp);CHKERRQ(ierr);
  ierr = PetscRandomSetSeed(r,1);CHKERRQ(ierr);
  ierr = PetscRandomSeed(r);CHKERRQ(ierr);
  ierr = VecSetRandom(Vp,r);CHKERRQ(ierr);

  ierr = PetscRandomSetInterval(r,Vs0-delta_Vs,Vs0+delta_Vs);CHKERRQ(ierr);
  ierr = PetscRandomSetSeed(r,2);CHKERRQ(ierr);
  ierr = PetscRandomSeed(r);CHKERRQ(ierr);
  ierr = VecSetRandom(Vs,r);CHKERRQ(ierr);
  
  ierr = PetscRandomSetInterval(r,rho0-delta_rho,rho0+delta_rho);CHKERRQ(ierr);
  ierr = PetscRandomSetSeed(r,3);CHKERRQ(ierr);
  ierr = PetscRandomSeed(r);CHKERRQ(ierr);
  ierr = VecSetRandom(rho,r);CHKERRQ(ierr);

  ierr = PetscRandomDestroy(&r);CHKERRQ(ierr);

  ierr = VecGetArrayRead(Vp,&LA_Vp);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Vs,&LA_Vs);CHKERRQ(ierr);
  ierr = VecGetArrayRead(rho,&LA_rho);CHKERRQ(ierr);
  for (e=0; e<c->ne; e++) {
    PetscReal mu,lambda;

    mu     = LA_Vs[e] * LA_Vs[e] * LA_rho[e];
    lambda = LA_Vp[e] * LA_Vp[e] * LA_rho[e] - 2.0 * mu;
    
    c->cell_data[e].lambda = lambda;
    c->cell_data[e].mu     = mu;
    c->cell_data[e].rho    = LA_rho[e];
  }
  ierr = VecRestoreArrayRead(rho,&LA_rho);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Vs,&LA_Vs);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Vp,&LA_Vp);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"  [material]     Vp0 = %1.8e : delta = %+1.8e\n",Vp0,delta_Vp);
  PetscPrintf(PETSC_COMM_WORLD,"  [material]     Vs0 = %1.8e : delta = %+1.8e\n",Vs0,delta_Vs);
  PetscPrintf(PETSC_COMM_WORLD,"  [material]    rho0 = %1.8e : delta = %+1.8e\n",rho0,delta_rho);
  
  ierr = VecDestroy(&rho);CHKERRQ(ierr);
  ierr = VecDestroy(&Vs);CHKERRQ(ierr);
  ierr = VecDestroy(&Vp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void ElementEvaluateGeometry_CellWiseConstant2d(PetscInt npe,PetscReal el_coords[],
                                                PetscInt nbasis,PetscReal *detJ)
{
  PetscReal J00,J11;
  PetscReal dx,dy;
  
  dx = el_coords[2*(nbasis-1)+0] - el_coords[2*0+0];
  dy = el_coords[2*(npe-1)+1]    - el_coords[2*0+1];
  
  J00 = 0.5 * dx;
  J11 = 0.5 * dy;
  
  *detJ = J00*J11;
}

void ElementEvaluateDerivatives_CellWiseConstant2d(PetscInt nqp,PetscInt npe,PetscReal el_coords[],
                                                   PetscInt nbasis,PetscReal **dN_dxi,PetscReal **dN_deta,
                                                   PetscReal **dN_dx,PetscReal **dN_dy)
{
  PetscInt k,q;
  PetscReal J00,J11,iJ00,iJ11;
  PetscReal dx,dy;
  
  dx = el_coords[2*(nbasis-1)+0] - el_coords[2*0+0];
  dy = el_coords[2*(npe-1)+1]    - el_coords[2*0+1];
  
  J00 = 0.5 * dx;
  J11 = 0.5 * dy;
  
  for (q=0; q<nqp; q++) {
    
    iJ00 = 1.0/J00;
    iJ11 = 1.0/J11;
    
    /* shape function derivatives */
    for (k=0; k<npe; k++) {
      dN_dx[q][k] = iJ00 * dN_dxi[q][k];
      dN_dy[q][k] = iJ11 * dN_deta[q][k];
    }
  }
}

/*
 Assemble rhs
 L(u) = - \int B^T D B u dV
 */
PetscErrorCode AssembleLinearForm_ElastoDynamics2d(SpecFECtx c,Vec u,Vec F)
{
  PetscErrorCode ierr;
  PetscInt  e,nqp,q,i,nbasis,ndof;
  PetscReal e_vec[3],sigma_vec[3];
  PetscInt  *element,*elnidx,*eldofs;
  PetscReal *fe,*ux,*uy,*elcoords,detJ,*field;
  Vec       coor,ul,fl;
  const PetscReal *LA_coor,*LA_u;
  QPntIsotropicElastic *celldata;
  
  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  
  eldofs   = c->elbuf_dofs;
  elcoords = c->elbuf_coor;
  nbasis   = c->npe;
  nqp      = c->nqp;
  ndof     = c->dofs;
  fe       = c->elbuf_field;
  element  = c->element;
  field    = c->elbuf_field2;
  
  ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(c->dm,&ul);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(c->dm,u,INSERT_VALUES,ul);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,u,INSERT_VALUES,ul);CHKERRQ(ierr);
  ierr = VecGetArrayRead(ul,&LA_u);CHKERRQ(ierr);

  ierr = DMGetLocalVector(c->dm,&fl);CHKERRQ(ierr);
  ierr = VecZeroEntries(fl);CHKERRQ(ierr);

  ux = &field[0];
  uy = &field[nbasis];
  
  for (e=0; e<c->ne; e++) {
    /* get element -> node map */
    elnidx = &element[nbasis*e];
    
    
    /* generate dofs */
    for (i=0; i<nbasis; i++) {
      eldofs[2*i  ] = 2*elnidx[i];
      eldofs[2*i+1] = 2*elnidx[i]+1;
    }
    
    /* get element coordinates */
    for (i=0; i<nbasis; i++) {
      PetscInt nidx = elnidx[i];
      elcoords[2*i  ] = LA_coor[2*nidx  ];
      elcoords[2*i+1] = LA_coor[2*nidx+1];
    }
    
    /* get element displacements */
    for (i=0; i<nbasis; i++) {
      PetscInt nidx = elnidx[i];
      ux[i] = LA_u[2*nidx  ];
      uy[i] = LA_u[2*nidx+1];
    }
    
    /* compute derivatives */
    ElementEvaluateGeometry_CellWiseConstant2d(nbasis,elcoords,c->npe_1d,&detJ);
    ElementEvaluateDerivatives_CellWiseConstant2d(nqp,nbasis,elcoords,
                                                  c->npe_1d,c->dN_dxi,c->dN_deta,
                                                  c->dN_dx,c->dN_dy);
    
    ierr = PetscMemzero(fe,sizeof(PetscReal)*nbasis*ndof);CHKERRQ(ierr);
    
    /* get access to element->quadrature points */
    celldata = &c->cell_data[e];
    
    for (q=0; q<c->nqp; q++) {
      PetscReal            fac;
      PetscReal            c11,c12,c21,c22,c33,lambda_qp,mu_qp;
      PetscReal            *dNidx,*dNidy;
      
      
      dNidx = c->dN_dx[q];
      dNidy = c->dN_dy[q];
      
      /* compute strain @ quadrature point */
      /*
       e = Bu = [ d/dx  0    ][ u v ]^T
       [ 0     d/dy ]
       [ d/dy  d/dx ]
       */
      e_vec[0] = e_vec[1] = e_vec[2] = 0.0;
      for (i=0; i<nbasis; i++) {
        e_vec[0] += dNidx[i] * ux[i];
        e_vec[1] += dNidy[i] * uy[i];
        e_vec[2] += (dNidx[i] * uy[i] + dNidy[i] * ux[i]);
      }
      
      /* evaluate constitutive model */
      lambda_qp  = celldata->lambda;
      mu_qp      = celldata->mu;
      
      /*
       coeff = E_qp * (1.0 + nu_qp)/(1.0 - 2.0*nu_qp);
       c11 = coeff*(1.0 - nu_qp);
       c12 = coeff*(nu_qp);
       c21 = coeff*(nu_qp);
       c22 = coeff*(1.0 - nu_qp);
       c33 = coeff*(0.5 * (1.0 - 2.0 * nu_qp));
       */
      c11 = 2.0*mu_qp + lambda_qp;
      c12 = lambda_qp;
      c21 = lambda_qp;
      c22 = 2.0*mu_qp + lambda_qp;
      c33 = mu_qp;
      
      /* compute stress @ quadrature point */
      sigma_vec[TENS2D_XX] = c11 * e_vec[0] + c12 * e_vec[1];
      sigma_vec[TENS2D_YY] = c21 * e_vec[0] + c22 * e_vec[1];
      sigma_vec[TENS2D_XY] = c33 * e_vec[2];
      //printf("s = %1.4e %1.4e %1.4e \n",sigma_vec[0],sigma_vec[1],sigma_vec[2]);
      /*
       a(u,v) = B^T s
       = [ d/dx  0    d/dy ][ sxx syy sxy ]^T
       [ 0     d/dy d/dx ]
       */
      
      fac = detJ * c->w[q];
      
      for (i=0; i<nbasis; i++) {
        fe[2*i  ] += -fac * (dNidx[i] * sigma_vec[TENS2D_XX] + dNidy[i] * sigma_vec[TENS2D_XY]);
        fe[2*i+1] += -fac * (dNidy[i] * sigma_vec[TENS2D_YY] + dNidx[i] * sigma_vec[TENS2D_XY]);
      }
      
    }
    //ierr = VecSetValuesLocal(F,nbasis*ndof,eldofs,fe,ADD_VALUES);CHKERRQ(ierr);
    ierr = VecSetValues(fl,nbasis*ndof,eldofs,fe,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(fl);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(fl);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(c->dm,fl,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(c->dm,fl,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(ul,&LA_u);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&ul);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&fl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SpecFECtxGetDRCellData(SpecFECtx c,PetscInt e_index,DRVar **data)
{
  *data = &c->dr_qp_data[c->nqp * e_index];
  PetscFunctionReturn(0);
}

PetscErrorCode EvaluateVelocityAtPoint(SpecFECtx c,const PetscReal LA_v[],PetscReal xr[],PetscReal vr[])
{
  PetscReal      gmin[3],gmax[3],dx,dy;
  PetscInt       k,ei,ej,eid,*element,*elbasis;
  PetscReal      N[400];
  PetscErrorCode ierr;

  
  if (c->size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Needs updating to support MPI");
  
  /* get containing element */
  ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
  dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
  ei = (xr[0] - gmin[0])/dx; /* todo - needs to be sub-domain gmin */
  
  dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);
  ej = (xr[1] - gmin[1])/dy;
  
  eid = ei + ej * c->mx;
  
  /* get element -> node map */
  element = c->element;
  elbasis = &element[c->npe*eid];
  
  {
    PetscInt  nbasis,i,j;
    PetscReal **N_s1,**N_s2,xi,eta,x0,y0;
    
    x0 = gmin[0] + ei*dx; /* todo - needs to be sub-domain gmin */
    y0 = gmin[1] + ej*dy;
    
    // (xi - (-1))/2 = (x - x0)/dx
    xi = 2.0*(xr[0] - x0)/dx - 1.0;
    eta = 2.0*(xr[1] - y0)/dy - 1.0;
    
    /* compute basis */
    ierr = TabulateBasis1d_CLEGENDRE(1,&xi,c->basisorder,&nbasis,&N_s1);CHKERRQ(ierr);
    ierr = TabulateBasis1d_CLEGENDRE(1,&eta,c->basisorder,&nbasis,&N_s2);CHKERRQ(ierr);
    
    k = 0;
    for (j=0; j<c->npe_1d; j++) {
      for (i=0; i<c->npe_1d; i++) {
        N[k] = N_s1[0][i] * N_s2[0][j];
        k++;
      }
    }
    
    ierr = PetscFree(N_s1[0]);CHKERRQ(ierr);
    ierr = PetscFree(N_s1);CHKERRQ(ierr);
    ierr = PetscFree(N_s2[0]);CHKERRQ(ierr);
    ierr = PetscFree(N_s2);CHKERRQ(ierr);
  }
  
  vr[0] = vr[1] = 0.0;
  for (k=0; k<c->npe; k++) {
    PetscInt nid = elbasis[k];
    
    vr[0] += N[k] * LA_v[2*nid+0];
    vr[1] += N[k] * LA_v[2*nid+1];
  }
  
  PetscFunctionReturn(0);
}


PetscErrorCode PointLocation_v2(SpecFECtx c,const PetscReal xr[],PetscInt *_eid,PetscReal **N1,PetscReal **N2)
{
  static PetscBool beenhere = PETSC_FALSE;
  static PetscReal      gmin[3],gmax[3];
  PetscReal      dx,dy;
  PetscInt       ei,ej,eid;
  PetscInt       nbasis;
  PetscReal      **N_s1,**N_s2,xi,eta,x0,y0;
  PetscErrorCode ierr;
  
  
  if (c->size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Needs updating to support MPI");
  
  /* get containing element */
  if (!beenhere) {
    ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
    beenhere = PETSC_TRUE;
  }
  dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
  ei = (xr[0] - gmin[0])/dx; /* todo - needs to be sub-domain gmin */
  
  dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);
  ej = (xr[1] - gmin[1])/dy;
  
  eid = ei + ej * c->mx;
  
  x0 = gmin[0] + ei*dx; /* todo - needs to be sub-domain gmin */
  y0 = gmin[1] + ej*dy;
  
  // (xi - (-1))/2 = (x - x0)/dx
  xi = 2.0*(xr[0] - x0)/dx - 1.0;
  eta = 2.0*(xr[1] - y0)/dy - 1.0;
  
  /* compute basis */
  ierr = TabulateBasis1d_CLEGENDRE(1,&xi,c->basisorder,&nbasis,&N_s1);CHKERRQ(ierr);
  ierr = TabulateBasis1d_CLEGENDRE(1,&eta,c->basisorder,&nbasis,&N_s2);CHKERRQ(ierr);
  
  *_eid = eid;
  *N1 = N_s1[0];
  *N2 = N_s2[0];
  
  //ierr = PetscFree(N_s1[0]);CHKERRQ(ierr);
  ierr = PetscFree(N_s1);CHKERRQ(ierr);
  //ierr = PetscFree(N_s2[0]);CHKERRQ(ierr);
  ierr = PetscFree(N_s2);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


PetscErrorCode FaultSDFInit_v2(SpecFECtx c)
{
  PetscErrorCode  ierr;
  PetscInt        e,i,q,nbasis,nqp,ndof;
  Vec             coor;
  const PetscReal *LA_coor;
  PetscInt        *element,*elnidx,*eldofs;
  PetscReal       *elcoords;
  DRVar           *dr_celldata;
  PetscReal       factor;
  PetscInt        factor_i;
  
  eldofs   = c->elbuf_dofs;
  elcoords = c->elbuf_coor;
  nbasis   = c->npe;
  nqp      = c->nqp;
  ndof     = c->dofs;
  element  = c->element;
  

  factor = ((PetscReal)(c->ne)) * 0.1;
  factor_i = (PetscInt)factor;
  if (factor_i == 0) { factor_i = 1; }
  
  ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  for (e=0; e<c->ne; e++) {
    /* get element -> node map */
    elnidx = &element[nbasis*e];
    
    
    /* generate dofs */
    for (i=0; i<nbasis; i++) {
      eldofs[2*i  ] = 2*elnidx[i];
      eldofs[2*i+1] = 2*elnidx[i]+1;
    }
    
    /* get element coordinates */
    for (i=0; i<nbasis; i++) {
      PetscInt nidx = elnidx[i];
      elcoords[2*i  ] = LA_coor[2*nidx  ];
      elcoords[2*i+1] = LA_coor[2*nidx+1];
    }
    
    ierr = SpecFECtxGetDRCellData(c,e,&dr_celldata);CHKERRQ(ierr);
    
    for (q=0; q<c->nqp; q++) {
      PetscReal coor_qp[2];
      PetscBool modify_stress_state;
      
      coor_qp[0] = elcoords[2*q  ];
      coor_qp[1] = elcoords[2*q+1];
      
      modify_stress_state = PETSC_FALSE;
      dr_celldata[q].eid[0] = -1;
      dr_celldata[q].eid[1] = -1;
      
      ierr = FaultSDFQuery(coor_qp,c->delta,NULL,&modify_stress_state);CHKERRQ(ierr);
      
      if (modify_stress_state) {
        PetscReal x_plus[2],x_minus[2];
        
        //printf("[e %d , q %d] x_qp %+1.4e , %+1.4e\n",e,q,coor_qp[0],coor_qp[1]);
        
        ierr = FaultSDFGetPlusMinusCoor(coor_qp,c->delta,NULL,x_plus,x_minus);CHKERRQ(ierr);
        
        ierr = PointLocation_v2(c,(const PetscReal*)x_plus, &dr_celldata[q].eid[0],&dr_celldata[q].N1_plus,&dr_celldata[q].N2_plus);CHKERRQ(ierr);
        ierr = PointLocation_v2(c,(const PetscReal*)x_minus,&dr_celldata[q].eid[1],&dr_celldata[q].N1_minus,&dr_celldata[q].N2_minus);CHKERRQ(ierr);
        //printf("  [e %d,q %d] x_qp -> x+ %+1.4e , %+1.4e [eid %d]\n",e,q,x_plus[0],x_plus[1],dr_celldata[q].eid[0]);
        //printf("  [e %d,q %d] x_qp -> x- %+1.4e , %+1.4e [eid %d]\n",e,q,x_minus[0],x_minus[1],dr_celldata[q].eid[1]);
      }
      
    }
    if (e%factor_i == 0) {
      printf("[Fault point location] Done element %d of %d\n",e,c->ne);
    }
  }
  printf("[Fault point location-v2] Finished\n");
  
  ierr = VecRestoreArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/* Use tabulated basis at delta(+,-) to interpolate velocity */
PetscErrorCode FaultSDFTabulateInterpolation_v2(SpecFECtx c,const PetscReal LA_v[],DRVar *dr_celldata_q,
                                                PetscReal v_plus[],PetscReal v_minus[])
{
  PetscInt  eid_delta,*element_delta;
  PetscReal vx_d_e,vy_d_e,Ni;
  PetscInt i,ii,jj,nbasis,*element,nidx;
  
  nbasis   = c->npe;
  element  = c->element;
  
  if (!dr_celldata_q->N1_plus || !dr_celldata_q->N2_plus) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"[+] N1,N2 not allocated");
  eid_delta = dr_celldata_q->eid[0];
  element_delta = &element[nbasis*eid_delta];
  
  v_plus[0] = v_plus[1] = 0.0;
  i = 0;
  for (jj=0; jj<c->npe_1d; jj++) {
    for (ii=0; ii<c->npe_1d; ii++) {
      nidx = element_delta[i];
      
      Ni = dr_celldata_q->N1_plus[ii] * dr_celldata_q->N2_plus[jj];
      
      vx_d_e = LA_v[2*nidx  ];
      vy_d_e = LA_v[2*nidx+1];
      v_plus[0] += Ni * vx_d_e;
      v_plus[1] += Ni * vy_d_e;
      
      i++;
    }
  }

  if (!dr_celldata_q->N1_minus || !dr_celldata_q->N2_minus) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"[-] N1,N2 not allocated");
  eid_delta = dr_celldata_q->eid[1];
  element_delta = &element[nbasis*eid_delta];
  
  v_minus[0] = v_minus[1] = 0.0;
  i = 0;
  for (jj=0; jj<c->npe_1d; jj++) {
    for (ii=0; ii<c->npe_1d; ii++) {
      nidx = element_delta[i];
      
      Ni = dr_celldata_q->N1_minus[ii] * dr_celldata_q->N2_minus[jj];
      
      vx_d_e = LA_v[2*nidx  ];
      vy_d_e = LA_v[2*nidx+1];
      v_minus[0] += Ni * vx_d_e;
      v_minus[1] += Ni * vy_d_e;
      
      i++;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VoigtTensorContract_ai_Tij_bj(PetscReal a[],PetscReal t[2][2],PetscReal b[],PetscReal *r)
{
  PetscReal s=0;
  PetscInt i,j;
  
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      s += a[i] * t[i][j] * b[j];
    }
  }
  *r = s;
  PetscFunctionReturn(0);
}

PetscErrorCode VoigtTensorConvert(PetscReal T[],PetscReal t[2][2])
{
  t[0][0] = T[TENS2D_XX];
  t[0][1] = T[TENS2D_XY];
  t[1][0] = T[TENS2D_XY];
  t[1][1] = T[TENS2D_YY];
  PetscFunctionReturn(0);
}

PetscErrorCode TensorConvertToVoigt(PetscReal t[2][2],PetscReal T[])
{
  T[TENS2D_XX] = t[0][0];
  T[TENS2D_XY] = t[0][1];
  T[TENS2D_YY] = t[1][1];
  if (fabs(t[1][0] - t[0][1]) > 1.0e-12) {
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Cannot convert non-symmetric tensor into Voigt format");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TensorZeroEntries(PetscReal t[2][2])
{
  t[0][0] = 0;
  t[0][1] = 0;
  t[1][0] = 0;
  t[1][1] = 0;
  PetscFunctionReturn(0);
}

PetscErrorCode TensorScale(PetscReal y[2][2],PetscReal a)
{
  y[0][0] = a*y[0][0];
  y[0][1] = a*y[0][1];
  y[1][0] = a*y[1][0];
  y[1][1] = a*y[1][1];
  PetscFunctionReturn(0);
}

PetscErrorCode TensorAXPY(PetscReal y[2][2],PetscReal a,PetscReal x[2][2])
{
  y[0][0] += a*x[0][0];
  y[0][1] += a*x[0][1];
  y[1][0] += a*x[1][0];
  y[1][1] += a*x[1][1];
  PetscFunctionReturn(0);
}


PetscErrorCode VectorContract_ai_bj(PetscReal a[],PetscReal b[],PetscReal t[2][2])
{
  PetscInt i,j;
  
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      t[i][j] = a[i] * b[j];
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VectorContractAdd_ai_bj(PetscReal a[],PetscReal b[],PetscReal t[2][2])
{
  PetscInt i,j;
  
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      t[i][j] += a[i] * b[j];
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode VectorContractAbsAdd_ai_bj(PetscReal a[],PetscReal b[],PetscReal t[2][2])
{
  PetscInt i,j;
  
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      t[i][j] += fabs(a[i] * b[j]);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TensorRtAR(PetscReal R[2][2],PetscReal T[2][2],PetscReal Tr[2][2])
{
  PetscInt i,j,k,l;
  
  // Tr[i][j] = Rt[i][k]T[k][l]R[l][j] = Rt[k][i]T[k][l]R[l][j]
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      Tr[i][j] = 0.0;
      for (k=0; k<2; k++) {
        for (l=0; l<2; l++) {
          Tr[i][j] += R[k][i] * T[k][l] * R[l][j];
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TensorTransform(PetscReal e1[],PetscReal e2[],PetscReal T[2][2],PetscReal Tr[2][2])
{
  PetscReal R[2][2];
  
  R[0][0] = e1[0]; R[0][1] = e2[0];
  R[1][0] = e1[1]; R[1][1] = e2[1];
  TensorRtAR(R,T,Tr);
  PetscFunctionReturn(0);
}

PetscErrorCode TensorInverseTransform(PetscReal e1[],PetscReal e2[],PetscReal T[2][2],PetscReal Tr[2][2])
{
  PetscReal R[2][2],iR[2][2],det;
  
  R[0][0] = e1[0]; R[0][1] = e2[0];
  R[1][0] = e1[1]; R[1][1] = e2[1];
  det = R[0][0] * R[1][1] - R[0][1] * R[1][0];
  iR[0][0] =  R[1][1]/det;
  iR[0][1] = -R[0][1]/det;
  iR[1][0] = -R[1][0]/det;
  iR[1][1] =  R[0][0]/det;
  TensorRtAR(iR,T,Tr);
  PetscFunctionReturn(0);
}

PetscErrorCode AssembleLinearForm_ElastoDynamics_StressGlut2d_tpv(SpecFECtx c,Vec u,Vec v,PetscReal dt,PetscReal time,PetscReal gamma,Vec F)
{
  PetscErrorCode ierr;
  PetscInt  e,nqp,q,i,nbasis,ndof;
  PetscReal e_vec[3],edot_vec[3],sigma_vec[3],sigma_trial[3],gradu[4],gradv[4],gradv_q[9*9][4];
  PetscInt  *element,*elnidx,*eldofs;
  PetscReal *fe,*ux,*uy,*vx,*vy,*elcoords,detJ,*fieldU,*fieldV;
  Vec       coor,ul,vl,fl;
  const PetscReal *LA_coor,*LA_u,*LA_v;
  QPntIsotropicElastic *celldata;
  DRVar                *dr_celldata;

  // 0.55 is a good fit to matching sem2pack's displacement field...near the point (0,500)
  PetscReal sigma_n_0 = 120.0 * 1.0e6 * 1.0;
  PetscReal sigma_t_0 = 70.0  * 1.0e6 * 1.0;
  PetscReal sigma_n_1  = 120.0 * 1.0e6 * 1.0;
  PetscReal sigma_t_1  = 81.6  * 1.0e6 * 1.0;
  static PetscBool beenhere = PETSC_FALSE;
  static PetscReal gmin[3],gmax[3];
  PetscReal dx,dy;
  if (!beenhere) {
    ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
    beenhere = PETSC_TRUE;
  }
  dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
  dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);

  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  
  eldofs   = c->elbuf_dofs;
  elcoords = c->elbuf_coor;
  nbasis   = c->npe;
  nqp      = c->nqp;
  ndof     = c->dofs;
  fe       = c->elbuf_field;
  element  = c->element;
  fieldU   = c->elbuf_field2;
  fieldV   = c->elbuf_field3;
  
  ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(c->dm,&ul);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(c->dm,u,INSERT_VALUES,ul);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,u,INSERT_VALUES,ul);CHKERRQ(ierr);
  ierr = VecGetArrayRead(ul,&LA_u);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(c->dm,&vl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(c->dm,v,INSERT_VALUES,vl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,v,INSERT_VALUES,vl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(vl,&LA_v);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(c->dm,&fl);CHKERRQ(ierr);
  ierr = VecZeroEntries(fl);CHKERRQ(ierr);
  
  ux = &fieldU[0];
  uy = &fieldU[nbasis];
  
  vx = &fieldV[0];
  vy = &fieldV[nbasis];
  
  
  for (e=0; e<c->ne; e++) {
    ierr = SpecFECtxGetDRCellData(c,e,&dr_celldata);CHKERRQ(ierr);
  }
  
  for (e=0; e<c->ne; e++) {
    PetscReal x_cell[] = {0,0};
    
    /* get element -> node map */
    elnidx = &element[nbasis*e];
    
    /* generate dofs */
    for (i=0; i<nbasis; i++) {
      eldofs[2*i  ] = 2*elnidx[i];
      eldofs[2*i+1] = 2*elnidx[i]+1;
    }
    
    /* get element coordinates */
    for (i=0; i<nbasis; i++) {
      PetscInt nidx = elnidx[i];
      elcoords[2*i  ] = LA_coor[2*nidx  ];
      elcoords[2*i+1] = LA_coor[2*nidx+1];
      x_cell[0] += elcoords[2*i  ];
      x_cell[1] += elcoords[2*i+1];
    }
    x_cell[0] = x_cell[0] / ((PetscReal)nbasis);
    x_cell[1] = x_cell[1] / ((PetscReal)nbasis);
    
    /* get element displacements & velocities */
    for (i=0; i<nbasis; i++) {
      PetscInt nidx = elnidx[i];
      ux[i] = LA_u[2*nidx  ];
      uy[i] = LA_u[2*nidx+1];
      
      vx[i] = LA_v[2*nidx  ];
      vy[i] = LA_v[2*nidx+1];
    }
    
    /* compute derivatives */
    ElementEvaluateGeometry_CellWiseConstant2d(nbasis,elcoords,c->npe_1d,&detJ);
    ElementEvaluateDerivatives_CellWiseConstant2d(nqp,nbasis,elcoords,
                                                  c->npe_1d,c->dN_dxi,c->dN_deta,
                                                  c->dN_dx,c->dN_dy);
    
    ierr = PetscMemzero(fe,sizeof(PetscReal)*nbasis*ndof);CHKERRQ(ierr);
    
    /* get access to element->quadrature points */
    celldata = &c->cell_data[e];
    
    ierr = SpecFECtxGetDRCellData(c,e,&dr_celldata);CHKERRQ(ierr);
    

    
    for (q=0; q<c->nqp; q++) {
      PetscReal *dNidx,*dNidy;
      
      dNidx = c->dN_dx[q];
      dNidy = c->dN_dy[q];
      
      gradv[0] = gradv[1] = gradv[2] = gradv[3] = 0.0;
      for (i=0; i<nbasis; i++) {
        gradv[0] += dNidx[i] * vx[i];
        gradv[1] += dNidy[i] * vx[i];
        gradv[2] += dNidx[i] * vy[i];
        gradv[3] += dNidy[i] * vy[i];
      }
      gradv_q[q][0] = gradv[0];
      gradv_q[q][1] = gradv[1];
      gradv_q[q][2] = gradv[2];
      gradv_q[q][3] = gradv[3];
    }

    
    for (q=0; q<c->nqp; q++) {
      PetscReal fac;
      PetscReal c11,c12,c21,c22,c33,lambda_qp,mu_qp;
      PetscReal *dNidx,*dNidy;
      PetscReal coor_qp[2];
      PetscBool inside_fault_region;
      PetscBool sliding_active;
      
      dNidx = c->dN_dx[q];
      dNidy = c->dN_dy[q];
      
      /* compute strain @ quadrature point */
      /*
       e = Bu = [ d/dx  0    ][ u v ]^T
       [ 0     d/dy ]
       [ d/dy  d/dx ]
       */
      e_vec[0] = e_vec[1] = e_vec[2] = 0.0;
      for (i=0; i<nbasis; i++) {
        e_vec[0] += dNidx[i] * ux[i];
        e_vec[1] += dNidy[i] * uy[i];
        e_vec[2] += (dNidx[i] * uy[i] + dNidy[i] * ux[i]);
      }

      edot_vec[0] = edot_vec[1] = edot_vec[2] = 0.0;
      for (i=0; i<nbasis; i++) {
        edot_vec[0] += dNidx[i] * vx[i];
        edot_vec[1] += dNidy[i] * vy[i];
        edot_vec[2] += (dNidx[i] * vy[i] + dNidy[i] * vx[i]);
      }

      gradu[0] = gradu[1] = gradu[2] = gradu[3] = 0.0;
      gradv[0] = gradv[1] = gradv[2] = gradv[3] = 0.0;
      for (i=0; i<nbasis; i++) {
        gradu[0] += dNidx[i] * ux[i];
        gradu[1] += dNidy[i] * ux[i];
        gradu[2] += dNidx[i] * uy[i];
        gradu[3] += dNidy[i] * uy[i];
        
        gradv[0] += dNidx[i] * vx[i];
        gradv[1] += dNidy[i] * vx[i];
        gradv[2] += dNidx[i] * vy[i];
        gradv[3] += dNidy[i] * vy[i];
      }

      coor_qp[0] = elcoords[2*q  ];
      coor_qp[1] = elcoords[2*q+1];

      /* evaluate constitutive model */
      lambda_qp = celldata->lambda;
      mu_qp     = celldata->mu;
      
      /*
       coeff = E_qp * (1.0 + nu_qp)/(1.0 - 2.0*nu_qp);
       c11 = coeff*(1.0 - nu_qp);
       c12 = coeff*(nu_qp);
       c21 = coeff*(nu_qp);
       c22 = coeff*(1.0 - nu_qp);
       c33 = coeff*(0.5 * (1.0 - 2.0 * nu_qp));
       */
      c11 = 2.0 * mu_qp + lambda_qp;
      c12 = lambda_qp;
      c21 = lambda_qp;
      c22 = 2.0 * mu_qp + lambda_qp;
      c33 = mu_qp;
      
      /* compute stress @ quadrature point */
      sigma_vec[TENS2D_XX] = c11 * e_vec[0] + c12 * e_vec[1];
      sigma_vec[TENS2D_YY] = c21 * e_vec[0] + c22 * e_vec[1];
      sigma_vec[TENS2D_XY] = c33 * e_vec[2];
            
      /*
       From
       Day and Ely "Effect of a Shallow Weak Zone on Fault Rupture: Numerical Simulation of Scale-Model Experiments",
       BSSA, 2002
       
       alpha = cp
       beta = cs
       volumetric terms; rho (cp^2 - 2 cs^2) gamma [div(v)]
       shear terms; rho cs^2 gamma [v_{i,j} + v_{j,i}]
       
       */
      /*
      //printf("lambda_qp * gamma %+1.4e : mu_qp * gamma %+1.4e\n",lambda_qp * gamma,mu_qp * gamma);
      sigma_vec[TENS2D_XX] += lambda_qp * gamma * edot_vec[0];
      sigma_vec[TENS2D_YY] += lambda_qp * gamma * edot_vec[1];
      sigma_vec[TENS2D_XY] += mu_qp * gamma * edot_vec[2];
      */
      
      {
        PetscReal factor = 1.0;
        
        //ierr = FaultSDFMollifer(coor_qp,2.0 * c->delta,NULL,&factor);CHKERRQ(ierr);
        //printf("%+1.4e %+1.4e %+1.4e\n",coor_qp[0],coor_qp[1],factor);

        c11 = factor * (2.0 * mu_qp + lambda_qp) * gamma;
        c12 = factor * (lambda_qp) * gamma;
        c21 = factor * (lambda_qp) * gamma;
        c22 = factor * (2.0 * mu_qp + lambda_qp) * gamma;
        c33 = factor * (mu_qp) * gamma;
      }

      /* compute stress @ quadrature point */
      sigma_vec[TENS2D_XX] += c11 * edot_vec[0] + c12 * edot_vec[1];
      sigma_vec[TENS2D_YY] += c21 * edot_vec[0] + c22 * edot_vec[1];
      sigma_vec[TENS2D_XY] += c33 * edot_vec[2];
      
      sigma_trial[TENS2D_XX] = sigma_vec[TENS2D_XX];
      sigma_trial[TENS2D_YY] = sigma_vec[TENS2D_YY];
      sigma_trial[TENS2D_XY] = sigma_vec[TENS2D_XY];

      
      
      inside_fault_region = PETSC_FALSE;
      
      ierr = FaultSDFQuery(coor_qp,c->delta,NULL,&inside_fault_region);CHKERRQ(ierr);
      if (fabs(x_cell[1]) > c->delta) { inside_fault_region = PETSC_FALSE; }
      
      inside_fault_region = PETSC_FALSE;
      if (fabs(x_cell[1]) < c->delta && fabs(x_cell[0]) < 15.0e3) { inside_fault_region = PETSC_TRUE; }

      /* NOTE - Not sure how to generalize the notion of an off-fault normal stress for non-planar geometries */
      /* NOTE - I'm not sure it is even well defined... */
      if (inside_fault_region) { /* add the initial stress state on fault */
        if (fabs(coor_qp[0]) < 1.5*1.0e3+1.0) {
          sigma_trial[TENS2D_XY] += sigma_t_1;
          sigma_trial[TENS2D_YY] += (-sigma_n_1); /* negative in compression */
        } else {
          sigma_trial[TENS2D_XY] += sigma_t_0;
          sigma_trial[TENS2D_YY] += (-sigma_n_0); /* negative in compression */
        }
      } else {
        sigma_trial[TENS2D_XY] += sigma_t_0;
        sigma_trial[TENS2D_YY] += (-sigma_n_0); /* negative in compression */
      }
      
      /* Make stress glut corrections here */
      if (inside_fault_region) {
        PetscReal x_plus[2],x_minus[2],v_plus[2],v_minus[2];
        PetscReal normal[2],tangent[2],Vplus,Vminus,slip,slip_k,slip_rate;
        PetscReal sigma_n,sigma_t,phi_p;
        PetscReal e_inelastic_xy = 0.0;
        PetscReal tau,mu_s,mu_d,D_c,mu_friction,T, ttau;

        evaluate_sdf(NULL,coor_qp,&phi_p);
        //if (phi_p < 0) printf("[e %d , q %d] x_qp %+1.4e , %+1.4e : phi %+1.4e \n",e,q,coor_qp[0],coor_qp[1],phi_p);
        //printf("  x_qp -> phi %+1.4e\n",phi_p);
        
        ierr = FaultSDFGetPlusMinusCoor(coor_qp,c->delta,NULL,x_plus,x_minus);CHKERRQ(ierr);
        //printf("  x_qp -> x+ %+1.4e , %+1.4e\n",x_plus[0],x_plus[1]);
        //printf("  x_qp -> x- %+1.4e , %+1.4e\n",x_minus[0],x_minus[1]);
        

        
#if 1
        /* ================================================================ */
        ierr = FaultSDFTabulateInterpolation_v2(c,LA_u,&dr_celldata[q],v_plus,v_minus);CHKERRQ(ierr);
        
        ierr = FaultSDFNormal(coor_qp,NULL,normal);CHKERRQ(ierr);
        ierr = FaultSDFTangent(coor_qp,NULL,tangent);CHKERRQ(ierr);
        
        /* Resolve velocities at delta(+,-) onto fault */
        /* [option 2] Removal of normal component method */
        /* I like this approach as it does not require a tangenet vector */
        {
          PetscReal mag_vdotn;
          
          mag_vdotn =  (v_plus[0] * normal[0] +  v_plus[1] * normal[1]);
          v_plus[0] = v_plus[0] - mag_vdotn * normal[0];
          v_plus[1] = v_plus[1] - mag_vdotn * normal[1];
          
          /* Error checking is not generalized to non-planar faults */
          if (fabs(v_plus[1]) > 1.0e-10) { /* error checking to ensure that the y component is completely removed */
            printf("  phi %+1.8e : v+_y > 0 (%+1.8e)\n",phi_p,v_plus[1]);
            printf("  |v+| %+1.8e\n",mag_vdotn);
            printf("  x_qp -> v+ %+1.8e , %+1.8e\n",v_plus[0],v_plus[1]);
            exit(1);
          }
          
          mag_vdotn =  (v_minus[0] * normal[0] +  v_minus[1] * normal[1]);
          v_minus[0] = v_minus[0] - mag_vdotn * normal[0];
          v_minus[1] = v_minus[1] - mag_vdotn * normal[1];
          
          /* Error checking is not generalized to non-planar faults */
          if (fabs(v_minus[1]) > 1.0e-10) {
            printf("  phi %+1.8e : v-_y > 0 (%+1.8e)\n",phi_p,v_minus[1]);
            printf("  |v-| %+1.8e\n",mag_vdotn);
            printf("  x_qp -> v+ %+1.8e , %+1.8e\n",v_minus[0],v_minus[1]);
            exit(1);
          }
        }
        
        /* I believe that extract the first component of the vector is perfectly valid for non-planar faults */
        Vplus  = v_plus[0];
        Vminus = v_minus[0];
        slip = Vplus - Vminus;
        
        slip = dr_celldata[q].slip;
        slip_rate = dr_celldata[q].slip_rate;
        
        //slip = (0.5*Vplus*c->delta - 0.5*Vminus*c->delta)/(2.0 * c->delta);
        
        /* ================================================================ */
#endif
        //slip = dr_celldata[q].slip;
        
        // n_i s_ij n_j
        // = ni si0 n0 + ni si1 n1
        // = n0 s00 n0 + n1 s10 n0 + n0 s01 n1 + n1 s11 n1
        sigma_n = 0.0;
        /* [option 1] */
        sigma_n += normal[0] * sigma_trial[TENS2D_XX] * normal[0];
        sigma_n += normal[1] * sigma_trial[TENS2D_XY] * normal[0];
        sigma_n += normal[0] * sigma_trial[TENS2D_XY] * normal[1];
        sigma_n += normal[1] * sigma_trial[TENS2D_YY] * normal[1];
        /* [option 2] - only valid for a horizontal fault! */
        /* option 1 and 2 are identical for horizontal fault */
        sigma_n = sigma_trial[TENS2D_YY];
        
        sigma_t = 0.0;
        /* [option 1] - uses the tangent vectors - yuck */
        sigma_t += tangent[0] * sigma_trial[TENS2D_XX] * normal[0];
        sigma_t += tangent[1] * sigma_trial[TENS2D_XY] * normal[0];
        sigma_t += tangent[0] * sigma_trial[TENS2D_XY] * normal[1];
        sigma_t += tangent[1] * sigma_trial[TENS2D_YY] * normal[1];
        /* [option 2] - only valid for a horizontal fault! */
        sigma_t = sigma_trial[TENS2D_XY];

        sigma_t = sigma_t - mu_qp * gradu[2]; // new tweak [March 16]


        dr_celldata[q].mu = 0;
        e_inelastic_xy = 0.0;
        sliding_active = PETSC_FALSE;

        if (sigma_n < 0) { /* only consider inelastic corrections if in compression */
          T = sqrt(sigma_t * sigma_t);
          
          /* Hard code linear slip weakening */
          mu_s = c->mu_s;
          mu_d = c->mu_d;
          D_c  = c->D_c;
          FricSW(&mu_friction, mu_s, mu_d, D_c, fabs(slip)); /* note the inclusion of making slip always positive */
          
          dr_celldata[q].mu = mu_friction;
          
          tau = -mu_friction * sigma_n;
          if (tau < 0) {
            printf("-mu sigma_n < 0 error\n");
            exit(1);
          }

          
          if (T > tau) {
            PetscReal factor;
            
            ierr = FaultSDFMollifer(coor_qp,dy,NULL,&factor);CHKERRQ(ierr);

            e_inelastic_xy = c->delta * (sigma_t - tau) / (mu_qp); // = du / dy
            

            /**Antiparallel condition between slip rate and critical shear */
            ttau = tau;
            if ( sigma_t < 0.0) //slip_rate=v(+)-v(-) defined following Dalguer
            {
              ttau = -tau;
              //printf(" [phi: %f, T: %f, tau: %f, Coor_x: %f,Coor_y: %f, slipe_rate: %f, sigma_t (PreBlend): %f, slip: %f, e_vec[2]: %f]\n", phi_p, T, tau, coor_qp[0], coor_qp[1], slip_rate, sigma_t, slip, e_vec[2]);
            } 
            sigma_trial[TENS2D_XY] = ttau;
            /**TP3 - Smoothing for p > 1 */
            if(c->basisorder > 1)
            {
              //ierr = PetscTanHWeighting( &sigma_t,  sigma_t, tau, phi_p , 4.*(c->basisorder)/c->delta,  0.65*c->delta); CHKERRQ(ierr);
              ierr = PetscTanHWeighting( &sigma_t,  sigma_t, ttau, phi_p , 6.5*(c->basisorder)/c->delta,  0.85*c->delta); CHKERRQ(ierr);
              sigma_trial[TENS2D_XY] = sigma_t;
            }

            //printf("  sigma_xy %+1.8e\n",sigma_vec[TENS2D_XY]);
            sliding_active = PETSC_TRUE;
            dr_celldata[q].sliding = PETSC_TRUE;
          } else {
            e_inelastic_xy = 0.0;
            slip_rate = 0;
            sliding_active = PETSC_FALSE;
          }

          
          
          // Error checking / verification that consistency conditions are approximately satisfied
          /*
          if (T > fabs(tau)) {
            if (fabs(sigma_trial[TENS2D_XY]) - fabs(tau) > 1e-7) {
              printf("  [1] |T_t| - mu |T_n| %+1.12e < = ?\n",fabs(sigma_trial[TENS2D_XY]) - fabs(tau));
            }
            {
              double a = slip_rate * fabs(sigma_trial[TENS2D_XY]) - fabs(slip_rate) * sigma_trial[TENS2D_XY];
              
              if (fabs(a) > 1.0e-10) {
                printf("  [3] dot s |T_t| - |dot s| T_t %+1.12e = 0 ?\n",slip_rate * fabs(sigma_trial[TENS2D_XY]) - fabs(slip_rate) * sigma_trial[TENS2D_XY]);
              }
            }
          }
          */
        }
      }

      /* Remove weird non-generalizable background stress state */

      if (inside_fault_region) { /* remove the initial stress state on fault */
        if (fabs(coor_qp[0]) < 1.5*1.0e3+1.0) {
          sigma_trial[TENS2D_XY] -= sigma_t_1;
          sigma_trial[TENS2D_YY] -= (-sigma_n_1); /* negative in compression */
        } else {
          sigma_trial[TENS2D_XY] -= sigma_t_0;
          sigma_trial[TENS2D_YY] -= (-sigma_n_0); /* negative in compression */
        }
      } else {
        sigma_trial[TENS2D_XY] -= sigma_t_0;
        sigma_trial[TENS2D_YY] -= (-sigma_n_0); /* negative in compression */
      }

      /* These components weren't modified in the horizontal fault case - but they might be in general */
      sigma_vec[TENS2D_XX] = sigma_trial[TENS2D_XX];
      sigma_vec[TENS2D_YY] = sigma_trial[TENS2D_YY];
      /* This component was modified in the horizontal fault case - it's likely it might also be modified in the general case */
      sigma_vec[TENS2D_XY] = sigma_trial[TENS2D_XY];
      
      fac = detJ * c->w[q];
      for (i=0; i<nbasis; i++) {
        fe[2*i  ] += -fac * (dNidx[i] * sigma_vec[TENS2D_XX] + dNidy[i] * sigma_vec[TENS2D_XY]);
        fe[2*i+1] += -fac * (dNidy[i] * sigma_vec[TENS2D_YY] + dNidx[i] * sigma_vec[TENS2D_XY]);
      }
      
    }
    ierr = VecSetValues(fl,nbasis*ndof,eldofs,fe,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(fl);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(fl);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(c->dm,fl,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(c->dm,fl,ADD_VALUES,F);CHKERRQ(ierr);
  
  ierr = VecRestoreArrayRead(vl,&LA_v);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&vl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(ul,&LA_u);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&ul);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&fl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(coor,&LA_coor);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}



PetscErrorCode Update_StressGlut2d(SpecFECtx c,Vec u,Vec v,PetscReal dt)
{
  PetscErrorCode ierr;
  PetscInt  e,nqp,q,i,nbasis,ndof;
  PetscReal e_vec[3],gradu[4],gradv[4],gradv_q[9*9][4];
  PetscInt  *element,*elnidx,*eldofs;
  PetscReal *ux,*uy,*vx,*vy,*elcoords,detJ,*fieldU,*fieldV;
  Vec       coor,ul,vl;
  const PetscReal *LA_coor,*LA_u,*LA_v;
  QPntIsotropicElastic *celldata;
  DRVar                *dr_celldata;
  
  static PetscBool beenhere = PETSC_FALSE;
  static PetscReal gmin[3],gmax[3];
  PetscReal dx,dy;
  
  if (!beenhere) {
    ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
    beenhere = PETSC_TRUE;
  }
  dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
  dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);
  
  
  
  eldofs   = c->elbuf_dofs;
  elcoords = c->elbuf_coor;
  nbasis   = c->npe;
  nqp      = c->nqp;
  ndof     = c->dofs;
  element  = c->element;
  fieldU   = c->elbuf_field2;
  fieldV   = c->elbuf_field3;
  
  ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(c->dm,&ul);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(c->dm,u,INSERT_VALUES,ul);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,u,INSERT_VALUES,ul);CHKERRQ(ierr);
  ierr = VecGetArrayRead(ul,&LA_u);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(c->dm,&vl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(c->dm,v,INSERT_VALUES,vl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,v,INSERT_VALUES,vl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(vl,&LA_v);CHKERRQ(ierr);
  
  ux = &fieldU[0];
  uy = &fieldU[nbasis];
  
  vx = &fieldV[0];
  vy = &fieldV[nbasis];
  
  for (e=0; e<c->ne; e++) {
    PetscReal x_cell[] = {0,0};
    
    /* get element -> node map */
    elnidx = &element[nbasis*e];
    
    /* generate dofs */
    for (i=0; i<nbasis; i++) {
      eldofs[2*i  ] = 2*elnidx[i];
      eldofs[2*i+1] = 2*elnidx[i]+1;
    }
    
    /* get element coordinates */
    for (i=0; i<nbasis; i++) {
      PetscInt nidx = elnidx[i];
      elcoords[2*i  ] = LA_coor[2*nidx  ];
      elcoords[2*i+1] = LA_coor[2*nidx+1];
      x_cell[0] += elcoords[2*i  ];
      x_cell[1] += elcoords[2*i+1];
    }
    x_cell[0] = x_cell[0] / ((PetscReal)nbasis);
    x_cell[1] = x_cell[1] / ((PetscReal)nbasis);
    
    /* get element displacements & velocities */
    for (i=0; i<nbasis; i++) {
      PetscInt nidx = elnidx[i];
      ux[i] = LA_u[2*nidx  ];
      uy[i] = LA_u[2*nidx+1];
      
      vx[i] = LA_v[2*nidx  ];
      vy[i] = LA_v[2*nidx+1];
    }
    
    /* compute derivatives */
    ElementEvaluateGeometry_CellWiseConstant2d(nbasis,elcoords,c->npe_1d,&detJ);
    ElementEvaluateDerivatives_CellWiseConstant2d(nqp,nbasis,elcoords,
                                                  c->npe_1d,c->dN_dxi,c->dN_deta,
                                                  c->dN_dx,c->dN_dy);
    
    /* get access to element->quadrature points */
    celldata = &c->cell_data[e];
    
    ierr = SpecFECtxGetDRCellData(c,e,&dr_celldata);CHKERRQ(ierr);
    
    for (q=0; q<c->nqp; q++) {
      PetscReal *dNidx,*dNidy;
      
      dNidx = c->dN_dx[q];
      dNidy = c->dN_dy[q];
      
      gradv[0] = gradv[1] = gradv[2] = gradv[3] = 0.0;
      for (i=0; i<nbasis; i++) {
        gradv[0] += dNidx[i] * vx[i];
        gradv[1] += dNidy[i] * vx[i];
        gradv[2] += dNidx[i] * vy[i];
        gradv[3] += dNidy[i] * vy[i];
      }
      gradv_q[q][0] = gradv[0];
      gradv_q[q][1] = gradv[1];
      gradv_q[q][2] = gradv[2];
      gradv_q[q][3] = gradv[3];
    }
    
    
    for (q=0; q<c->nqp; q++) {
      PetscReal *dNidx,*dNidy;
      PetscReal coor_qp[2];
      PetscBool inside_fault_region;
      
      dNidx = c->dN_dx[q];
      dNidy = c->dN_dy[q];
      
      /* compute strain @ quadrature point */
      /*
       e = Bu = [ d/dx  0    ][ u v ]^T
       [ 0     d/dy ]
       [ d/dy  d/dx ]
       */
      e_vec[0] = e_vec[1] = e_vec[2] = 0.0;
      for (i=0; i<nbasis; i++) {
        e_vec[0] += dNidx[i] * ux[i];
        e_vec[1] += dNidy[i] * uy[i];
        e_vec[2] += (dNidx[i] * uy[i] + dNidy[i] * ux[i]);
      }
      
      gradu[0] = gradu[1] = gradu[2] = gradu[3] = 0.0;
      gradv[0] = gradv[1] = gradv[2] = gradv[3] = 0.0;
      for (i=0; i<nbasis; i++) {
        gradu[0] += dNidx[i] * ux[i];
        gradu[1] += dNidy[i] * ux[i];
        gradu[2] += dNidx[i] * uy[i];
        gradu[3] += dNidy[i] * uy[i];
        
        gradv[0] += dNidx[i] * vx[i];
        gradv[1] += dNidy[i] * vx[i];
        gradv[2] += dNidx[i] * vy[i];
        gradv[3] += dNidy[i] * vy[i];
      }
      
      coor_qp[0] = elcoords[2*q  ];
      coor_qp[1] = elcoords[2*q+1];
      
      inside_fault_region = PETSC_FALSE;
      
      ierr = FaultSDFQuery(coor_qp,c->delta,NULL,&inside_fault_region);CHKERRQ(ierr);
      if (fabs(x_cell[1]) > c->delta) { inside_fault_region = PETSC_FALSE; }
      
      inside_fault_region = PETSC_FALSE;
      if (fabs(x_cell[1]) < c->delta && fabs(x_cell[0]) < 15.0e3) { inside_fault_region = PETSC_TRUE; }

      //if (e == 720849) {
      //  printf("e %d : xcell %+1.8e %+1.8e inside? %d\n",e,x_cell[0],x_cell[1],(int)inside_fault_region);
      //}
      
      /* Make stress glut corrections here */
      if (inside_fault_region) {
        PetscReal x_plus[2],x_minus[2],plus[2],minus[2];
        PetscReal normal[2],tangent[2],Uplus,Uminus,Vplus,Vminus,slip,slip_rate,phi_p;
        
        evaluate_sdf(NULL,coor_qp,&phi_p);
        
        ierr = FaultSDFGetPlusMinusCoor(coor_qp,c->delta,NULL,x_plus,x_minus);CHKERRQ(ierr);
        
        /* ================================================================ */
        ierr = FaultSDFTabulateInterpolation_v2(c,LA_v,&dr_celldata[q],plus,minus);CHKERRQ(ierr);
        ierr = FaultSDFNormal(coor_qp,NULL,normal);CHKERRQ(ierr);
        ierr = FaultSDFTangent(coor_qp,NULL,tangent);CHKERRQ(ierr);
        
        /* Resolve velocities at delta(+,-) onto fault */
        {
          PetscReal mag_vdotn;
          
          mag_vdotn = plus[0] * normal[0] +  plus[1] * normal[1];
          plus[0] = plus[0] - mag_vdotn * normal[0];
          plus[1] = plus[1] - mag_vdotn * normal[1];
          
          mag_vdotn = minus[0] * normal[0] +  minus[1] * normal[1];
          minus[0] = minus[0] - mag_vdotn * normal[0];
          minus[1] = minus[1] - mag_vdotn * normal[1];
        }

        Vplus  = plus[0];
        Vminus = minus[0];
        slip_rate = Vplus - Vminus;

        /* mid-point quadrature rule */
        //slip_rate = (0.5*Vplus*c->delta - 0.5*Vminus*c->delta)/(2.0*c->delta);

        
        /* ================================================================ */
        ierr = FaultSDFTabulateInterpolation_v2(c,LA_u,&dr_celldata[q],plus,minus);CHKERRQ(ierr);
        ierr = FaultSDFNormal(coor_qp,NULL,normal);CHKERRQ(ierr);
        ierr = FaultSDFTangent(coor_qp,NULL,tangent);CHKERRQ(ierr);
        
        /* Resolve displacement at delta(+,-) onto fault */
        {
          PetscReal mag_vdotn;
          
          mag_vdotn = plus[0] * normal[0] +  plus[1] * normal[1];
          plus[0] = plus[0] - mag_vdotn * normal[0];
          plus[1] = plus[1] - mag_vdotn * normal[1];
          
          mag_vdotn = minus[0] * normal[0] +  minus[1] * normal[1];
          minus[0] = minus[0] - mag_vdotn * normal[0];
          minus[1] = minus[1] - mag_vdotn * normal[1];
        }
        
        Uplus  = plus[0];
        Uminus = minus[0];
        slip = Uplus - Uminus;
        
        /* mid-point quadrature rule */
        //slip = (0.5*Uplus*c->delta - 0.5*Uminus*c->delta)/(2.0*c->delta);

        /* ================================================================ */
        
        
        if (dr_celldata[q].sliding) {
          /*
          dr_celldata[q].slip_rate = slip_rate;
          dr_celldata[q].slip      += slip_rate * dt;

          dr_celldata[q].slip_rate  = (slip - dr_celldata[q].slip)/dt;
          dr_celldata[q].slip       = slip;
          */

          dr_celldata[q].slip_rate = slip_rate;
          dr_celldata[q].slip      = slip;
        }
        
      }
      
    }
  }
  
  ierr = VecRestoreArrayRead(vl,&LA_v);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&vl);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(ul,&LA_u);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&ul);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


PetscErrorCode AssembleBilinearForm_Mass2d(SpecFECtx c,Vec A)
{
  PetscErrorCode ierr;
  PetscInt  e,index,q,i,nbasis,ndof;
  PetscInt  *element,*elnidx,*eldofs;
  PetscReal *elcoords,*Me,detJ;
  Vec       coor;
  const PetscReal *LA_coor;
  QPntIsotropicElastic *celldata;
  
  ierr = VecZeroEntries(A);CHKERRQ(ierr);
  
  eldofs   = c->elbuf_dofs;
  elcoords = c->elbuf_coor;
  nbasis   = c->npe;
  ndof     = c->dofs;
  Me       = c->elbuf_field;
  element  = c->element;
  
  ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  for (e=0; e<c->ne; e++) {
    /* get element -> node map */
    elnidx = &element[nbasis*e];
    
    /* generate dofs */
    for (i=0; i<nbasis; i++) {
      
      eldofs[2*i  ] = 2*elnidx[i];
      eldofs[2*i+1] = 2*elnidx[i]+1;
    }
    
    /* get element coordinates */
    for (i=0; i<nbasis; i++) {
      PetscInt nidx = elnidx[i];
      elcoords[2*i  ] = LA_coor[2*nidx  ];
      elcoords[2*i+1] = LA_coor[2*nidx+1];
    }
    
    ElementEvaluateGeometry_CellWiseConstant2d(nbasis,elcoords,c->npe_1d,&detJ);
    
    /* get access to element->quadrature points */
    celldata = &c->cell_data[e];

    for (q=0; q<nbasis; q++) {
      PetscReal            fac,Me_ii;
      
      fac = detJ * c->w[q];
      
      Me_ii = fac * (celldata->rho);
      
      /* \int u0v0 dV */
      index = 2*q;
      Me[index] = Me_ii;
      
      /* \int u1v1 dV */
      index = 2*q + 1;
      Me[index] = Me_ii;
    }
    ierr = VecSetValuesLocal(A,nbasis*ndof,eldofs,Me,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(A);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(A);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(coor,&LA_coor);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode ElastoDynamicsConvertLame2Velocity(PetscReal rho,PetscReal mu,PetscReal lambda,PetscReal *Vs,PetscReal *Vp)
{
  if (Vs) { *Vs = PetscSqrtReal(mu/rho); }
  if (Vp) { *Vp = PetscSqrtReal( (lambda + 2.0*mu)/rho); }
  
  PetscFunctionReturn(0);
}

PetscErrorCode ElastoDynamicsComputeTimeStep_2d(SpecFECtx ctx,PetscReal *_dt)
{
  PetscInt e,q,order;
  PetscReal dt_min,dt_min_g,polynomial_fac;
  //QPntIsotropicElastic *qpdata;
  PetscReal gmin[3],gmax[3],min_el_r,dx,dy;
  PetscErrorCode ierr;
  QPntIsotropicElastic *celldata;
  
  *_dt = PETSC_MAX_REAL;
  dt_min = PETSC_MAX_REAL;
  
  order = ctx->basisorder;
  polynomial_fac = 1.0 / (2.0 * (PetscReal)order + 1.0);
  
  ierr = DMGetBoundingBox(ctx->dm,gmin,gmax);CHKERRQ(ierr);
  dx = (gmax[0] - gmin[0])/((PetscReal)ctx->mx_g);
  dy = (gmax[1] - gmin[1])/((PetscReal)ctx->my_g);
  
  min_el_r = dx;
  min_el_r = PetscMin(min_el_r,dy);
  
  /* find smallest dx across the element in local coordinates */
  {
    PetscInt  n;
    PetscReal sep2min,sep2;
    
    sep2min = 1.0e32;
    for (n=0; n<ctx->npe_1d-1; n++) {
      sep2 = PetscAbsReal(ctx->xi1d[n+1] - ctx->xi1d[n]);
      /*printf(" xi %+1.4e [n] : xi %+1.4e [n+1] : delta_xi %+1.6e\n",ctx->xi1d[n],ctx->xi1d[n+1],sep2); */
      if (sep2 < sep2min) {
        sep2min = sep2;
      }
    }
    
    polynomial_fac = 1.0;
    min_el_r = min_el_r * ( sep2min / 2.0 ); /* the factor 2.0 here is associated with the size of the element in the local coordinate system xi \in [-1,+1] */
  }
  
  

  
  for (e=0; e<ctx->ne; e++) {
    PetscReal max_el_Vp,value;
    
    /* get max Vp for element */
    max_el_Vp = PETSC_MIN_REAL;
    
    /* get access to element->quadrature points */
    celldata = &ctx->cell_data[e];
    
    for (q=0; q<ctx->nqp; q++) {
      PetscReal qp_rho,qp_mu,qp_lambda,qp_Vp;
      
      qp_rho    = celldata->rho;
      qp_mu     = celldata->mu;
      qp_lambda = celldata->lambda;
      
      ierr = ElastoDynamicsConvertLame2Velocity(qp_rho,qp_mu,qp_lambda,0,&qp_Vp);CHKERRQ(ierr);
      
      max_el_Vp = PetscMax(max_el_Vp,qp_Vp);
    }
    
    value = polynomial_fac * 1.0 * min_el_r / max_el_Vp;
    
    dt_min = PetscMin(dt_min,value);
  }
  ierr = MPI_Allreduce(&dt_min,&dt_min_g,1,MPIU_REAL,MPIU_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  *_dt = dt_min_g;
  
  PetscFunctionReturn(0);
}

PetscErrorCode RecordUV(SpecFECtx c,PetscReal time,PetscReal xr[],Vec u,Vec v)
{
  FILE *fp = NULL;
  PetscReal gmin[3],gmax[3],dx,dy,sep2min,sep2;
  const PetscReal *LA_u,*LA_v,*LA_c;
  Vec coor;
  static PetscBool beenhere = PETSC_FALSE;
  PetscErrorCode ierr;
  PetscInt ei,ej,n,nid,eid,*element,*elbasis;
  static char filename[PETSC_MAX_PATH_LEN];
  
  if (c->size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Needs updating to support MPI");
  if (!beenhere) {
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"receiverCP-%Dx%D-p%D.dat",c->mx_g,c->my_g,c->basisorder);CHKERRQ(ierr);
  }
  
  ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
  dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
  ei = (xr[0] - gmin[0])/dx; /* todo - needs to be sub-domain gmin */
  
  dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);
  ej = (xr[1] - gmin[1])/dy; /* todo - needs to be sub-domain gmin */
  
  eid = ei + ej * c->mx;
  
  /* get element -> node map */
  element = c->element;
  elbasis = &element[c->npe*eid];
  
  ierr = DMGetCoordinates(c->dm,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_c);CHKERRQ(ierr);
  
  // find closest //
  sep2min = 1.0e32;
  nid = -1;
  for (n=0; n<c->npe; n++) {
    sep2  = (xr[0]-LA_c[2*elbasis[n]])*(xr[0]-LA_c[2*elbasis[n]]);
    sep2 += (xr[1]-LA_c[2*elbasis[n]+1])*(xr[1]-LA_c[2*elbasis[n]+1]);
    if (sep2 < sep2min) {
      nid = elbasis[n];
      sep2min = sep2;
    }
  }
  
  if (!beenhere) {
    fp = fopen(filename,"w");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
    fprintf(fp,"# SpecFECtx meta data\n");
    fprintf(fp,"#   mx %d : my %d : basis order %d\n",c->mx_g,c->my_g,c->basisorder);
    fprintf(fp,"# Receiver meta data\n");
    fprintf(fp,"#   + receiver location: x,y %+1.8e %+1.8e\n",xr[0],xr[1]);
    fprintf(fp,"#   + takes displ/velo from basis nearest to requested receiver location\n");
    fprintf(fp,"#   + receiver location: x,y %+1.8e %+1.8e --mapped to nearest node --> %+1.8e %+1.8e\n",xr[0],xr[1],LA_c[2*nid],LA_c[2*nid+1]);
    fprintf(fp,"# Time series header\n");
    fprintf(fp,"#   time ux uy vx vy\n");
    beenhere = PETSC_TRUE;
  } else {
    fp = fopen(filename,"a");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
  }
  
  ierr = VecGetArrayRead(u,&LA_u);CHKERRQ(ierr);
  ierr = VecGetArrayRead(v,&LA_v);CHKERRQ(ierr);
  
  fprintf(fp,"%1.4e %+1.8e %+1.8e %+1.8e %+1.8e\n",time,LA_u[2*nid],LA_u[2*nid+1],LA_v[2*nid],LA_v[2*nid+1]);
  
  ierr = VecRestoreArrayRead(v,&LA_v);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(u,&LA_u);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(coor,&LA_c);CHKERRQ(ierr);
  
  fclose(fp);
  
  PetscFunctionReturn(0);
}

PetscErrorCode RecordUV_interp(SpecFECtx c,PetscReal time,PetscReal xr[],Vec u,Vec v)
{
  FILE *fp = NULL;
  PetscReal gmin[3],gmax[3],dx,dy,ur[2],vr[2];
  const PetscReal *LA_u,*LA_v;
  static PetscBool beenhere = PETSC_FALSE;
  PetscErrorCode ierr;
  PetscInt k,ei,ej,eid,*element,*elbasis;
  static PetscReal N[400];
  static char filename[PETSC_MAX_PATH_LEN];
  
  if (c->size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Needs updating to support MPI");
  if (!beenhere) {
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"receiver-%Dx%D-p%D.dat",c->mx_g,c->my_g,c->basisorder);CHKERRQ(ierr);
  }
  
  if (!beenhere) {
    fp = fopen(filename,"w");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
    fprintf(fp,"# SpecFECtx meta data\n");
    fprintf(fp,"#   mx %d : my %d : basis order %d\n",c->mx_g,c->my_g,c->basisorder);
    fprintf(fp,"# Receiver meta data\n");
    fprintf(fp,"#   + receiver location: x,y %+1.8e %+1.8e\n",xr[0],xr[1]);
    fprintf(fp,"#   + records displ/velo at requested receiver location through interpolating the FE solution\n");
    fprintf(fp,"# Time series header\n");
    fprintf(fp,"#   time ux uy vx vy\n");
  } else {
    fp = fopen(filename,"a");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
  }
  
  /* get containing element */
  ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
  dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
  ei = (xr[0] - gmin[0])/dx; /* todo - needs to be sub-domain gmin */
  
  dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);
  ej = (xr[1] - gmin[1])/dy;
  
  eid = ei + ej * c->mx;
  
  /* get element -> node map */
  element = c->element;
  elbasis = &element[c->npe*eid];
  
  if (!beenhere) {
    PetscInt nbasis,i,j;
    PetscReal **N_s1,**N_s2,xri[2],xi,eta,x0,y0;
    const PetscReal *LA_c;
    Vec coor;
    
    /* compute xi,eta */
    ierr = DMGetCoordinates(c->dm,&coor);CHKERRQ(ierr);
    
    x0 = gmin[0] + ei*dx; /* todo - needs to be sub-domain gmin */
    y0 = gmin[1] + ej*dy;
    
    // (xi - (-1))/2 = (x - x0)/dx
    xi = 2.0*(xr[0] - x0)/dx - 1.0;
    eta = 2.0*(xr[1] - y0)/dy - 1.0;
    
    /* compute basis */
    ierr = TabulateBasis1d_CLEGENDRE(1,&xi,c->basisorder,&nbasis,&N_s1);CHKERRQ(ierr);
    ierr = TabulateBasis1d_CLEGENDRE(1,&eta,c->basisorder,&nbasis,&N_s2);CHKERRQ(ierr);
    
    k = 0;
    for (j=0; j<c->npe_1d; j++) {
      for (i=0; i<c->npe_1d; i++) {
        N[k] = N_s1[0][i] * N_s2[0][j];
        k++;
      }
    }
    
    ierr = VecGetArrayRead(coor,&LA_c);CHKERRQ(ierr);
    
    xri[0] = xri[1] = 0.0;
    for (k=0; k<c->npe; k++) {
      PetscInt nid = elbasis[k];
      
      xri[0] += N[k] * LA_c[2*nid+0];
      xri[1] += N[k] * LA_c[2*nid+1];
    }
    
    
    PetscPrintf(PETSC_COMM_SELF,"# receiver location: x,y %+1.8e %+1.8e -- interpolated coordinate --> %+1.8e %+1.8e\n",xr[0],xr[1],xri[0],xri[1]);
    
    ierr = VecRestoreArrayRead(coor,&LA_c);CHKERRQ(ierr);
    ierr = PetscFree(N_s1[0]);CHKERRQ(ierr);
    ierr = PetscFree(N_s1);CHKERRQ(ierr);
    ierr = PetscFree(N_s2[0]);CHKERRQ(ierr);
    ierr = PetscFree(N_s2);CHKERRQ(ierr);
  }
  
  
  ierr = VecGetArrayRead(u,&LA_u);CHKERRQ(ierr);
  ierr = VecGetArrayRead(v,&LA_v);CHKERRQ(ierr);
  
  ur[0] = ur[1] = vr[0] = vr[1] = 0.0;
  for (k=0; k<c->npe; k++) {
    PetscInt nid = elbasis[k];
    
    ur[0] += N[k] * LA_u[2*nid+0];
    ur[1] += N[k] * LA_u[2*nid+1];
    
    vr[0] += N[k] * LA_v[2*nid+0];
    vr[1] += N[k] * LA_v[2*nid+1];
  }
  
  fprintf(fp,"%1.4e %+1.8e %+1.8e %+1.8e %+1.8e\n",time,ur[0],ur[1],vr[0],vr[1]);
  
  ierr = VecRestoreArrayRead(v,&LA_v);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(u,&LA_u);CHKERRQ(ierr);
  
  beenhere = PETSC_TRUE;
  fclose(fp);
  
  PetscFunctionReturn(0);
}

PetscErrorCode RecordUVA_MultipleStations_NearestGLL_SEQ(SpecFECtx c,PetscReal time,PetscInt nr,PetscReal xr[],Vec u,Vec v,Vec a)
{
  FILE             *fp = NULL;
  const PetscReal  *LA_u,*LA_v,*LA_a;
  static PetscBool beenhere = PETSC_FALSE;
  static char      filename[PETSC_MAX_PATH_LEN];
  static PetscInt  *nid_list = NULL;
  PetscInt         r;
  PetscErrorCode   ierr;

  
  if (c->size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Supports sequential only");
  if (!beenhere) {
    const PetscReal *LA_c;
    Vec coor;
    PetscReal gmin[3],gmax[3],dx,dy,sep2min,sep2;
    PetscInt ei,ej,n,nid,eid,*element,*elbasis;
    
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"closestqpsource-receiverCP-uva-%Dx%D-p%D.dat",c->mx_g,c->my_g,c->basisorder);CHKERRQ(ierr);
    ierr = PetscMalloc1(nr,&nid_list);CHKERRQ(ierr);
    
    ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
    ierr = DMGetCoordinates(c->dm,&coor);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coor,&LA_c);CHKERRQ(ierr);
    
    for (r=0; r<nr; r++) {
      
      if (xr[2*r+0] < gmin[0]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, x-coordinate (%+1.4e) < min(domain).x (%+1.4e)",r,xr[2*r+0],gmin[0]);
      if (xr[2*r+1] < gmin[1]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, y-coordinate (%+1.4e) < min(domain).y (%+1.4e)",r,xr[2*r+1],gmin[1]);
      if (xr[2*r+0] > gmax[0]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, x-coordinate (%+1.4e) > max(domain).x (%+1.4e)",r,xr[2*r+0],gmax[0]);
      if (xr[2*r+1] > gmax[1]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, y-coordinate (%+1.4e) > max(domain).y (%+1.4e)",r,xr[2*r+1],gmax[1]);
      
      dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
      ei = (xr[2*r+0] - gmin[0])/dx;
      if (ei == c->mx_g) ei--;
      
      dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);
      ej = (xr[2*r+1] - gmin[1])/dy;
      if (ej == c->my_g) ej--;
      
      eid = ei + ej * c->mx_g;
    
      /* get element -> node map */
      element = c->element;
      elbasis = &element[c->npe*eid];
    
      // find closest //
      sep2min = 1.0e32;
      nid = -1;
      for (n=0; n<c->npe; n++) {
        sep2  = (xr[2*r+0]-LA_c[2*elbasis[n]])*(xr[2*r+0]-LA_c[2*elbasis[n]]);
        sep2 += (xr[2*r+1]-LA_c[2*elbasis[n]+1])*(xr[2*r+1]-LA_c[2*elbasis[n]+1]);
        if (sep2 < sep2min) {
          nid = elbasis[n];
          sep2min = sep2;
        }
      }
      nid_list[r] = nid;
    }
  
    fp = fopen(filename,"w");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
    fprintf(fp,"# SpecFECtx meta data\n");
    fprintf(fp,"#   mx %d : my %d : basis order %d\n",c->mx_g,c->my_g,c->basisorder);
    fprintf(fp,"# Receiver meta data\n");
    fprintf(fp,"#   + number receiver locations: %d\n",nr);
    fprintf(fp,"#   + takes displ/velo/accel from basis nearest to requested receiver location\n");
    for (r=0; r<nr; r++) {
      fprintf(fp,"#   + receiver location [%d]: x,y %+1.8e %+1.8e\n",r,xr[2*r+0],xr[2*r+1]);
      fprintf(fp,"#   +   mapped to nearest node --> %+1.8e %+1.8e\n",LA_c[2*nid_list[r]],LA_c[2*nid_list[r]+1]);
    }
    fprintf(fp,"# Time series header <field>(<column index>)\n");
    fprintf(fp,"#   time(1)\n");
    for (r=0; r<nr; r++) {
      PetscInt offset = 1 + r*6; /* 1 is for time */
      
      fprintf(fp,"#     ux(%d) uy(%d) vx(%d) vy(%d) ax(%d) ay(%d) -> station [%d]\n",offset+1,offset+2,offset+3,offset+4,offset+5,offset+6,r);
    }
    ierr = VecRestoreArrayRead(coor,&LA_c);CHKERRQ(ierr);
    beenhere = PETSC_TRUE;
  } else {
    fp = fopen(filename,"a");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
  }
  
  ierr = VecGetArrayRead(u,&LA_u);CHKERRQ(ierr);
  ierr = VecGetArrayRead(v,&LA_v);CHKERRQ(ierr);
  ierr = VecGetArrayRead(a,&LA_a);CHKERRQ(ierr);
  
  fprintf(fp,"%1.4e",time);
  for (r=0; r<nr; r++) {
    fprintf(fp," %+1.8e %+1.8e %+1.8e %+1.8e %+1.8e %+1.8e",LA_u[2*nid_list[r]],LA_u[2*nid_list[r]+1],LA_v[2*nid_list[r]],LA_v[2*nid_list[r]+1],LA_a[2*nid_list[r]],LA_a[2*nid_list[r]+1]);
  }
  fprintf(fp,"\n");
  
  ierr = VecRestoreArrayRead(a,&LA_a);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(v,&LA_v);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(u,&LA_u);CHKERRQ(ierr);
  
  fclose(fp);
  
  PetscFunctionReturn(0);
}

PetscErrorCode RecordUVA_MultipleStations_NearestGLL_MPI(SpecFECtx c,PetscReal time,PetscInt nr,PetscReal xr[],Vec u,Vec v,Vec a)
{
  FILE             *fp = NULL;
  const PetscReal  *LA_u,*LA_v,*LA_a,*LA_c;
  static PetscBool beenhere = PETSC_FALSE;
  static char      filename[PETSC_MAX_PATH_LEN];
  static PetscInt  *nid_list = NULL;
  static PetscInt  *eid_list = NULL;
  static PetscInt  *gll_list = NULL;
  static PetscInt  nr_local = 0;
  PetscInt         r,k;
  Vec              lu,lv,la,coor;
  PetscErrorCode   ierr;
  
  
  if (!beenhere) {
    PetscReal       gmin[3],gmax[3],gmin_domain[3],gmax_domain[3],dx,dy,sep2min,sep2;
    PetscInt        ei,ej,n,nid,gllid,eid,*element,*elbasis;
    
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"closestqpsource-receiverCP-uva-%Dx%D-p%D-rank%d.dat",c->mx_g,c->my_g,c->basisorder,(int)c->rank);CHKERRQ(ierr);
    ierr = PetscMalloc1(nr,&nid_list);CHKERRQ(ierr);
    ierr = PetscMalloc1(nr,&eid_list);CHKERRQ(ierr);
    ierr = PetscMalloc1(nr,&gll_list);CHKERRQ(ierr);
    for (r=0; r<nr; r++) {
      nid_list[r] = -1;
      eid_list[r] = -1;
      gll_list[r] = -1;
    }

    ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
    ierr = SpecFECtxGetLocalBoundingBox(c,gmin_domain,gmax_domain);CHKERRQ(ierr);
    
    ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coor,&LA_c);CHKERRQ(ierr);
    
    for (r=0; r<nr; r++) {
      int count,recv_count;
      PetscBool receiver_found = PETSC_TRUE;
      int rank,rank_min_g;
      
      if (xr[2*r+0] < gmin[0]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, x-coordinate (%+1.4e) < min(domain).x (%+1.4e)",r,xr[2*r+0],gmin[0]);
      if (xr[2*r+1] < gmin[1]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, y-coordinate (%+1.4e) < min(domain).y (%+1.4e)",r,xr[2*r+1],gmin[1]);
      if (xr[2*r+0] > gmax[0]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, x-coordinate (%+1.4e) > max(domain).x (%+1.4e)",r,xr[2*r+0],gmax[0]);
      if (xr[2*r+1] > gmax[1]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, y-coordinate (%+1.4e) > max(domain).y (%+1.4e)",r,xr[2*r+1],gmax[1]);
      
      if (xr[2*r+0] < gmin_domain[0]) receiver_found = PETSC_FALSE;
      if (xr[2*r+1] < gmin_domain[1]) receiver_found = PETSC_FALSE;
      if (xr[2*r+0] > gmax_domain[0]) receiver_found = PETSC_FALSE;
      if (xr[2*r+1] > gmax_domain[1]) receiver_found = PETSC_FALSE;
      
      dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
      ei = (xr[2*r+0] - gmin_domain[0])/dx;
      if (ei == c->mx) ei--;
      
      dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);
      ej = (xr[2*r+1] - gmin_domain[1])/dy;
      if (ej == c->my) ej--;
      
      if (ei < 0) receiver_found = PETSC_FALSE;
      if (ej < 0) receiver_found = PETSC_FALSE;
      
      if (ei > c->mx) receiver_found = PETSC_FALSE;
      if (ej > c->my) receiver_found = PETSC_FALSE;
      
      nid = -1;
      gllid = -1;
      if (receiver_found) {
        eid = ei + ej * c->mx;
        
        /* get element -> node map */
        element = c->element;
        elbasis = &element[c->npe*eid];
        
        // find closest //
        sep2min = 1.0e32;
        for (n=0; n<c->npe; n++) {
          sep2  = (xr[2*r+0]-LA_c[2*elbasis[n]])*(xr[2*r+0]-LA_c[2*elbasis[n]]);
          sep2 += (xr[2*r+1]-LA_c[2*elbasis[n]+1])*(xr[2*r+1]-LA_c[2*elbasis[n]+1]);
          if (sep2 < sep2min) {
            nid = elbasis[n];
            gllid = n;
            sep2min = sep2;
          }
        }
      }
      
      /* check for duplicates */
      count = 0;
      if (receiver_found) {
        count = 1;
      }
      ierr = MPI_Allreduce(&count,&recv_count,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
      
      if (recv_count == 0) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"A receiver was defined but no rank claimed it");
      
      if (recv_count > 1) {
        /* resolve duplicates */
        
        rank = (int)c->rank;
        if (!receiver_found) {
          rank = (int)c->size;
        }
        ierr = MPI_Allreduce(&rank,&rank_min_g,1,MPI_INT,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
        if (rank == rank_min_g) {
          PetscPrintf(PETSC_COMM_SELF,"[RecordUVA]  + Multiple ranks located receiver (%+1.4e,%+1.4e) - rank %d claiming ownership\n",xr[2*r+0],xr[2*r+1],rank_min_g);
        }
        
        /* mark non-owning ranks as not claiming source */
        if (rank != rank_min_g) {
          receiver_found = PETSC_FALSE;
        }
      }

      if (receiver_found) {
        nid_list[r] = nid;
        eid_list[r] = eid;
        gll_list[r] = gllid;
        nr_local++;
      }
    }
    
    fp = fopen(filename,"w");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
    fprintf(fp,"# SpecFECtx meta data\n");
    fprintf(fp,"#   mx %d : my %d : basis order %d\n",c->mx_g,c->my_g,c->basisorder);
    fprintf(fp,"# Receiver meta data\n");
    fprintf(fp,"#   + number receiver locations: %d\n",nr);
    fprintf(fp,"#   + number receiver locations <local>: %d\n",nr_local);
    fprintf(fp,"#   + takes displ/velo/accel from basis nearest to requested receiver location\n");
    for (r=0; r<nr; r++) {
      if (nid_list[r] == -1) { continue; }
      fprintf(fp,"#   + receiver location [%d]: x,y %+1.8e %+1.8e\n",r,xr[2*r+0],xr[2*r+1]);
      fprintf(fp,"#   +   mapped to nearest node --> %+1.8e %+1.8e\n",LA_c[2*nid_list[r]],LA_c[2*nid_list[r]+1]);
    }

    if (nr_local != 0) {
      PetscInt count = 0;

      fprintf(fp,"# Time series header <field>(<column index>)\n");
      fprintf(fp,"#   time(1)\n");

    
      for (r=0; r<nr; r++) {
        PetscInt offset;
        
        if (nid_list[r] == -1) { continue; }

        offset = 1 + count*7; /* 1 is for time */
        
        fprintf(fp,"#     ux(%d) uy(%d) vx(%d) vy(%d) ax(%d) ay(%d) curl(v) (%d)-> station [%d]\n",offset+1,offset+2,offset+3,offset+4,offset+5,offset+6,offset+7,r);
        count++;
      }
    } else {
      fprintf(fp,"# <note> No receivers found on this sub-domain\n");
      fprintf(fp,"# <note> This file will remain empty\n");
    }

    ierr = VecRestoreArrayRead(coor,&LA_c);CHKERRQ(ierr);
    
    fclose(fp);
    fp = NULL;
  }
  
  if (!beenhere) {
    char metafname[PETSC_MAX_PATH_LEN];
    FILE *fp_meta = NULL;
    int *owned,*owned_g;
    
    ierr = PetscSNPrintf(metafname,PETSC_MAX_PATH_LEN-1,"closestqpsource-receiverCP-uva-%Dx%D-p%D.mpimeta",c->mx_g,c->my_g,c->basisorder);CHKERRQ(ierr);
    
    ierr = PetscMalloc1(nr,&owned);CHKERRQ(ierr);
    ierr = PetscMalloc1(nr,&owned_g);CHKERRQ(ierr);
    for (r=0; r<nr; r++) {
      owned[r] = -1;
      if (nid_list[r] != -1) { owned[r] = (int)c->rank; }
    }
    ierr = MPI_Allreduce(owned,owned_g,nr,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
    
    if (c->rank == 0) {
      fp_meta = fopen(metafname,"w");
      if (!fp_meta) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",metafname);
    
      fprintf(fp_meta,"# SpecFECtx parallel/MPI meta data\n");
      fprintf(fp_meta,"#   mx %d : my %d : basis order %d\n",c->mx_g,c->my_g,c->basisorder);
      fprintf(fp_meta,"# Receiver meta data\n");
      fprintf(fp_meta,"#   + number receiver locations: %d\n",nr);
      for (r=0; r<nr; r++) {
        fprintf(fp_meta,"#   + receiver [%d]: mapped to MPI rank %d\n",r,owned_g[r]);
      }
      
      fclose(fp_meta);
    }
    
    ierr = PetscFree(owned_g);CHKERRQ(ierr);
    ierr = PetscFree(owned);CHKERRQ(ierr);
  }
  
  beenhere = PETSC_TRUE;

  ierr = DMGetCoordinatesLocal(c->dm,&coor);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coor,&LA_c);CHKERRQ(ierr);
  
  ierr = DMGetLocalVector(c->dm,&lu);CHKERRQ(ierr);
  ierr = DMGetLocalVector(c->dm,&lv);CHKERRQ(ierr);
  ierr = DMGetLocalVector(c->dm,&la);CHKERRQ(ierr);
  
  ierr = DMGlobalToLocalBegin(c->dm,u,INSERT_VALUES,lu);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,u,INSERT_VALUES,lu);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(c->dm,v,INSERT_VALUES,lv);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,v,INSERT_VALUES,lv);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(c->dm,a,INSERT_VALUES,la);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(c->dm,a,INSERT_VALUES,la);CHKERRQ(ierr);
  
  ierr = VecGetArrayRead(lu,&LA_u);CHKERRQ(ierr);
  ierr = VecGetArrayRead(lv,&LA_v);CHKERRQ(ierr);
  ierr = VecGetArrayRead(la,&LA_a);CHKERRQ(ierr);
  
  if (nr_local != 0) {
    fp = fopen(filename,"a");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
  
    fprintf(fp,"%1.4e",time);
  
    for (r=0; r<nr; r++) {
      PetscInt eidx,gllidx;
      PetscReal *elcoor,*elvelocity;
      PetscReal dvxdy,dvydx,curl;
      PetscReal *grad_N_xi[2];
      PetscReal *grad_N_x[2];
      
      if (eid_list[r] == -1) { continue; }
      
      /* write the components of u,v,a */
      fprintf(fp," %+1.8e %+1.8e %+1.8e %+1.8e %+1.8e %+1.8e",LA_u[2*nid_list[r]],LA_u[2*nid_list[r]+1],LA_v[2*nid_list[r]],LA_v[2*nid_list[r]+1],LA_a[2*nid_list[r]],LA_a[2*nid_list[r]+1]);
      
      /* compute and write the k^th component of the curl(v) */
      eidx   = eid_list[r];
      gllidx = gll_list[r];
      
      grad_N_xi[0] = c->dN_dxi[gllidx];
      grad_N_xi[1] = c->dN_deta[gllidx];
      
      grad_N_x[0] = c->dN_dx[gllidx];
      grad_N_x[1] = c->dN_dy[gllidx];
      
      elcoor = c->elbuf_field;
      elvelocity = c->elbuf_field2;
      
      for (k=0; k<c->npe; k++) {
        PetscInt basisid = c->element[c->npe*eidx + k];
        
        elcoor[2*k+0] = LA_c[2*basisid + 0];
        elcoor[2*k+1] = LA_c[2*basisid + 1];
        
        elvelocity[2*k+0] = LA_v[2*basisid + 0];
        elvelocity[2*k+1] = LA_v[2*basisid + 1];
      }
      
      ElementEvaluateDerivatives_CellWiseConstant2d(1,c->npe,elcoor,c->npe_1d,&grad_N_xi[0],&grad_N_xi[1],&grad_N_x[0],&grad_N_x[1]);
      
      dvxdy = 0.0;
      dvydx = 0.0;
      for (k=0; k<c->npe; k++) {
        PetscReal vx,vy;
        
        vx = elvelocity[2*k+0];
        vy = elvelocity[2*k+1];
        
        dvxdy += grad_N_x[1][k] * vx;
        dvydx += grad_N_x[0][k] * vy;
      }
      curl = dvydx - dvxdy;
      fprintf(fp," %+1.8e",curl);
    }
    
    fprintf(fp,"\n");

    if (fp) {
      fclose(fp);
      fp = NULL;
    }
  }
  
  ierr = VecRestoreArrayRead(a,&LA_a);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(v,&LA_v);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(u,&LA_u);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(coor,&LA_c);CHKERRQ(ierr);

  ierr = DMRestoreLocalVector(c->dm,&la);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&lv);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(c->dm,&lu);CHKERRQ(ierr);

  
  PetscFunctionReturn(0);
}

PetscErrorCode RecordUVA_MultipleStations_NearestGLL(SpecFECtx c,PetscReal time,PetscInt nr,PetscReal xr[],Vec u,Vec v,Vec a)
{
  PetscErrorCode ierr;
  
  //if (c->size == 1) {
  //ierr = RecordUVA_MultipleStations_NearestGLL_SEQ(c,time,nr,xr,u,v,a);CHKERRQ(ierr);
  //} else {
  ierr = RecordUVA_MultipleStations_NearestGLL_MPI(c,time,nr,xr,u,v,a);CHKERRQ(ierr);
  //}
  PetscFunctionReturn(0);
}

PetscErrorCode RecordDRVar_MultipleStations_NearestGLL_SEQ(SpecFECtx c,PetscReal time,PetscInt nr,PetscReal xr[])
{
  FILE             *fp = NULL;
  static PetscBool beenhere = PETSC_FALSE;
  static char      filename[PETSC_MAX_PATH_LEN];
  static PetscInt  *qid_list = NULL;
  static PetscInt  *nid_list = NULL;
  static PetscInt  *eid_list = NULL;
  PetscInt         r;
  PetscErrorCode   ierr;
  
  
  if (c->size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Supports sequential only");
  if (!beenhere) {
    const PetscReal *LA_c;
    Vec coor;
    PetscReal gmin[3],gmax[3],dx,dy,sep2min,sep2;
    PetscInt ei,ej,n,qid,nid,eid,*element,*elbasis;
    
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"receiverCP-dr-%Dx%D-p%D.dat",c->mx_g,c->my_g,c->basisorder);CHKERRQ(ierr);
    ierr = PetscMalloc1(nr,&qid_list);CHKERRQ(ierr);
    ierr = PetscMalloc1(nr,&nid_list);CHKERRQ(ierr);
    ierr = PetscMalloc1(nr,&eid_list);CHKERRQ(ierr);
    
    ierr = DMGetBoundingBox(c->dm,gmin,gmax);CHKERRQ(ierr);
    ierr = DMGetCoordinates(c->dm,&coor);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coor,&LA_c);CHKERRQ(ierr);
    
    for (r=0; r<nr; r++) {
      
      if (xr[2*r+0] < gmin[0]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, x-coordinate (%+1.4e) < min(domain).x (%+1.4e)",r,xr[2*r+0],gmin[0]);
      if (xr[2*r+1] < gmin[1]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, y-coordinate (%+1.4e) < min(domain).y (%+1.4e)",r,xr[2*r+1],gmin[1]);
      if (xr[2*r+0] > gmax[0]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, x-coordinate (%+1.4e) > max(domain).x (%+1.4e)",r,xr[2*r+0],gmax[0]);
      if (xr[2*r+1] > gmax[1]) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_USER,"Receiver %D, y-coordinate (%+1.4e) > max(domain).y (%+1.4e)",r,xr[2*r+1],gmax[1]);
      
      dx = (gmax[0] - gmin[0])/((PetscReal)c->mx_g);
      ei = (xr[2*r+0] - gmin[0])/dx;
      if (ei == c->mx_g) ei--;
      
      dy = (gmax[1] - gmin[1])/((PetscReal)c->my_g);
      ej = (xr[2*r+1] - gmin[1])/dy;
      if (ej == c->my_g) ej--;
      
      eid = ei + ej * c->mx_g;
      eid_list[r] = eid;
      
      
      {
        PetscBool inside;
        PetscInt e,i;
        
        sep2min = 1.0e32;
        element = c->element;
        for (e=0; e<c->ne; e++) {
          PetscReal xcell[] = {0,0};
          
          elbasis = &element[c->npe*e];
          for (i=0; i<c->npe; i++) {
            xcell[0] += LA_c[2*elbasis[i]+0];
            xcell[1] += LA_c[2*elbasis[i]+1];
          }
          xcell[0] = xcell[0] / ((PetscReal)c->npe);
          xcell[1] = xcell[1] / ((PetscReal)c->npe);

          ierr = FaultSDFQuery(xcell,c->delta,NULL,&inside);CHKERRQ(ierr);
          if (inside) {
            sep2  = (xr[2*r+0]-xcell[0])*(xr[2*r+0]-xcell[0]);
            sep2 += (xr[2*r+1]-xcell[1])*(xr[2*r+1]-xcell[1]);
            if (sep2 < sep2min) {
              eid = e;
              sep2min = sep2;
            }
          }
        }
      }
      eid_list[r] = eid;
      
      
      /* get element -> node map */
      element = c->element;
      elbasis = &element[c->npe*eid];
      
      // find closest //
      sep2min = 1.0e32;
      nid = -1;
      qid = -1;
      for (n=0; n<c->npe; n++) {
        sep2  = (xr[2*r+0]-LA_c[2*elbasis[n]])*(xr[2*r+0]-LA_c[2*elbasis[n]]);
        sep2 += (xr[2*r+1]-LA_c[2*elbasis[n]+1])*(xr[2*r+1]-LA_c[2*elbasis[n]+1]);
        if (sep2 < sep2min) {
          nid = elbasis[n];
          qid = n;
          sep2min = sep2;
        }
      }
      nid_list[r] = nid;
      qid_list[r] = qid;
    }
    
    fp = fopen(filename,"w");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
    fprintf(fp,"# SpecFECtx meta data\n");
    fprintf(fp,"#   mx %d : my %d : basis order %d\n",c->mx_g,c->my_g,c->basisorder);
    fprintf(fp,"# Receiver meta data\n");
    fprintf(fp,"#   + number receiver locations: %d\n",nr);
    fprintf(fp,"#   + takes DR variables from quadrature point nearest to requested receiver location\n");
    for (r=0; r<nr; r++) {
      PetscReal xcell[] = {0,0};
      PetscInt i;
      
      fprintf(fp,"#   + receiver location [%d]: x,y %+1.8e %+1.8e\n",r,xr[2*r+0],xr[2*r+1]);
      fprintf(fp,"#   +   mapped to nearest node --> %+1.8e %+1.8e\n",LA_c[2*nid_list[r]],LA_c[2*nid_list[r]+1]);
      fprintf(fp,"#   +   mapped to nearest element/quad-point --> %d %d\n",eid_list[r],qid_list[r]);
      elbasis = &element[c->npe*eid_list[r]];
      for (i=0; i<c->npe; i++) {
        xcell[0] += LA_c[2*elbasis[i]+0];
        xcell[1] += LA_c[2*elbasis[i]+1];
      }
      xcell[0] = xcell[0] / ((PetscReal)c->npe);
      xcell[1] = xcell[1] / ((PetscReal)c->npe);
      fprintf(fp,"#   +   mapped to nearest element --> %d with centroid %+1.8e %+1.8e\n",eid_list[r],xcell[0],xcell[1]);
    }
    fprintf(fp,"# Time series header <field>(<column index>)\n");
    fprintf(fp,"#   time(1)\n");
    for (r=0; r<nr; r++) {
      PetscInt offset = 1 + r*4; /* 1 is for time */
      
      fprintf(fp,"#     slip(%d) sliprate(%d) mu(%d) sliding(%d) -> station [%d]\n",offset+1,offset+2,offset+3,offset+4,r);
    }
    ierr = VecRestoreArrayRead(coor,&LA_c);CHKERRQ(ierr);
    beenhere = PETSC_TRUE;
  } else {
    fp = fopen(filename,"a");
    if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",filename);
  }
  
  fprintf(fp,"%1.4e",time);
  for (r=0; r<nr; r++) {
    DRVar *cell_data;
    PetscInt e_index,q_index;
    
    e_index = eid_list[r];
    q_index = qid_list[r];
    ierr = SpecFECtxGetDRCellData(c,e_index,&cell_data);CHKERRQ(ierr);
    
    fprintf(fp," %+1.8e %+1.8e %+1.8e %+1.2e",cell_data[q_index].slip,cell_data[q_index].slip_rate,cell_data[q_index].mu,(double)cell_data[q_index].sliding);
  }
  fprintf(fp,"\n");
  
  fclose(fp);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SE2WaveViewer_JSON(SpecFECtx ctx,PetscInt step,PetscReal time,
                                  const char data_description[],PetscInt len,const char *fieldname[],const Vec field[],
                                  const char pbin[],const char jfilename[])
{
  PetscErrorCode ierr;
  char str[PETSC_MAX_PATH_LEN];
  FILE *fp = NULL;
  PetscMPIInt commrank;
  PetscInt k;
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&commrank);CHKERRQ(ierr);
  if (commrank != 0) PetscFunctionReturn(0);
  
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"%s",jfilename);CHKERRQ(ierr);
  fp = fopen(str,"w");
  if (!fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Failed to open file \"%s\"",str);
  
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"{\"se2wave\":{");CHKERRQ(ierr); fprintf(fp,"%s\n",str);
  
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"time\": %1.12e",time);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"step\": %D",step);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"spatial_dimension\": %D",ctx->dim);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"mx\": %D",ctx->mx_g);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"my\": %D",ctx->my_g);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"nx\": %D",ctx->nx_g);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"ny\": %D",ctx->ny_g);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"basis_degree\": %D",ctx->basisorder);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"fields\": [ ");CHKERRQ(ierr); fprintf(fp,"%s",str);
  for (k=0; k<len; k++) {
    ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"\"null\"");CHKERRQ(ierr);
    if (field[k]) {
      ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"\"%s\"",fieldname[k]);CHKERRQ(ierr);
    }
    
    if (k != (len-1)) { fprintf(fp,"%s, ",str);
    } else { fprintf(fp,"%s",str); }
  }
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1," ]");CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  
  //ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"datafile\": \"%s\"",pbin);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"data\":{");CHKERRQ(ierr); fprintf(fp,"%s\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"    \"description\": \"%s\"",data_description);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"    \"fields\": [ ");CHKERRQ(ierr); fprintf(fp,"%s",str);
  for (k=0; k<len; k++) {
    ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"\"null\"");CHKERRQ(ierr);
    if (field[k]) {
      ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"\"%s\"",fieldname[k]);CHKERRQ(ierr);
    }
    
    if (k != (len-1)) { fprintf(fp,"%s, ",str);
    } else { fprintf(fp,"%s",str); }
  }
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1," ]");CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  
  
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"    \"writer\": \"petsc_binary\"");CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"    \"type\": \"Vec\"");CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"    \"filename\": \"%s\"",pbin);CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  /* petsc always writes binary in big endian ordering (even on a small endian machine) */
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"    \"endian\": \"big\"");CHKERRQ(ierr); fprintf(fp,"%s\n",str);
  /*
   #ifdef WORDSIZE_BIGENDIAN
   ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"    \"endian\": \"big\"");CHKERRQ(ierr); fprintf(fp,"%s\n",str);
   #else
   ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"    \"endian\": \"little\"");CHKERRQ(ierr); fprintf(fp,"%s\n",str);
   #endif
   */
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  }");CHKERRQ(ierr); fprintf(fp,"%s,\n",str);
  
  
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"  \"version\": [1,0,0]");CHKERRQ(ierr); fprintf(fp,"%s\n",str);
  ierr = PetscSNPrintf(str,PETSC_MAX_PATH_LEN-1,"}}");CHKERRQ(ierr); fprintf(fp,"%s\n",str);
  
  fclose(fp);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SE2WaveCoordinateViewerViewer(SpecFECtx ctx,PetscInt step,PetscReal time,const char prefix[])
{
  PetscErrorCode ierr;
  PetscViewer vu;
  Vec coor = NULL;
  char fname[PETSC_MAX_PATH_LEN];
  
  ierr = DMGetCoordinates(ctx->dm,&coor);CHKERRQ(ierr);
  if (!coor) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Must have a valid coordinate vector");
  
  ierr = PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"%s_coor.pbin",prefix);CHKERRQ(ierr);
  
  {
    char jname[PETSC_MAX_PATH_LEN];
    Vec input[] = {NULL};
    const char *fieldname[] = { "coor" };
    
    input[0] = coor;
    ierr = PetscSNPrintf(jname,PETSC_MAX_PATH_LEN-1,"%s_coor.json",prefix);CHKERRQ(ierr);
    ierr = SE2WaveViewer_JSON(ctx,step,time,"coordinates",1,fieldname,input,fname,jname);CHKERRQ(ierr);
  }
  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fname,FILE_MODE_WRITE,&vu);CHKERRQ(ierr);
  
  ierr = PetscViewerBinaryWrite(vu,(void*)&ctx->mx_g,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscViewerBinaryWrite(vu,(void*)&ctx->my_g,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscViewerBinaryWrite(vu,(void*)&ctx->nx_g,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscViewerBinaryWrite(vu,(void*)&ctx->ny_g,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  
  ierr = VecView(coor,vu);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&vu);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode SE2WaveWaveFieldViewer(SpecFECtx ctx,PetscInt step,PetscReal time,Vec u,Vec v,const char prefix[])
{
  PetscErrorCode ierr;
  PetscViewer vu;
  char fname[PETSC_MAX_PATH_LEN];
  
  if (!u && !v) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"At least one of the displacement or velocity vectors must be non-NULL");
  
  ierr = PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"%s_wavefield.pbin",prefix);CHKERRQ(ierr);
  
  {
    char jname[PETSC_MAX_PATH_LEN];
    Vec input[] = {NULL,NULL};
    const char *fieldname[] = { "u", "v" };
    input[0] = u;
    input[1] = v;
    
    ierr = PetscSNPrintf(jname,PETSC_MAX_PATH_LEN-1,"%s_wavefield.json",prefix);CHKERRQ(ierr);
    ierr = SE2WaveViewer_JSON(ctx,step,time,"wavefield",2,fieldname,input,fname,jname);CHKERRQ(ierr);
  }
  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fname,FILE_MODE_WRITE,&vu);CHKERRQ(ierr);
  
  ierr = PetscViewerBinaryWrite(vu,(void*)&ctx->mx_g,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscViewerBinaryWrite(vu,(void*)&ctx->my_g,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscViewerBinaryWrite(vu,(void*)&ctx->nx_g,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscViewerBinaryWrite(vu,(void*)&ctx->ny_g,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  
  ierr = PetscViewerBinaryWrite(vu,(void*)&step,1,PETSC_INT,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscViewerBinaryWrite(vu,(void*)&time,1,PETSC_REAL,PETSC_FALSE);CHKERRQ(ierr);
  
  if (u) { ierr = VecView(u,vu);CHKERRQ(ierr); }
  if (v) { ierr = VecView(v,vu);CHKERRQ(ierr); }
  ierr = PetscViewerDestroy(&vu);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode se2dr_demo(PetscInt mx,PetscInt my)
{
  PetscErrorCode ierr;
  SpecFECtx ctx;
  PetscInt p,k,nt,of;
  PetscViewer viewer;
  Vec u,v,a,f,g,Md;
  PetscReal time,dt,time_max;
  PetscInt nrecv;
  PetscReal *xr_list;
  PetscBool dump_ic_src_vts = PETSC_FALSE;
  PetscBool ignore_receiver_output = PETSC_FALSE;
  PetscReal nrm,max,min,dx,dy;
  char vts_fname[PETSC_MAX_PATH_LEN];
  
  
  ierr = PetscOptionsGetBool(NULL,NULL,"-dump_ic_src",&dump_ic_src_vts,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-ignore_receiver_output",&ignore_receiver_output,NULL);CHKERRQ(ierr);
  
  /*
    Create the structured mesh for the spectral element method.
   The default mesh is defined over the domain [0,1]^2.
  */
  ierr = SpecFECtxCreate(&ctx);CHKERRQ(ierr);
  p = 2;
  ierr = PetscOptionsGetInt(NULL,NULL,"-bdegree",&p,NULL);CHKERRQ(ierr);
  ierr = SpecFECtxCreateMesh(ctx,2,mx,my,PETSC_DECIDE,p,2);CHKERRQ(ierr);

  /*
   Define your domain by shifting and scaling the default [0,1]^2 domain
  */
  {
    PetscReal alpha = 20.0e3;
    PetscReal scale[] = {  3.0*alpha, 3.0*alpha };
    PetscReal shift[] = { -1.5*alpha,-1.5*alpha };
    
    ierr = SpecFECtxScaleMeshCoords(ctx,scale,shift);CHKERRQ(ierr);
  }

  /* 
   Specify fault dimensions
  */
  {
    PetscReal gmin[3],gmax[3];
    
    ierr = DMGetBoundingBox(ctx->dm,gmin,gmax);CHKERRQ(ierr);
    dx = (gmax[0] - gmin[0])/((PetscReal)ctx->mx_g);
    dy = (gmax[1] - gmin[1])/((PetscReal)ctx->my_g);
  }
  PetscPrintf(PETSC_COMM_WORLD,"[se2dr] cell sizes: dx = %1.4e, dy = %1.4e\n",dx,dy);
  
  
  ctx->delta = 25.0;
  ierr = PetscOptionsGetReal(NULL,NULL,"-delta",&ctx->delta,NULL);CHKERRQ(ierr);
  {
    PetscBool found;
    PetscReal delta_factor = 0.0;
    
    found = PETSC_FALSE;
    ierr = PetscOptionsGetReal(NULL,NULL,"-delta_cell_factor",&delta_factor,&found);CHKERRQ(ierr);
    if (found) {
      ctx->delta = dy * delta_factor;
    }
  }
  PetscPrintf(PETSC_COMM_WORLD,"[se2dr] using fault delta = %1.4e\n",ctx->delta);
  PetscPrintf(PETSC_COMM_WORLD,"[se2dr] elements across fault = %1.4e\n",2.0 * ctx->delta/dy);
  
  //ierr = FaultSDFInit_v1(ctx);CHKERRQ(ierr);
  ierr = FaultSDFInit_v2(ctx);CHKERRQ(ierr);
  
  /*
   Specify the material properties for the domain.
   This function sets constant material properties in every cell.
   More general methods can be easily added.
  */
  ierr = SpecFECtxSetConstantMaterialProperties_Velocity(ctx,6000.0 ,3464.0, 2670.0);CHKERRQ(ierr); // vp,vs,rho
  /* Linear slip weakening parameters */
  ctx->mu_s = 0.677;
  ctx->mu_d = 0.525;
  ctx->D_c = 0.40;

  
  ierr = DMCreateGlobalVector(ctx->dm,&u);CHKERRQ(ierr); ierr = PetscObjectSetName((PetscObject)u,"disp");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->dm,&v);CHKERRQ(ierr); ierr = PetscObjectSetName((PetscObject)v,"velo");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->dm,&a);CHKERRQ(ierr); ierr = PetscObjectSetName((PetscObject)a,"accl");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->dm,&f);CHKERRQ(ierr); ierr = PetscObjectSetName((PetscObject)f,"f");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->dm,&g);CHKERRQ(ierr); ierr = PetscObjectSetName((PetscObject)g,"g");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->dm,&Md);CHKERRQ(ierr);
  
  ierr = VecZeroEntries(u);CHKERRQ(ierr);
  
  /*
   Write out the mesh and intial values for the displacement, velocity and acceleration (u,v,a)
  */
  if (dump_ic_src_vts) {
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,"uva.vts",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(u,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
  ierr = AssembleBilinearForm_Mass2d(ctx,Md);CHKERRQ(ierr);
  
  /*
   Define the location of the receivers
  */
  nrecv = 8;
  ierr = PetscMalloc1(nrecv*2,&xr_list);CHKERRQ(ierr);
  xr_list[0] = 4.0e3;
  xr_list[1] = ctx->delta;

  xr_list[2] = 4.0e3;
  xr_list[3] = -ctx->delta;

  xr_list[4] = 6.0e3;
  xr_list[5] = ctx->delta;

  xr_list[6] = 6.0e3;
  xr_list[7] = -ctx->delta;

  xr_list[8] = 8.0e3;
  xr_list[9] = ctx->delta;

  xr_list[10] = 8.0e3;
  xr_list[11] = -ctx->delta;

  xr_list[12] = 10.0e3;
  xr_list[13] = ctx->delta;

  xr_list[14] = 10.0e3;
  xr_list[15] = -ctx->delta;


  /* Initialize time loop */
  k = 0;
  time = 0.0;
  
  time_max = 0.4;
  ierr = PetscOptionsGetReal(NULL,NULL,"-tmax",&time_max,NULL);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"[se2dr] Requested time period: %1.4e\n",time_max);
  
  ierr = ElastoDynamicsComputeTimeStep_2d(ctx,&dt);CHKERRQ(ierr);
  dt = dt * 0.5;
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"[se2dr] Using time step size: %1.4e\n",dt);
  
  nt = 1000000;
  nt = (PetscInt)(time_max / dt ) + 4;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nt",&nt,NULL);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"[se2dr] Estimated number of time steps: %D\n",nt);
  
  of = 5000;
  ierr = PetscOptionsGetInt(NULL,NULL,"-of",&of,NULL);CHKERRQ(ierr);
  
  ierr = SE2WaveCoordinateViewerViewer(ctx,0,0.0,"default_mesh");CHKERRQ(ierr);
  {
    char prefix[PETSC_MAX_PATH_LEN];
    
    ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step-%.4D",0);CHKERRQ(ierr);
    ierr = SE2WaveWaveFieldViewer(ctx,k,time,u,v,prefix);CHKERRQ(ierr);
  }
  
  if (k%of == 0) {
    ierr = PetscSNPrintf(vts_fname,PETSC_MAX_PATH_LEN-1,"step-%.4D.vts",0);CHKERRQ(ierr);
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,vts_fname,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = VecView(u,viewer);CHKERRQ(ierr);
    ierr = VecView(v,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
  /* Perform time stepping */
  for (k=1; k<=nt; k++) {
    
    time = time + dt;
    
    ierr = VecAXPY(u,dt,v);CHKERRQ(ierr); /* u_{n+1} = u_{n} + dt.v_{n} */
    
    ierr = VecAXPY(u,0.5*dt*dt,a);CHKERRQ(ierr); /* u_{n+1} = u_{n+1} + 0.5.dt^2.a_{n} */
    
    ierr = VecAXPY(v,0.5*dt,a);CHKERRQ(ierr); /* v' = v_{n} + 0.5.dt.a_{n} */
    
    /* Compute f = -F^{int}( u_{n+1} ) */
    //printf("[time:%f] ",time);
    ierr = AssembleLinearForm_ElastoDynamics_StressGlut2d_tpv(ctx,u,v,dt,time,0.5*dt,f);CHKERRQ(ierr);
    //ierr = AssembleLinearForm_ElastoDynamics_StressGlut2d_tpv_cellwise(ctx,u,v,dt,time,1.1*dt,f);CHKERRQ(ierr);
    //ierr = AssembleLinearForm_ElastoDynamics_StressGlut2d_tpv_cellwise_v2(ctx,u,v,dt,time,1.1*dt,f);CHKERRQ(ierr);
    
    /* Update force; F^{ext}_{n+1} = f + S(t_{n+1}) g(x) */
    ierr = VecAXPY(f,1.0,g);CHKERRQ(ierr);
    
    /* "Solve"; a_{n+1} = M^{-1} f */
    ierr = VecPointwiseDivide(a,f,Md);CHKERRQ(ierr);
    
    /* Update velocity */
    ierr = VecAXPY(v,0.5*dt,a);CHKERRQ(ierr); /* v_{n+1} = v' + 0.5.dt.a_{n+1} */
    
    /* Update slip-rate & slip */
    ierr = Update_StressGlut2d(ctx,u,v,dt);CHKERRQ(ierr);

    
    if (k%10 == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"[step %9D] time = %1.4e : dt = %1.4e \n",k,time,dt);
      VecNorm(u,NORM_2,&nrm);
      VecMin(u,0,&min);
      VecMax(u,0,&max); PetscPrintf(PETSC_COMM_WORLD,"  [displacement] max = %+1.4e : min = %+1.4e : l2 = %+1.4e \n",max,min,nrm);
      VecNorm(v,NORM_2,&nrm);
      VecMin(v,0,&min);
      VecMax(v,0,&max); PetscPrintf(PETSC_COMM_WORLD,"  [velocity]     max = %+1.4e : min = %+1.4e : l2 = %+1.4e \n",max,min,nrm);
    }

    /*
      Write out the u,v,a values at each receiver
    */
    if (!ignore_receiver_output) {
      ierr = RecordUVA_MultipleStations_NearestGLL(ctx,time,nrecv,xr_list,u,v,a);CHKERRQ(ierr);
      ierr = RecordDRVar_MultipleStations_NearestGLL_SEQ(ctx,time,nrecv,xr_list);CHKERRQ(ierr);
    }
    
    if (k%of == 0) {
      ierr = PetscSNPrintf(vts_fname,PETSC_MAX_PATH_LEN-1,"step-%.4D.vts",k);CHKERRQ(ierr);
      ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,vts_fname,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(u,viewer);CHKERRQ(ierr);
      ierr = VecView(v,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
    if (k%of == 0) {
      char prefix[PETSC_MAX_PATH_LEN];
      ierr = PetscSNPrintf(prefix,PETSC_MAX_PATH_LEN-1,"step-%.4D",k);CHKERRQ(ierr);
      ierr = SE2WaveWaveFieldViewer(ctx,k,time,u,v,prefix);CHKERRQ(ierr);
    }
    
    if (time >= time_max) {
      break;
    }
  }
  PetscPrintf(PETSC_COMM_WORLD,"[step %9D] time = %1.4e : dt = %1.4e \n",k,time,dt);
  VecNorm(u,NORM_2,&nrm);
  VecMin(u,0,&min);
  VecMax(u,0,&max); PetscPrintf(PETSC_COMM_WORLD,"  [displacement] max = %+1.4e : min = %+1.4e : l2 = %+1.4e \n",max,min,nrm);
  VecNorm(v,NORM_2,&nrm);
  VecMin(v,0,&min);
  VecMax(v,0,&max); PetscPrintf(PETSC_COMM_WORLD,"  [velocity]     max = %+1.4e : min = %+1.4e : l2 = %+1.4e \n",max,min,nrm);
  
  /* plot last snapshot */
  ierr = PetscSNPrintf(vts_fname,PETSC_MAX_PATH_LEN-1,"step-%.4D.vts",k);CHKERRQ(ierr);
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,vts_fname,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(u,viewer);CHKERRQ(ierr);
  ierr = VecView(v,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  
  ierr = PetscFree(xr_list);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&v);CHKERRQ(ierr);
  ierr = VecDestroy(&a);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr);
  ierr = VecDestroy(&Md);CHKERRQ(ierr);
  ierr = VecDestroy(&g);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  PetscInt       mx,my;
  PetscMPIInt    size;
  
  ierr = PetscInitialize(&argc,&args,(char*)0,NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  
  mx = my = 8;
  ierr = PetscOptionsGetInt(NULL,NULL,"-mx",&mx,NULL);CHKERRQ(ierr);
  my = mx;
  ierr = PetscOptionsGetInt(NULL,NULL,"-my",&my,NULL);CHKERRQ(ierr);

  ierr = se2dr_demo(mx,my);CHKERRQ(ierr);
  
  ierr = PetscFinalize();
  return(ierr);
}
