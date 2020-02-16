#include "fe.h"

PetscErrorCode FEcreate(PetscInt P, PetscInt Q, FE *fe){

  PetscErrorCode ierr;

  FE feSpace;
  PetscFunctionBeginUser;

  ierr = PetscNew(&feSpace);CHKERRQ(ierr);
  feSpace->P = P;
  feSpace->Q = Q;
  ierr = PetscCalloc4(P*Q,&feSpace->B1d, P*Q,&feSpace->D1d, Q,&feSpace->qref1d, Q,&feSpace->qweight1d);CHKERRQ(ierr);
  *fe = feSpace;
  PetscFunctionReturn(0);
}

PetscErrorCode FEdestroy(FE *fe){

  PetscErrorCode ierr;

  FE feSpace;

  PetscFunctionBeginUser;
  if(!fe) PetscFunctionReturn(0);
  feSpace = *fe;
  ierr = PetscFree(feSpace->B1d);CHKERRQ(ierr);
  ierr = PetscFree(feSpace->D1d);CHKERRQ(ierr);
  ierr = PetscFree(feSpace->qref1d);CHKERRQ(ierr);
  ierr = PetscFree(feSpace->qweight1d);CHKERRQ(ierr);
  ierr = PetscFree(feSpace);CHKERRQ(ierr);
  *fe = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode FEsetup(FE fe, QuadMode qmode){

  PetscErrorCode ierr;
  PetscInt P, Q;
  PetscScalar *nodes;

  PetscFunctionBeginUser;
  P = fe->P;
  Q = fe->Q;
  ierr = PetscCalloc1(P, &nodes);
  ierr = LobattoQuadrature(P, nodes, NULL); CHKERRQ(ierr);
  switch (qmode) {
      case GAUSS:
        ierr = PetscDTGaussQuadrature(Q,-1,1,fe->qref1d, fe->qweight1d);CHKERRQ(ierr);
        break;
      case LGL:
        ierr = LobattoQuadrature(Q, fe->qref1d, fe->qweight1d);CHKERRQ(ierr);
        break;
    }
    ierr = FEBasisEval(P, Q, nodes, fe->qref1d, fe->B1d, fe->D1d);CHKERRQ(ierr);

  PetscFree(nodes);
  PetscFunctionReturn(0);
}

PetscErrorCode LobattoQuadrature(PetscInt Q, PetscScalar *qref1d, PetscScalar *qweight1d){

  PetscScalar P0, P1, P2, dP2, d2P2, xi, wi;

  PetscFunctionBeginUser;

  // Build qref1d, qweight1d
  // Set endpoints
  wi = 2.0/((PetscScalar)(Q*(Q-1)));
  if (qweight1d) {
    qweight1d[0] = wi;
    qweight1d[Q-1] = wi;
  }
  qref1d[0] = -1.0;
  qref1d[Q-1] = 1.0;
  // Interior
  for (PetscInt i = 1; i <= (Q-1)/2; i++) {
    // Guess
    xi = PetscCosReal(PETSC_PI*(PetscScalar)(i)/(PetscScalar)(Q-1));

    // Pn(xi)
    P0 = 1.0;
    P1 = xi;
    P2 = 0.0;
    for (PetscInt j = 2; j < Q; j++) {
      P2 = (((PetscScalar)(2*j-1))*xi*P1-((PetscScalar)(j-1))*P0)/((PetscScalar)(j));
      P0 = P1;
      P1 = P2;
    }
    // First Newton step
    dP2 = (xi*P2 - P0)*(PetscScalar)Q/(xi*xi-1.0);
    d2P2 = (2*xi*dP2 - (PetscScalar)(Q*(Q-1))*P2)/(1.0-xi*xi);
    xi = xi-dP2/d2P2;
    // Newton to convergence
    for (PetscInt k=0; k<100 && PetscAbsReal(dP2)>1e-15; k++) {
      P0 = 1.0;
      P1 = xi;
      for (PetscInt j = 2; j < Q; j++) {
        P2 = (((PetscScalar)(2*j-1))*xi*P1-((PetscScalar)(j-1))*P0)/((PetscScalar)(j));
        P0 = P1;
        P1 = P2;
      }
      dP2 = (xi*P2 - P0)*(PetscScalar)Q/(xi*xi-1.0);
      d2P2 = (2*xi*dP2 - (PetscScalar)(Q*(Q-1))*P2)/(1.0-xi*xi);
      xi = xi-dP2/d2P2;
    }
    // Save xi, wi
    wi = 2.0/(((PetscScalar)(Q*(Q-1)))*P2*P2);
    if (qweight1d) {
      qweight1d[i] = wi;
      qweight1d[Q-1-i] = wi;
    }
    qref1d[i] = -xi;
    qref1d[Q-1-i]= xi;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FEBasisEval(PetscInt P, PetscInt Q, PetscScalar nodes[], PetscScalar qref1d[],PetscScalar *B1d, PetscScalar *D1d){

  PetscInt i,j,k;
  PetscScalar c1, c2, c3, c4, dx;

  PetscFunctionBeginUser;

  // Build B, D matrix
  // Fornberg, 1998
  for (i = 0; i  < Q; i++) {
    c1 = 1.0;
    c3 = nodes[0] - qref1d[i];
    B1d[i*P+0] = 1.0;
    for (j = 1; j < P; j++) {
      c2 = 1.0;
      c4 = c3;
      c3 = nodes[j] - qref1d[i];
      for (k = 0; k < j; k++) {
        dx = nodes[j] - nodes[k];
        c2 *= dx;
        if (k == j - 1) {
          D1d[i*P + j] = c1*(B1d[i*P + k] - c4*D1d[i*P + k]) / c2;
          B1d[i*P + j] = - c1*c4*B1d[i*P + k] / c2;
        }
        D1d[i*P + k] = (c3*D1d[i*P + k] - B1d[i*P + k]) / dx;
        B1d[i*P + k] = c3*B1d[i*P + k] / dx;
      }
      c1 = c2;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FEgetComponents(FE fe, const PetscScalar **W1d, const PetscScalar **B1d, const PetscScalar **D1d){

PetscFunctionBeginUser;
if (W1d) *W1d = fe->qweight1d;
if (B1d) *B1d = fe->B1d;
if (D1d) *D1d = fe->D1d;
PetscFunctionReturn(0);
}

PetscErrorCode FEView(FE fe){
    
  PetscInt         i;
  PetscErrorCode   ierr;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_SELF,"\nFE Space:\n");CHKERRQ(ierr);
  //print P and Q
  ierr = PetscPrintf(PETSC_COMM_SELF,"  fe->P: %d\n", fe->P);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"  fe->Q: %d\n", fe->Q);CHKERRQ(ierr);
  //print quadrature points
  ierr = PetscPrintf(PETSC_COMM_SELF,"  fe->qref1d: [ ");CHKERRQ(ierr);
  for(i=0; i<fe->Q;i++)
    ierr = PetscPrintf(PETSC_COMM_SELF,"%g ", fe->qref1d[i]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"]\n");CHKERRQ(ierr);
  //print weights associated with quadrature points
  ierr = PetscPrintf(PETSC_COMM_SELF,"  fe->qweight1d: [ ");CHKERRQ(ierr);
  for(i=0; i<fe->Q;i++)
    ierr = PetscPrintf(PETSC_COMM_SELF,"%g ", fe->qweight1d[i]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"]\n");CHKERRQ(ierr);
  //print B matrix in array form
  ierr = PetscPrintf(PETSC_COMM_SELF,"  fe->B1d: [ ");CHKERRQ(ierr);
  for(i=0; i<(fe->P)*(fe->Q);i++)
    ierr = PetscPrintf(PETSC_COMM_SELF,"%g ", fe->B1d[i]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"]\n");CHKERRQ(ierr);
  //print D matrix in array form
  ierr = PetscPrintf(PETSC_COMM_SELF,"  fe->D1d: [ ");CHKERRQ(ierr);
  for(i=0; i<(fe->P)*(fe->Q);i++)
    ierr = PetscPrintf(PETSC_COMM_SELF,"%g ", fe->D1d[i]);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"]\n");CHKERRQ(ierr);

 PetscFunctionReturn(0);
}
