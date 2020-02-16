#ifndef __FE__H
#define __FE__H

#include<petscdt.h>
#include<petscmath.h>

typedef enum {GAUSS,LGL} QuadMode; // LGL: Legendre-Gauss-Lobatto quadrature
                                   // GAUSS: Gauss quadrature (Also known as Gauss-Legendre)

typedef struct FE_private *FE;

struct FE_private
{
  PetscInt P;            // Finite Element Polynomial degree + 1 (number of nodes) in 1D
  PetscInt Q;            // Number of quadrature points in 1D
  PetscReal *qref1d;     // Quadrature points in 1D  in reference element
  PetscReal *qweight1d;  // weights of quadrature points in 1D  in reference element
  PetscReal *B1d;        // Evaluated Basis functions at quadrature points in 1D (reference element)
  PetscReal *D1d;        // Evaluated Derivative of basis functions at quadrature points in 1D (reference element)
};

PetscErrorCode FEcreate(PetscInt P, PetscInt Q, FE *fe);
PetscErrorCode FEsetup(FE fe, QuadMode qmod);
PetscErrorCode FEdestroy(FE *fe);
PetscErrorCode LobattoQuadrature(PetscInt Q, PetscScalar *qref1d, PetscScalar *qweight1d); //From CEED 2018 and Forenberg 1998
PetscErrorCode FEBasisEval(PetscInt P, PetscInt Q, PetscScalar nodes[], PetscScalar qref1d[], PetscScalar *interp1d, PetscScalar *grad1d);
PetscErrorCode FEView(FE fe);
PetscErrorCode FEgetComponents(FE fe, const PetscScalar **W1d, const PetscScalar **B1d, const PetscScalar **D1d);

#endif //end of __FE__H
