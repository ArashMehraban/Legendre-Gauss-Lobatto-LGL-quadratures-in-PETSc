static char help[] = "Gauss-Legendre-Lobatto quadrature";

#include "fe.h"
int main(int argc, char **argv)
{
  PetscErrorCode    ierr;
  PetscBool         degreeFalg = PETSC_FALSE;
  PetscInt          degree;
  PetscInt          P, Q;
  FE                feGauss,feLGL;
  const PetscScalar *W, *B,*D;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "LGL Quadratures in PETSc", NULL); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-degree", "Polynomial degree of tensor product basis", NULL, degree, &degree, &degreeFalg); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  if(!degreeFalg) {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "-degree option needed\n\n"); CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "degree ERROR!");
    }

  P = degree + 1;
  Q = P;

  // Create and view LGL space
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Guass:\n"); CHKERRQ(ierr);
  ierr = FEcreate(P,Q, &feGauss);
  ierr = FEsetup(feGauss, GAUSS);
  ierr = FEgetComponents(feGauss, &W, &B, &D);
  ierr = FEView(feGauss);

  // Create and view LGL space
  ierr = PetscPrintf(PETSC_COMM_WORLD, "LGL:\n"); CHKERRQ(ierr);
  ierr = FEcreate(P,Q, &feLGL);
  ierr = FEsetup(feLGL, LGL);
  ierr = FEgetComponents(feLGL, &W, &B, &D);
  ierr = FEView(feLGL);

  //Destroy the FE space
  ierr = FEdestroy(&feGauss);
  ierr = FEdestroy(&feLGL);

  ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}//end of main
