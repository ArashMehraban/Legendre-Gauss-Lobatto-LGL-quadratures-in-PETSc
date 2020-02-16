### Legendre-Gauss-Lobatto (LGL) quadratures in PETSc
This code provides the Legendre-Gauss-Lobatto (LGL) quadratures in PETSc in 1D. This is the PETSc version of LGL found in libCEED [libCEED](https://github.com/CEED/libCEED).

### How to compile:
1) set ``PETSC_DIR`` and ``PETSC_ARCH``
2) make all

### How to run:
./main -degree [n]

Example:\
``./main -degree 4``

``degree`` corresponds to the desired polynomial degree

```
Guass:
FE Space:
  fe->P: 4
  fe->Q: 4
  fe->qref1d: [ -0.861136 -0.339981 0.339981 0.861136 ]
  fe->qweight1d: [ 0.347855 0.652145 0.652145 0.347855 ]
  fe->B1d: [ 0.629943 0.472559 -0.149503 0.0470015 -0.0706948 0.972976 0.13254 -0.0348213 -0.0348213 0.13254 0.972976 -0.0706948 0.0470015 -0.149503 0.472559 0.629943 ]
  fe->D1d: [ -2.34184 2.78794 -0.635104 0.188997 -0.516702 -0.487952 1.33791 -0.33325 0.33325 -1.33791 0.487952 0.516702 -0.188997 0.635104 -2.78794 2.34184 ]

LGL:
FE Space:
  fe->P: 4
  fe->Q: 4
  fe->qref1d: [ -1. -0.447214 0.447214 1. ]
  fe->qweight1d: [ 0.166667 0.833333 0.833333 0.166667 ]
  fe->B1d: [ 1. -0. 0. -0. 0. 1. -0. 0. -0. 0. 1. -0. 0. -0. 0. 1. ]
  fe->D1d: [ -3. 4.04508 -1.54508 0.5 -0.809017 -1.53429e-16 1.11803 -0.309017 0.309017 -1.11803 0. 0.809017 -0.5 1.54508 -4.04508 3. ]
```

``fe->qref1d``    : Quadrature Points\
``fe->qweight1d`` : Quadrature Weights\
``fe->B1d``       : Evaluation of Shape (Basis) functions at Quadrature Points\
``fe->D1d``       : Evaluation of Derivative of Shape (Basis) functions at Quadrature Points
