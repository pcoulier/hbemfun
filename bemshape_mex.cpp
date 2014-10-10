/*BEMSHAPE   Boundary element shape function.
 *
 *   N = BEMSHAPE(ShapeID,xi) computes the boundary element shape function
 *   of the type ShapeID in the points with natural coordinates xi.
 *
 *   ShapeType Shape function type:
 *             1: Constant element shape function (nNod=1)
 *             2: 3D 3-node triangular shape function (nNod=3)
 *             3: 3D 6-node triangular shape function (nNod=6)
 *             4: 3D 4-node quadrilateral shape function (nNod=4)
 *             5: 3D 8-node quadrilateral shape function (nNod=8)
 *             6: 3D 9-node quadrilateral shape function (nNod=9)
 *             7: 2D 2-node linear shape function (nNod=2)
 *             8: 2D 3-node quadratic shape function (nNod=3)
 *             9: 2D 4-node cubic shape function (nNod=4)
 *   xi        Natural coordinates (nXi * nDim). nDim equals 1 for 2D elements
 *             and 2 for 3D elements.
 *   N         Boundary element shape function (nNod * nXi).
 */

/* $Make: mex -O -output bemshape bemshape_mex.cpp shapefun.cpp $ */

#include "mex.h"
#include <string.h>
#include "shapefun.h"
#include "checklicense.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
     checklicense();
    
     if (nrhs<2) throw("Not enough input arguments.");
     if (nrhs>2) throw("Too many input arguments.");
     if (nlhs>1) throw("Too many output arguments.");

     if (!mxIsNumeric(prhs[0])) throw("Input argument 'ShapeType' must be numeric.");
     if (mxIsSparse(prhs[0])) throw("Input argument 'ShapeType' must not be sparse.");
     if (mxIsComplex(prhs[0])) throw("Input argument 'ShapeType' must be real.");
     if (!(mxGetNumberOfElements(prhs[0])==1)) throw("Input argument 'ShapeType' should be a scalar.");
     const double* const Inarg=mxGetPr(prhs[0]);
     const int shapeType = int(Inarg[0]);

     int nRow;
     if (shapeType==1) nRow=1;
     if (shapeType==2) nRow=3;
     if (shapeType==3) nRow=6;
     if (shapeType==4) nRow=4;
     if (shapeType==5) nRow=8;
     if (shapeType==6) nRow=9;
     if (shapeType==7) nRow=2;
     if (shapeType==8) nRow=3;
     if (shapeType==9) nRow=4;

     int xiDim;
     if (shapeType==1) xiDim=0;
     if (shapeType>1)  xiDim=2;
     if (shapeType>6)  xiDim=1;

     if (!mxIsNumeric(prhs[1])) throw("Input argument 'Xi' must be numeric.");
     if (mxIsSparse(prhs[1])) throw("Input argument 'Xi' must not be sparse.");
     if (mxIsComplex(prhs[1])) throw("Input argument 'Xi' must be real.");
     const int nXi=mxGetM(prhs[1]);
     if ((xiDim>0)&&(!(mxGetN(prhs[1])==(unsigned)xiDim)))  throw("Wrong number of columns for input arument 'Xi'");
     const double* const xi=mxGetPr(prhs[1]);

     plhs[0]=mxCreateDoubleMatrix(nRow,nXi,mxREAL);
     double* const N=mxGetPr(plhs[0]);
     shapefun(shapeType,nXi,xi,N);
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
