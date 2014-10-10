#include "eltdef.h"

int bemDimension(const double* const Elt, const unsigned int& nElt,
                     unsigned int* const  TypeID, char* TypeName[], char* TypeKeyOpts[],
                     unsigned int* const nKeyOpt, const unsigned int& nEltType)
{
 /*  Return boundary element dimensions. The dimension equals 2 for a 2D
  *  problem, 3 for a 3D problem. The function also checks wether
  *  all elements have the same dimension.
  */
  int PreviousEltDim=-1;
  for (unsigned int iElt=0; iElt<nElt; iElt++)
  {
    const unsigned int EltType = (unsigned int)(Elt[nElt+iElt]);
    unsigned int Parent;
    unsigned int nNod;
    unsigned int nCol;
    unsigned int ShapeTypeN;
    unsigned int ShapeTypeM;
    unsigned int EltDim;
    unsigned int AxiSym;
    unsigned int Periodic;
    unsigned int nGauss;
    unsigned int nEltDiv;
    unsigned int nGaussSing;
    unsigned int nEltDivSing;
    eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,Parent, 
           nNod,nCol,ShapeTypeN,ShapeTypeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing);
    
    if ((!(PreviousEltDim==-1)) && (!(PreviousEltDim==(int)EltDim)))
       throw("Boundary element mesh contains elements with incompatible dimensions."); 
    PreviousEltDim=EltDim;
  }
  // Problem dimension: EltDim+1
  return(PreviousEltDim+1);
}
