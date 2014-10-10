#include "eltdef.h"

bool isAxisym(const double* const Elt, const unsigned int& nElt,
              unsigned int* const  TypeID, char* TypeName[], char* TypeKeyOpts[],
              unsigned int* const nKeyOpt, const unsigned int& nEltType)
/*
 * Check mesh geometry: returns true if geometry is axisymmetric.
 */
{
   int PreviousEltGeom=-1;
   for (unsigned int iElt=0; iElt<nElt; iElt++)
   {
     const unsigned int EltType = (unsigned int)(Elt[nElt+iElt]);
     unsigned int Parent;
     unsigned int nEltNod;
     unsigned int nEltColl;
     unsigned int ShapeTypeN;
     unsigned int ShapeTypeM;
     unsigned int EltDim;
     unsigned int AxiSym;
     unsigned int Periodic;
     unsigned int nGauss;
     unsigned int nEltDiv;
     unsigned int nGaussSing;
     unsigned int nEltDivSing;
     eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,Parent,nEltNod,
           nEltColl,ShapeTypeN,ShapeTypeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing);
           
     if ((!(PreviousEltGeom==-1)) && (!(PreviousEltGeom==(int)AxiSym)))
        throw("Boundary element mesh contains elements with axisymmetric and non-axisymmetric geometry.");
     PreviousEltGeom=AxiSym;
   }
     // Problem dimension: EltDim+1
     bool isAxisym = (PreviousEltGeom == 1);
     return(isAxisym);
}
