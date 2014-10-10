#include "eltdef.h"
#include "mex.h"
#include "shapefun.h"
#include <complex>

using namespace std;

void BemNodeIndex(const double* const Nod, const unsigned int& nNod,
                  const unsigned int& NodeID, int& index)
/* BemNodeIndex Lookup node index from NodeID.
 *
 * Nod        Node Array.
 * nNod       Number of nodes.
 * NodeID     NodeId for which index is requested
 * index      resulting node index : Nod[index]=NodeID.
 */
{
  bool found=false;
  index=-1;
//   mexPrintf("nNod: %d\n",nNod); // DEBUG
  while ((!found) && (index<((int)nNod)))
  {
    index++;
    if ((unsigned int)(Nod[index])==NodeID) found=true;
//     mexPrintf("index: %d\n",index); // DEBUG
  }
  if (index==(int)nNod) throw("Unknown node in element array.");
}


void BemCollPoints(const double* const Elt, const double* const Nod,
                   unsigned int* const  TypeID, unsigned int* const  nKeyOpt,
                   const char* const TypeName[],const char* const TypeKeyOpts[],
                   const unsigned int& nEltType, const unsigned int& nElt, const unsigned int& maxEltCol, 
                   const unsigned int& nNod, unsigned int* const NodalColl,
                   unsigned int* const CentroidColl, unsigned int& nNodalColl, 
                   unsigned int& nCentroidColl)
/* BemCollPoints: Retreive collocation points.
 * There are 2 types of collocation points: Nodal collocation
 * (type 1 collocation) and centroid collocation (type 2 collocation).
 *
 *   Elt                  Element array.
 *   Nod                  Node Array.
 *   nElt                 Number of elements.
 *   maxEltCol            Number of Columns in element array.
 *   nNod                 Number of nodes.
 *   NodalColl [nNod]     Vector indicating(=1) if node is
 *                        a collocation point.
 *   CentroidColl [nElt]  Vector indicating(=1) if the element
 *                        centroid is a collocation point.
 *   nNodalColl           Total number of nodal collocation points.
 *   nCentroidColl        Total number of centroid collocation points.
 */
{
  for (unsigned int iNod=0; iNod<nNod; iNod++) NodalColl[iNod]=0;
  for (unsigned int iElt=0; iElt<nElt; iElt++) CentroidColl[iElt]=0;

  unsigned int EltType;
  unsigned int EltParent;
  unsigned int nEltNod;
  unsigned int nEltColl;
  unsigned int EltShapeN;
  unsigned int EltShapeM;
  unsigned int EltDim;
  unsigned int AxiSym;
  unsigned int Periodic;
  unsigned int nGauss;
  unsigned int nEltDiv;
  unsigned int nGaussSing;
  unsigned int nEltDivSing;

  for (unsigned int iElt=0; iElt<nElt; iElt++)
  {
    EltType = (unsigned int)(Elt[nElt+iElt]);
    eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
           nEltType,EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,
           nGauss,nEltDiv,nGaussSing,nEltDivSing);
    
    if (maxEltCol<2+nEltNod) throw("Number of colums in the element array is incompatible with elements defined.");

    // Check for centroid collocation.
    if (nEltColl==1) CentroidColl[iElt]=1;

    // Check nodal collocation points.
    else
    {
      for (unsigned int iColl=0; iColl<nEltColl; iColl++)
      {
        unsigned int CollNode = (unsigned int)(Elt[(2+iColl)*nElt+iElt]);
        int CollIndex;
        BemNodeIndex(Nod,nNod,CollNode,CollIndex);
        NodalColl[CollIndex]=1;
      }
    }
  }

  // Count centroid collocation points.
  nCentroidColl=0;
  for (unsigned int iColl=0; iColl<nElt; iColl++)
  {
    if (CentroidColl[iColl]==1) ++nCentroidColl;
  }

  //mexPrintf("nCentroidColl: %d\n",nCentroidColl); // DEBUG
  // Count nodal collocation points.
  nNodalColl=0;
  for (unsigned int iColl=0; iColl<nNod; iColl++)
  {
    if (NodalColl[iColl]==1) ++nNodalColl;
  }
//  mexPrintf("nNodalColl: %d\n",nNodalColl); // DEBUG
}


void BemCollCoords(const double* const Elt, const double* const Nod,
                   unsigned int* const  TypeID, unsigned int* const  nKeyOpt,
                   const char* const TypeName[], const char* const TypeKeyOpts[],
                   const unsigned int& nEltType,
                   const unsigned int* const CentroidColl,
                   const unsigned int* const NodalColl,double* const CollPoints,
                   const unsigned int& nTotalColl, const unsigned int& nElt,const unsigned int& nNod)
/* BemCollCoords: Retreive collocation point coordinates.
 *
 *  CollPoints[0*nTotalColl+iColl]  Collocation point type (1 or 2)
 *  CollPoints[1*nTotalColl+iColl]  Element number or Node number
 *  CollPoints[2*nTotalColl+iColl]  x-coordinate
 *  CollPoints[3*nTotalColl+iColl]  y-coordinate
 *  CollPoints[4*nTotalColl+iColl]  z-coordinate
 *
 */
{
  // --- CENTROID COLLOCATION POINTS ---
  unsigned int iColl=0;
  unsigned int EltType;
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

  double* const xiCentroid = new(nothrow) double[2];
  if (xiCentroid==0) throw("Out of memory.");
  unsigned int nXiCentroid=1;

  for (unsigned int iElt=0; iElt<nElt; iElt++)
  {
    if (CentroidColl[iElt]==1)
    {
      // --- Compute element centroid ---
      EltType = (unsigned int)(Elt[nElt+iElt]);
      eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
             nEltType,Parent,nEltNod,nEltColl,ShapeTypeN,ShapeTypeM,EltDim,AxiSym,Periodic,
             nGauss,nEltDiv,nGaussSing,nEltDivSing);
      
      
      if (Parent==1)
      {
        xiCentroid[0]=3.333333333333333e-01;
        xiCentroid[1]=3.333333333333333e-01;
      }
      else
      {
        xiCentroid[0]=0.0;
        xiCentroid[1]=0.0;
      }
      
      double* const N = new(nothrow) double[nEltNod];
      if (N==0) throw("Out of memory.");
      shapefun(ShapeTypeN,nXiCentroid,xiCentroid,N);
      
      double* const CollCoord = new(nothrow) double[3];
      if (CollCoord==0) throw("Out of memory.");
      for (unsigned int iCoord=0; iCoord<3; iCoord++) CollCoord[iCoord]=0.0;
      
      for (unsigned int iNod=0; iNod<nEltNod; iNod++)
      {
        unsigned int NodeNumber = (unsigned int)(Elt[(2+iNod)*nElt+iElt]);
        int NodeIndex=0;
        BemNodeIndex(Nod,nNod,NodeNumber,NodeIndex);
//         mexPrintf("NodeIndex: %d\n",NodeIndex); // DEBUG
        for (unsigned int iCoord=0; iCoord<3; iCoord++)
        {
          CollCoord[iCoord]+= Nod[nNod*(1+iCoord)+NodeIndex]*N[iNod];
//           mexPrintf("nNod*(1+iCoord)+NodeIndex: %d\n",nNod*(1+iCoord)+NodeIndex); // DEBUG
        }
      }

      CollPoints[0*nTotalColl+iColl]=1.0;          // type 1
      CollPoints[1*nTotalColl+iColl]=Elt[iElt];    // Element number
      CollPoints[2*nTotalColl+iColl]=CollCoord[0]; // x-coordinate
      CollPoints[3*nTotalColl+iColl]=CollCoord[1]; // y-coordinate
      CollPoints[4*nTotalColl+iColl]=CollCoord[2]; // z-coordinate
      iColl++;
//       mexPrintf("iColl: %d\n",iColl); // DEBUG

      delete [] N;
      delete [] CollCoord;
    }
  }

  // --- NODAL COLLOCATION POINTS ---
  for (unsigned int iNod=0; iNod<nNod; iNod++)
  {
    if (NodalColl[iNod]==1)
    {
      CollPoints[0*nTotalColl+iColl]=2.0;               // type 2: nodal collocation
      CollPoints[1*nTotalColl+iColl]=Nod[0*nNod+iNod];  // Node number
      CollPoints[2*nTotalColl+iColl]=Nod[1*nNod+iNod];  // x-coordinate
      CollPoints[3*nTotalColl+iColl]=Nod[2*nNod+iNod];  // y-coordinate
      CollPoints[4*nTotalColl+iColl]=Nod[3*nNod+iNod];  // z-coordinate
      iColl++;
    }
  }
  delete [] xiCentroid;
}



void BemCoincNodes(const double* const Nod, const unsigned int& nNod,
                   double* const CoincNod, bool& SlavesExist)
/*
 *  Retreive coincident nodes from a given node list.
 */
{
  double diff;
  const double CoincEps=1.0e-10;
  SlavesExist =false;
  for (unsigned int iNod=0; iNod<2*nNod; iNod++) CoincNod[iNod]=0.0;
  for (unsigned int iNod=0; iNod<nNod-1; iNod++)
  {
    if (CoincNod[iNod]==0)
    {
      for (unsigned int jNod=iNod+1; jNod<nNod; jNod++)
      {
        diff=fabs(Nod[nNod+iNod]-Nod[nNod+jNod]);
        if (diff<CoincEps)
        {
          diff=fabs(Nod[2*nNod+iNod]-Nod[2*nNod+jNod]);
          if (diff<CoincEps)
          {
            diff=fabs(Nod[3*nNod+iNod]-Nod[3*nNod+jNod]);
            if (diff<CoincEps)
            {
              CoincNod[jNod]=1.0;
              CoincNod[nNod+jNod]=Nod[iNod];
              SlavesExist = true;
            }
          }
        }
      }
    }
  }
}

void BemEltCollIndex(const double* const Elt, const unsigned int& iElt, const unsigned int& nElt,
                     const double* const CollPoints, const unsigned int& nCentroidColl,
                     const unsigned int& nTotalColl, const unsigned int& nEltColl,
                     const unsigned int& nEltNod,unsigned int* const eltCollIndex)
/*
 *   BemEltCollIndex looks up the collocation point index "iColl"
 *   for all collocation points of the element with index "iElt".
 *
 *      => eltCollIndex[iEltColl] = iColl;
 *
 *   eltCollIndex is a (nEltColl,1) integer array.
 *
 */
{
  if (nEltColl==1) // Centroid collocation
  {
    for (unsigned int iColl=0; iColl<nCentroidColl; iColl++)
    {
      if (CollPoints[nTotalColl+iColl]==(unsigned int)(Elt[iElt])) eltCollIndex[0]=iColl;
    }
  }
  else  // Nodal collocation
  {
    for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
    {
      for (unsigned int iColl=nCentroidColl; iColl<nTotalColl; iColl++)
      {
        if (CollPoints[nTotalColl+iColl]==(unsigned int)(Elt[nElt*(2+iEltNod)+iElt])) eltCollIndex[iEltNod]=iColl;
      }
    }
  }
}

void BemRegularColl(const double* const Elt,const unsigned int& iElt,
                    const unsigned int& nElt, const double* const Nod,
                    const unsigned int& nNod, const double* const CoincNod,
                    const bool& SlavesExist, const double* const CollPoints,
                    const unsigned int& nCentroidColl, const unsigned int& nTotalColl,
                    unsigned int* const RegularColl,unsigned int& nRegularColl,unsigned int& nSingularColl,
                    const unsigned int* const  TypeID, const unsigned int* const  nKeyOpt,
                    const char* const TypeName[], const char* const TypeKeyOpts[],
                    const unsigned int& nEltType)
 /*
  * Checks, for one element iElt, which collocation points are regular,
  * i.e. do not belong to the element.
  */
{
  for (unsigned int iColl=0; iColl<nTotalColl; iColl++) RegularColl[iColl]=1;
  nSingularColl=0;

  // --- Check Centroids.
  for (unsigned int iColl=0; iColl<nCentroidColl; iColl++)
  {
    if ((unsigned int)(CollPoints[nTotalColl+iColl])==Elt[iElt])
    {
      RegularColl[iColl]=0;
      nSingularColl++;
    }
  }

  // --- Check Nodal collocation points.
  unsigned int NodID;
  int NodIndex;
  unsigned int EltParent;
  unsigned int nEltNod;
  unsigned int nEltColl;
  unsigned int EltShapeN;
  unsigned int EltShapeM;
  unsigned int EltDim;
  unsigned int AxiSym;
  unsigned int Periodic;
  unsigned int nGauss;
  unsigned int nEltDiv;
  unsigned int nGaussSing;
  unsigned int nEltDivSing;

  unsigned int masterNodID;

  unsigned int EltType = (unsigned int)(Elt[nElt+iElt]);
  eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
         nEltType,EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,
         nGauss,nEltDiv,nGaussSing,nEltDivSing);

  for (unsigned int iNod=0; iNod<nEltNod; iNod++)
  {
    NodID = (unsigned int)(Elt[(2+iNod)*nElt+iElt]);
    BemNodeIndex(Nod,nNod,NodID,NodIndex);
    for (unsigned int iColl=nCentroidColl; iColl<nTotalColl; iColl++)
    {
      if ( ((unsigned int)(CollPoints[nTotalColl+iColl])==NodID) && (!(RegularColl[iColl]==0)))
      {
        RegularColl[iColl]=0;
        RegularColl[nTotalColl+iColl]=iNod;
        nSingularColl++;
      }
    }
    // --- Check for other coincident nodes and remove if collocation point.
    if (SlavesExist)
    {
      if (CoincNod[NodIndex]==1)   // node is a SLAVE node.
      {
        masterNodID = int(CoincNod[nNod+NodIndex]);
      }
      else masterNodID = NodID;     // else it is a MASTER node.
      
      for (unsigned int jNod=0; jNod<nNod; jNod++)
      {
        if (((CoincNod[jNod]==1) && ((unsigned int)(CoincNod[nNod+jNod])==masterNodID)) || ((unsigned int)(Nod[jNod])==masterNodID))
        // The node is a slave node with the same master or the master itsself
        // ==> is singular as a collocation point.
        {
          for (unsigned int iColl=nCentroidColl; iColl<nTotalColl; iColl++)
          {
            if (((unsigned int)(CollPoints[nTotalColl+iColl])==(unsigned int)(Nod[jNod])) && (!(RegularColl[iColl]==0)))
            {
              RegularColl[iColl]=0;
              RegularColl[nTotalColl+iColl]=iNod;
              nSingularColl++;
            }
          }
        }
      }
    }
  }
  nRegularColl=nTotalColl-nSingularColl;
}
