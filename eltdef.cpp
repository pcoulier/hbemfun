#include <string.h>

#ifndef __GNUC__
#define strcasecmp _strcmpi
#endif

void eltdef(const unsigned int& EltType, const unsigned int* const  TypeID,
            const char* const TypeName[], const char* const TypeKeyOpts[],
            const unsigned int* const  nKeyOpt, const unsigned int& nEltType, unsigned int& Parent,
            unsigned int& nNod, unsigned int& nColl, unsigned int& ShapeTypeN, unsigned int& ShapeTypeM,
            unsigned int& EltDim, unsigned int& AxiSym, unsigned int& Periodic, unsigned int& nGauss, unsigned int& nEltDiv,
            unsigned int& nGaussSing, unsigned int& nEltDivSing)
/*  ELTDEF Element properties.
 *    EltType    (int,in)  Element type reference number.
 *    Parent     (int,out) Parent mapping: 1 for triangle, 2 for quadrilateral,
 *                         0 for line element.
 *    nNod       (int,out) Number of nodes.
 *    nColl      (int,out) Number of collocation points.
 *    ShapeTypeN (int,out) Shape type for the interpolation of the geometry. See
 *                         shapefun.cpp for definitions.
 *    ShapeTypeM (int,out) Shape type for the interpolation of the field
 *                         variables (displacements, tractions). See
 *                         shapefun.cpp for definitions.
 *    EltDim     (int,out) Element dimension. 1 for a 1D element, 2 for a
 *                         2D element.
 *    AxiSym     (int,out) Axisymmetric key. 1 for axisymmetric, 0 otherwise
 */
{
  int TypeInd=-1;
  for (unsigned int iType=0; iType< nEltType; iType++) if (TypeID[iType] == EltType)  TypeInd=iType;
  if (TypeInd==-1) throw("Element type not found in type cell array.");

  if (strcasecmp(TypeName[TypeInd],"quad4")==0)
  {
    Parent=2;
    nNod=4;
    ShapeTypeN=4;
    EltDim=2;
  }
  else if (strcasecmp(TypeName[TypeInd],"quad8")==0)
  {
    Parent=2;
    nNod=8;
    ShapeTypeN=5;
    EltDim=2;
  }
  else if (strcasecmp(TypeName[TypeInd],"quad9")==0)
  {
    Parent=2;
    nNod=9;
    ShapeTypeN=6;
    EltDim=2;
  }
  else if (strcasecmp(TypeName[TypeInd],"tria3")==0)
  {
    Parent=1;
    nNod=3;
    ShapeTypeN=2;
    EltDim=2;
  }
  else if (strcasecmp(TypeName[TypeInd],"tria6")==0)
  {
    Parent=1;
    nNod=6;
    ShapeTypeN=3;
    EltDim=2;
  }
  else if (strcasecmp(TypeName[TypeInd],"line2")==0)
  {
    Parent=0;
    nNod=2;
    ShapeTypeN=7;
    EltDim=1;
  }
  else if (strcasecmp(TypeName[TypeInd],"line3")==0)
  {
    Parent=0;
    nNod=3;
    ShapeTypeN=8;
    EltDim=1;
  }
  else if (strcasecmp(TypeName[TypeInd],"line4")==0)
  {
    Parent=0;
    nNod=4;
    nColl=1;
    ShapeTypeN=9;
    ShapeTypeM=1;
    EltDim=1;
  }
  else throw("unknown element type.");

  // PROCESS KEY OPTIONS
  nColl=1;
  ShapeTypeM=1;
  AxiSym=0;
  Periodic=0;
  nGauss=6;
  nGaussSing=6;
  nEltDiv=1;
  nEltDivSing=1;

  for (unsigned int iKeyOpt=0; iKeyOpt<nKeyOpt[TypeInd]; iKeyOpt++)
  {
    if (strcasecmp(TypeKeyOpts[TypeInd+nEltType*iKeyOpt],"nodcol")==0)
    {
      nColl=nNod;
      ShapeTypeM=ShapeTypeN;
    }

    if (strcasecmp(TypeKeyOpts[TypeInd+nEltType*iKeyOpt],"eltcol")==0)
    {
      nColl=1;
      ShapeTypeM=1;
    }

    // ONLY ALLOW AXISYMMETRIC KEYOPTION FOR 2D ELEMENTS.
    char strutil[512];
    strncpy(strutil,TypeName[TypeInd],4);
    strutil[4]='\0';
    if (strcasecmp(strutil,"line")==0)
    {
      if (strcasecmp(TypeKeyOpts[TypeInd+nEltType*iKeyOpt],"axisym")==0) AxiSym=1;
      if (strcasecmp(TypeKeyOpts[TypeInd+nEltType*iKeyOpt],"periodic")==0) throw("Periodic keyoption not allowed for this element type.");
    }
    else
    {
      if (strcasecmp(TypeKeyOpts[TypeInd+nEltType*iKeyOpt],"axisym")==0)  throw("Axisymmetric keyoption not allowed for this element type.");
      if (strcasecmp(TypeKeyOpts[TypeInd+nEltType*iKeyOpt],"periodic")==0) Periodic=1;
    }

    // NUMBER OF GAUSSIAN POINTS nGaussXX
    int optLen=strlen(TypeKeyOpts[TypeInd+nEltType*iKeyOpt]);
    if (optLen>512) throw("keyoption string argument too large");
    if (optLen>6)
    {
      strncpy(strutil,TypeKeyOpts[TypeInd+nEltType*iKeyOpt],6);
      strutil[6]='\0';

      if (strcasecmp(strutil,"nGauss")==0)
      {
        strncpy(strutil,TypeKeyOpts[TypeInd+nEltType*iKeyOpt],optLen);
        strutil[optLen]='\0';

        // CHECK IF FOLLOWING CHARACTER IS NUMERIC
        int intutil1=((int)(strutil[6])-48);
        if ((intutil1>=0)&&(intutil1<10))
        {
          nGauss=0;
          for (int iChar=6;iChar<optLen;iChar++)
          {
             const int intutil=(int)(strutil[iChar]) - 48;
             if ((intutil>=0)&&(intutil<10)) nGauss=10*nGauss + intutil;
          }
        }
      }
    }
    if  (nGauss<1) throw("Error reading number of Gaussian points");

    if (optLen>10)
    {
      strncpy(strutil,TypeKeyOpts[TypeInd+nEltType*iKeyOpt],10);
      strutil[10]='\0';

      if (strcasecmp(strutil,"nGaussSing")==0)
      {
        strncpy(strutil,TypeKeyOpts[TypeInd+nEltType*iKeyOpt],optLen);
        strutil[optLen]='\0';

        // CHECK IF FOLLOWING CHARACTER IS NUMERIC
        int intutil1=((int)(strutil[10])-48);
        if ((intutil1>=0)&&(intutil1<10))
        {
          nGaussSing=0;
          for (int iChar=10;iChar<optLen;iChar++)
          {
             const int intutil=(int)(strutil[iChar]) - 48;
             if ((intutil>=0)&&(intutil<10)) nGaussSing=10*nGaussSing + intutil;
          }
        }

      }
    }
    if  (nGaussSing<1) throw("Error reading number of singular Gaussian points");

    if (optLen>7)
    {
      strncpy(strutil,TypeKeyOpts[TypeInd+nEltType*iKeyOpt],7);
      strutil[7]='\0';
      if (strcasecmp(strutil,"nEltDiv")==0)
      {
        strncpy(strutil,TypeKeyOpts[TypeInd+nEltType*iKeyOpt],optLen);
        strutil[optLen]='\0';

        // CHECK IF FOLLOWING CHARACTER IS NUMERIC
        int intutil1=((int)(strutil[7])-48);
        if ((intutil1>=0)&&(intutil1<10))
        {
          nEltDiv=0;
          for (int iChar=7;iChar<optLen;iChar++)
          {
             const int intutil=(int)(strutil[iChar]) - 48;
             if ((intutil>=0)&&(intutil<10)) nEltDiv=10*nEltDiv + intutil;
          }
        }
      }
    }
    if  (nEltDiv<1) throw("Error reading number of element divisions");

    if (optLen>11)
    {
      strncpy(strutil,TypeKeyOpts[TypeInd+nEltType*iKeyOpt],11);
      strutil[11]='\0';
      if (strcasecmp(strutil,"nEltDivSing")==0)
      {
        strncpy(strutil,TypeKeyOpts[TypeInd+nEltType*iKeyOpt],optLen);
        strutil[optLen]='\0';

        int intutil1=((int)(strutil[11])-48);
        if ((intutil1>=0)&&(intutil1<10))
        {
          nEltDivSing=0;
          for (int iChar=11;iChar<optLen;iChar++)
          {
             const int intutil=(int)(strutil[iChar]) - 48;
             if ((intutil>=0)&&(intutil<10)) nEltDivSing=10*nEltDivSing + intutil;
          }
        }
      }
    }
    if  (nEltDivSing<1) throw("Error reading number of singular element divisions");
  }
}


void eltnoddef(const unsigned int& EltType,const unsigned int* const  TypeID,
               const char* const TypeName[], const unsigned int& nEltType,
               double* const eltNodXi)
/*  ELTNODDEF Node definitions per element type.
 *
 *    EltType    (unsigned int,in)  Element type reference number.
 *    EltNodXi   (dp[nNod*EltDim],out) Node coordinates in natural coordinates.
 *                For 2D elements, the coordinates of node iNod=0:1:nNod-1
 *                are (eta=eltNodXi[iNod],xi=eltNodXi[nNod+iNod]).
 *                For 1D elements, eta=eltNodXi[iNod].
 */
{
  int TypeInd=-1;
  for (unsigned int iType=0; iType< nEltType; iType++) if (TypeID[iType] == EltType)  TypeInd=iType;
  if (TypeInd==-1) throw("Element type not found in type cell array.");

  if (strcasecmp(TypeName[TypeInd],"tria3")==0)
  {
    eltNodXi[0]=0.0;
    eltNodXi[1]=1.0;
    eltNodXi[2]=0.0;
    eltNodXi[3]=0.0;
    eltNodXi[4]=0.0;
    eltNodXi[5]=1.0;
  }
  else if (strcasecmp(TypeName[TypeInd],"quad4")==0)
  {
    eltNodXi[0]=-1.0;
    eltNodXi[1]= 1.0;
    eltNodXi[2]= 1.0;
    eltNodXi[3]=-1.0;
    eltNodXi[4]=-1.0;
    eltNodXi[5]=-1.0;
    eltNodXi[6]= 1.0;
    eltNodXi[7]= 1.0;
  }
  else if (strcasecmp(TypeName[TypeInd],"tria6")==0)
  {
    eltNodXi[0] =0.0;
    eltNodXi[1] =1.0;
    eltNodXi[2] =0.0;
    eltNodXi[3] =0.5;
    eltNodXi[4] =0.5;
    eltNodXi[5] =0.0;
    eltNodXi[6] =0.0;
    eltNodXi[7] =0.0;
    eltNodXi[8] =1.0;
    eltNodXi[9] =0.0;
    eltNodXi[10]=0.5;
    eltNodXi[11]=0.5;
  }
  else if (strcasecmp(TypeName[TypeInd],"quad8")==0)
  {
    eltNodXi[0] =-1.0;
    eltNodXi[1] = 1.0;
    eltNodXi[2] = 1.0;
    eltNodXi[3] =-1.0;
    eltNodXi[4] = 0.0;
    eltNodXi[5] = 1.0;
    eltNodXi[6] = 0.0;
    eltNodXi[7] =-1.0;
    eltNodXi[8] =-1.0;
    eltNodXi[9] =-1.0;
    eltNodXi[10]= 1.0;
    eltNodXi[11]= 1.0;
    eltNodXi[12]=-1.0;
    eltNodXi[13]= 0.0;
    eltNodXi[14]= 1.0;
    eltNodXi[15]= 0.0;
  }
  else if (strcasecmp(TypeName[TypeInd],"quad9")==0)
  {
    eltNodXi[0] =-1.0;
    eltNodXi[1] = 1.0;
    eltNodXi[2] = 1.0;
    eltNodXi[3] =-1.0;
    eltNodXi[4] = 0.0;
    eltNodXi[5] = 1.0;
    eltNodXi[6] = 0.0;
    eltNodXi[7] =-1.0;
    eltNodXi[8] = 0.0;
    eltNodXi[9] =-1.0;
    eltNodXi[10]=-1.0;
    eltNodXi[11]= 1.0;
    eltNodXi[12]= 1.0;
    eltNodXi[13]=-1.0;
    eltNodXi[14]= 0.0;
    eltNodXi[15]= 1.0;
    eltNodXi[16]= 0.0;
    eltNodXi[17]= 0.0;
  }
  else if (strcasecmp(TypeName[TypeInd],"line2")==0)
  {
    eltNodXi[0] =-1.0;
    eltNodXi[1] = 1.0;
  }
  else if (strcasecmp(TypeName[TypeInd],"line3")==0)
  {
    eltNodXi[0] =-1.0;
    eltNodXi[1] = 0.0;
    eltNodXi[2] = 1.0;
  }
  else if (strcasecmp(TypeName[TypeInd],"line4")==0)
  {
    eltNodXi[0] =-1.0;
    eltNodXi[1] =-0.33333333333333;
    eltNodXi[2] = 0.33333333333333;
    eltNodXi[3] = 1.0;
  }
  else throw("Element type undefined");
}
