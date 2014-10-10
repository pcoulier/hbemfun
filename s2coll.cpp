#include <math.h>
#include <new>
#include "uniquecoll.h"
#include "mex.h"

#ifndef _int64_
typedef long long int int64;
typedef unsigned long long int uint64;
#endif

using namespace std;

inline double sign(const double& a)
{
  if (a==0.0) return 1.0;
  else return (a>0.0 ? 1.0 : -1.0);
}
inline double sqr(const double& a)
{
  return a*a;
}
//==============================================================================
void s2coll(const double* const s,
            const unsigned int& ms,
            const unsigned int& ns,
            const unsigned int& nDof,
            const unsigned int& probDim,
            unsigned int* scolli,
            unsigned int* scompi,
            unsigned int* uniquescolli,
            unsigned int* Nuniquescolli,
            unsigned int* nuniquescolli,
            unsigned int* uniquescolliind,
            unsigned int* scollj,
            unsigned int* scompj,
            unsigned int* uniquescollj,
            unsigned int* Nuniquescollj,
            unsigned int* nuniquescollj,
            unsigned int* uniquescolljind,
			bool* scolliOnDiag)
//==============================================================================
{          
           // Get indices of collocation points i and collocation points j
           for (unsigned int iColli=0; iColli<ms; iColli++)
           {
               for (unsigned int iCollj=0; iCollj<ns; iCollj++)
               {
//                    unsigned int a = (unsigned int) s[iColli+ms*iCollj]; 
//                    unsigned int b = (unsigned int) s[iColli+ms*iCollj]; 
                   uint64 a = (uint64) s[iColli+ms*iCollj]; 
                   uint64 b = (uint64) s[iColli+ms*iCollj]; 
                   
                   scolli[iColli+ms*iCollj] = (a-1) % nDof; 
                   scollj[iColli+ms*iCollj] = (b-1) / nDof; 

                   scompi[iColli+ms*iCollj] = scolli[iColli+ms*iCollj] % probDim;
                   scompj[iColli+ms*iCollj] = scollj[iColli+ms*iCollj] % probDim;
                   
                   scolli[iColli+ms*iCollj] = scolli[iColli+ms*iCollj] / probDim; 
                   scollj[iColli+ms*iCollj] = scollj[iColli+ms*iCollj] / probDim; 
				   
				   if (scolli[iColli+ms*iCollj]  == scollj[iColli+ms*iCollj])
				   {
						scolliOnDiag[iColli+ms*iCollj]=1;
				   }
				   else
				   {
						scolliOnDiag[iColli+ms*iCollj]=0;
				   }
               }
           }
           
           // Get unique collocation points
           uniquecoll(ms,ns,scolli,uniquescolli,Nuniquescolli,nuniquescolli,uniquescolliind);  
		   uniquecoll(ms,ns,scollj,uniquescollj,Nuniquescollj,nuniquescollj,uniquescolljind); 
}
