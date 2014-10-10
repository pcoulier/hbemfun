#include <math.h>
#include <new>

#include "mex.h"

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
void uniquecoll(const unsigned int& ms,
                const unsigned int& ns,        
                unsigned int* scolli,
                unsigned int* uniquescolli,
                unsigned int* Nuniquescolli,
                unsigned int* nuniquescolli,
                unsigned int* uniquescolliind)
//==============================================================================
{          

       Nuniquescolli[0] = 0;  

       // Get unique collocation points
       // Determine the total number of unique collocation points
       // Determine the number of each unique collocation point
       for (unsigned int iColli=0; iColli<ms*ns; iColli++)
       {

            unsigned int iuniquescolli=0;
            unsigned int keyUnique=0;

            while (keyUnique !=1 && iuniquescolli<Nuniquescolli[0])
            {
               if (scolli[iColli]==uniquescolli[iuniquescolli])
               {
                   nuniquescolli[iuniquescolli]++;                                      
                   keyUnique=1;                   
               }
               else
               {
                 iuniquescolli++;  
               }  
            }
            
            if (keyUnique !=1)
            {
               uniquescolli[Nuniquescolli[0]]=scolli[iColli];          
               nuniquescolli[Nuniquescolli[0]]=1;
               Nuniquescolli[0]++;
            }
        }
        
		// mexPrintf("ms: %d \n",ms);
		// mexPrintf("ns: %d \n",ns);
		
        
        // Determine the indices of the unique collocation points    
        unsigned int nuniquescollicumsum = 0;    
            
        for (unsigned int iuniquescolli=0;iuniquescolli<Nuniquescolli[0];iuniquescolli++)
        {
		   // mexPrintf("iuniquescolli: %d \n",iuniquescolli);
           unsigned int iscolli = 0;
           unsigned int inuniquescolli = 0;
           
           while (iscolli<ms*ns && inuniquescolli<nuniquescolli[iuniquescolli])
           {	
				// mexPrintf("iscolli: %d \n",iscolli);
                 if (scolli[iscolli] == uniquescolli[iuniquescolli])
                 {
                    uniquescolliind[nuniquescollicumsum+inuniquescolli]=iscolli; 
                    inuniquescolli++;                                   
                 }
                 // iscolli++;
				// iscolli = 
				if ((iscolli+ms) < (ms*ns)){iscolli+=ms;}
				else {iscolli=(iscolli + ms) % (ms*ns)+1;}
				
           }
          nuniquescollicumsum+=nuniquescolli[iuniquescolli];  
        }

}
