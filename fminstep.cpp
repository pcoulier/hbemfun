#include <new>
#include <valarray>
#include <algorithm>
using namespace std;

bool sortUtil(double* i, double* j)
{
  return(*i < *j);
}
void sort_fv(valarray<double>& v, valarray<double>& fv, const int& n)
{
  double** sortIndex=new(nothrow) double*[n+1];
  for (int i=0; i<n+1; i++) sortIndex[i]=&fv[i];
  int* sortIndex2=new(nothrow) int[n+1];
  sort(sortIndex,sortIndex+n+1,sortUtil);
  valarray<double> fvutil(0.0,n+1);
  valarray<double> vutil(0.0,n*(n+1));
  for (int i=0; i<n+1; i++)
  {
    fvutil[i]=*sortIndex[i];
    sortIndex2[i]=(int)((sortIndex[i]-&fv[0]));
    for (int j=0; j<n; j++)
    {
      vutil[n*i+j]=v[n*sortIndex2[i]+j];
    }
  }
  for (int i=0; i<n+1; i++) fv[i]=fvutil[i];
  for (int i=0; i<n*(n+1); i++) v[i]=vutil[i];
  delete [] sortIndex;
  delete [] sortIndex2;
}

double fminstep(double(*funToMinimize)(valarray<double>&, const void** const),
                valarray<double>& x, const valarray<double>& xRes,
                const valarray<double>& xTol, const int& nMaxIter,
                const void** const varargin)
{
  const int n=x.size();
  const double rho=1.0;
  const double chi=2.0;
  const double psi=0.5;
  const double sigma=0.5;
  
  const valarray<double> xin=x;
  valarray<double> v(0.0,n*(n+1));
  valarray<double> fv(0.0,n+1);
  v[slice(0,n,1)]=xin;
  
  fv[0]=(*funToMinimize)(x,varargin);
  for (int j=0; j<n; j++)
  {
    valarray<double> y=xin;
    y[j]=y[j]+xRes[j];
    v[slice(n*(j+1),n,1)]=y;
    x=y;
    fv[j+1]=(*funToMinimize)(x,varargin);
  }
  
  sort_fv(v,fv,n);
  
  bool converged=false;
  int nIter=0;
  
  // CONVERGENCE CHECK
  valarray<double> convAux1(0.0,n*n);
  for (int j=0; j<n; j++)  convAux1[slice(n*j,n,1)] = 
                     abs(((valarray<double>)v[slice(n*j+n,n,1)])
                                        -((valarray<double>)v[slice(0,n,1)]));
  bool anyTooLarge=false;                                      
  for (int j=0; j<n; j++)
  {
    double ConvCheck=((valarray<double>)convAux1[slice(j,n,n)]).max();
    if (ConvCheck>xTol[j]) anyTooLarge=true;
  }
  if (!(anyTooLarge)) converged=true;
  
  
  while (!(converged))
  {
    nIter++;
    valarray<double> xbar(0.0,n);
    for (int i=0; i<n; i++) xbar[i]=((valarray<double>)v[slice(i,n,n)]).sum()/n;
    
    valarray<double> xr=(1.0+rho)*xbar-rho*(valarray<double>)v[slice(n*n,n,1)];
    x=xr;
    double fxr=(*funToMinimize)(x,varargin);
    
    if (fxr<fv[0])
    {
      valarray<double> xe=(1+rho*chi)*xbar
                                   -rho*chi*(valarray<double>)v[slice(n*n,n,1)];
      x=xe;
      double fxe=(*funToMinimize)(x,varargin);
    
      if (fxe<fxr)
      {
        v[slice(n*n,n,1)]=xe;
        fv[n]=fxe;
      }
      else
      {
        v[slice(n*n,n,1)]=xr;
        fv[n]=fxr;
      }
    }
    else
    {
      if (fxr<fv[n-1])
      {
        v[slice(n*n,n,1)]=xr;
        fv[n]=fxr;
      }
      else
      {
        bool shrink=false;
        if (fxr<fv[n])
        {
          valarray<double> xc=(1+psi*rho)*xbar
                                   -psi*rho*(valarray<double>)v[slice(n*n,n,1)];
          x=xc;
          double fxc=(*funToMinimize)(x,varargin);
    
          if (fxc<=fxr)
          {
            v[slice(n*n,n,1)]=xc;
            fv[n]=fxc;
            shrink=false;
          }
          else shrink=true;
        }
       else
        {
          valarray<double> xcc=(1-psi)*xbar
                                       +psi*(valarray<double>)v[slice(n*n,n,1)];
          x=xcc;
          double fxcc=(*funToMinimize)(x,varargin);
    
          if (fxcc<fv[n])
          {
            v[slice(n*n,n,1)]=xcc;
            fv[n]=fxcc;
            shrink=false;
          }
          else shrink=true;
        }
        if (shrink)
        {
          for (int j=1; j<n+1; j++)
          {
            v[slice(n*j,n,1)]=(valarray<double>)v[slice(n*j,n,1)]
                                +sigma*((valarray<double>)v[slice(n*j,n,1)]
                                            -(valarray<double>)v[slice(0,n,1)]);
            x=v[slice(n,n,1)];
            fv[j]=(*funToMinimize)(x,varargin);
          }
        }
      }
    }    
    sort_fv(v,fv,n);
    
    // CONVERGENCE CHECK
    for (int j=0; j<n; j++)  convAux1[slice(n*j,n,1)] = 
                       abs(((valarray<double>)v[slice(n*j+n,n,1)])
                                          -((valarray<double>)v[slice(0,n,1)]));
    anyTooLarge=false;                                      
    for (int j=0; j<n; j++)
    {
      double ConvCheck=((valarray<double>)convAux1[slice(j,n,n)]).max();
      if (ConvCheck>xTol[j]) anyTooLarge=true;
    }
    if (!(anyTooLarge)) converged=true;
    if (nIter>=nMaxIter) converged=true;
  }
  
  // ASSIGN STUFF
  x=v[slice(0,n,1)];  // Location of the minimium
  return fv.min();    // Function value at minimum
}
