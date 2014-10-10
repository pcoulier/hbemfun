/* rsa.cpp
 *
 * Mattias Schevenels
 * May 2009
 */

#include <string>
#include <stdlib.h>
#include "bigint.h"
#include <cstdio>

using namespace std;

/******************************************************************************/

void modpow(BigUnsigned b, BigUnsigned e, const BigUnsigned& m, BigUnsigned& r)
{
  r = 1;
  BigUnsigned one = 1;
  BigUnsigned zero = 0;
  while (e > zero)
  {
    if ((e & one) == one) r = (r * b) % m;
    e >>= 1;
    b = (b * b) % m;
  }
}

/******************************************************************************/

void hex2asc(const string& hex, string& asc)
{
  const int n = hex.length()/2;
  asc.clear();
  asc.reserve(n);
  string s;
  for (int i=0; i<n; i++)
  {
    s=strtoul(hex.substr(2*i,2).c_str(),0,16);
    asc.append(s);
  }
}

/******************************************************************************/

void asc2hex(const string& asc, string& hex)
{
  const int n = asc.length();

  hex.clear();
  hex.reserve(2*n);
  char c[3];
  for (int i=0; i<n; i++)
  {
    int k = asc[i];
    if (k<0) k+=256;
    sprintf(c,"%02X",k);
    c[2]='\0';
    hex.append(c);
  }
}

/******************************************************************************/

void rsadec(const string& in, const string& exp, const string& mod, string& out)
{
  BigUnsigned inBU;
  BigUnsigned expBU(BigUnsignedInABase(exp.c_str(),16));
  BigUnsigned modBU(BigUnsignedInABase(mod.c_str(),16));
  BigUnsigned outBU;
  string outHex;

  const int modSize = mod.length();
  const int inSize = in.length();
  const int nChunk = (inSize-1)/modSize+1;

  out.clear();

  for (int iChunk=0; iChunk<nChunk; iChunk++)
  {
    inBU = BigUnsignedInABase(in.substr(modSize*iChunk,modSize).c_str(),16);
    modpow(inBU,expBU,modBU,outBU);
    outHex=string(BigUnsignedInABase(outBU,16));
    out.append(outHex);
  }
}

/******************************************************************************/

void rsaenc(const string& in, const string& exp, const string& mod, string& out)
{
  BigUnsigned inBU;
  BigUnsigned expBU(BigUnsignedInABase(exp.c_str(),16));
  BigUnsigned modBU(BigUnsignedInABase(mod.c_str(),16));
  BigUnsigned outBU;
  string inHex;
  string outHex;

  const int modSize = mod.length();
  const int inSize = in.length();
  const int outSize = inSize;
  const int nChunk = (outSize-1)/modSize+1;

  out.clear();

  for (int iChunk=0; iChunk<nChunk; iChunk++)
  {
    //asc2hex(in.substr(modSize/2*iChunk,modSize/2).c_str(),inHex);
    inHex=in.substr(modSize*iChunk,modSize).c_str();
    inHex.insert(0,modSize-inHex.length(),'0');
    inBU = BigUnsignedInABase(inHex,16);
    modpow(inBU,expBU,modBU,outBU);
    outHex=string(BigUnsignedInABase(outBU,16));
    outHex.insert(0,modSize-outHex.length(),'0');
    out.append(outHex);
  }
}
