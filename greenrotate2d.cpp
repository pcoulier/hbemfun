//==============================================================================
void greenrotate2d(const double* const normal,const unsigned int& iXi,
                   const unsigned int& nGrSet,const unsigned int& ntgComp,
                   const bool& tgCmplx,const bool& tg0Cmplx,
                   const double* const TgrRe, const double* const TgrIm,
                   const double* const Tgr0Re, const double* const Tgr0Im,
                   double* const TXiRe, double* const TXiIm, double* const TXi0Re,
                   double* const TXi0Im, const bool& TmatOut)
//==============================================================================
{
  const bool calcTg0=Tgr0Re!=0;
  if (TmatOut)
  {
    for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
    {
      if (ntgComp==2)      // 2D out-of-plane
      {
        TXiRe[1*iGrSet+0]=TgrRe[2*iGrSet+0]*normal[3*iXi+0]+TgrRe[2*iGrSet+1]*normal[3*iXi+2];     // tyy
        if(calcTg0) TXi0Re[1*iGrSet+0]=Tgr0Re[2*iGrSet+0]*normal[3*iXi+0]+Tgr0Re[2*iGrSet+1]*normal[3*iXi+2];

        if (tgCmplx)
        {
          TXiIm[1*iGrSet+0]=TgrIm[2*iGrSet+0]*normal[3*iXi+0]+TgrIm[2*iGrSet+1]*normal[3*iXi+2];
        }
        if(calcTg0)
        {
          if (tg0Cmplx)
          {
            TXi0Im[1*iGrSet+0]=Tgr0Im[2*iGrSet+0]*normal[3*iXi+0]+Tgr0Im[2*iGrSet+1]*normal[3*iXi+2];
          }
          else
          {
            TXi0Im[1*iGrSet+0]=0.0;
          }
        }
      }
      else if (ntgComp==6) // 2D in plane
      {
        TXiRe[4*iGrSet+0]=TgrRe[6*iGrSet+0]*normal[3*iXi+0]+TgrRe[6*iGrSet+2]*normal[3*iXi+2];  // txx
        TXiRe[4*iGrSet+1]=TgrRe[6*iGrSet+2]*normal[3*iXi+0]+TgrRe[6*iGrSet+1]*normal[3*iXi+2];  // txz
        TXiRe[4*iGrSet+2]=TgrRe[6*iGrSet+3]*normal[3*iXi+0]+TgrRe[6*iGrSet+5]*normal[3*iXi+2];  // tzx
        TXiRe[4*iGrSet+3]=TgrRe[6*iGrSet+5]*normal[3*iXi+0]+TgrRe[6*iGrSet+4]*normal[3*iXi+2];  // tzz
        if (tgCmplx)
        {
          TXiIm[4*iGrSet+0]=TgrIm[6*iGrSet+0]*normal[3*iXi+0]+TgrIm[6*iGrSet+2]*normal[3*iXi+2];  // txx
          TXiIm[4*iGrSet+1]=TgrIm[6*iGrSet+2]*normal[3*iXi+0]+TgrIm[6*iGrSet+1]*normal[3*iXi+2];  // txz
          TXiIm[4*iGrSet+2]=TgrIm[6*iGrSet+3]*normal[3*iXi+0]+TgrIm[6*iGrSet+5]*normal[3*iXi+2];  // tzx
          TXiIm[4*iGrSet+3]=TgrIm[6*iGrSet+5]*normal[3*iXi+0]+TgrIm[6*iGrSet+4]*normal[3*iXi+2];  // tzz
        }
        if(calcTg0)
        {
          TXi0Re[4*iGrSet+0]=Tgr0Re[6*iGrSet+0]*normal[3*iXi+0]+Tgr0Re[6*iGrSet+2]*normal[3*iXi+2];  // txx
          TXi0Re[4*iGrSet+1]=Tgr0Re[6*iGrSet+2]*normal[3*iXi+0]+Tgr0Re[6*iGrSet+1]*normal[3*iXi+2];  // txz
          TXi0Re[4*iGrSet+2]=Tgr0Re[6*iGrSet+3]*normal[3*iXi+0]+Tgr0Re[6*iGrSet+5]*normal[3*iXi+2];  // tzx
          TXi0Re[4*iGrSet+3]=Tgr0Re[6*iGrSet+5]*normal[3*iXi+0]+Tgr0Re[6*iGrSet+4]*normal[3*iXi+2];  // tzz
          if (tg0Cmplx)
          {
            TXi0Im[4*iGrSet+0]=Tgr0Im[6*iGrSet+0]*normal[3*iXi+0]+Tgr0Im[6*iGrSet+2]*normal[3*iXi+2];  // txx
            TXi0Im[4*iGrSet+1]=Tgr0Im[6*iGrSet+2]*normal[3*iXi+0]+Tgr0Im[6*iGrSet+1]*normal[3*iXi+2];  // txz
            TXi0Im[4*iGrSet+2]=Tgr0Im[6*iGrSet+3]*normal[3*iXi+0]+Tgr0Im[6*iGrSet+5]*normal[3*iXi+2];  // tzx
            TXi0Im[4*iGrSet+3]=Tgr0Im[6*iGrSet+5]*normal[3*iXi+0]+Tgr0Im[6*iGrSet+4]*normal[3*iXi+2];  // tzz
          }
          else
          {
            TXi0Im[4*iGrSet+0]=0.0;
            TXi0Im[4*iGrSet+1]=0.0;
            TXi0Im[4*iGrSet+2]=0.0;
            TXi0Im[4*iGrSet+3]=0.0;
          }
        }
      }
      else if (ntgComp==18) // 2.5D
      {
        TXiRe[9*iGrSet+0]=TgrRe[18*iGrSet+ 0]*normal[3*iXi+0]+TgrRe[18*iGrSet+ 5]*normal[3*iXi+2];  // txx
        TXiRe[9*iGrSet+1]=TgrRe[18*iGrSet+ 3]*normal[3*iXi+0]+TgrRe[18*iGrSet+ 4]*normal[3*iXi+2];  // txy
        TXiRe[9*iGrSet+2]=TgrRe[18*iGrSet+ 5]*normal[3*iXi+0]+TgrRe[18*iGrSet+ 2]*normal[3*iXi+2];  // txz
        TXiRe[9*iGrSet+3]=TgrRe[18*iGrSet+ 6]*normal[3*iXi+0]+TgrRe[18*iGrSet+11]*normal[3*iXi+2];  // tyx
        TXiRe[9*iGrSet+4]=TgrRe[18*iGrSet+ 9]*normal[3*iXi+0]+TgrRe[18*iGrSet+10]*normal[3*iXi+2];  // tyy
        TXiRe[9*iGrSet+5]=TgrRe[18*iGrSet+11]*normal[3*iXi+0]+TgrRe[18*iGrSet+ 8]*normal[3*iXi+2];  // tyz
        TXiRe[9*iGrSet+6]=TgrRe[18*iGrSet+12]*normal[3*iXi+0]+TgrRe[18*iGrSet+17]*normal[3*iXi+2];  // tzx
        TXiRe[9*iGrSet+7]=TgrRe[18*iGrSet+15]*normal[3*iXi+0]+TgrRe[18*iGrSet+16]*normal[3*iXi+2];  // tzy
        TXiRe[9*iGrSet+8]=TgrRe[18*iGrSet+17]*normal[3*iXi+0]+TgrRe[18*iGrSet+14]*normal[3*iXi+2];  // tzz
        if (tgCmplx)
        {
          TXiIm[9*iGrSet+0]=TgrIm[18*iGrSet+ 0]*normal[3*iXi+0]+TgrIm[18*iGrSet+ 5]*normal[3*iXi+2];  // txx
          TXiIm[9*iGrSet+1]=TgrIm[18*iGrSet+ 3]*normal[3*iXi+0]+TgrIm[18*iGrSet+ 4]*normal[3*iXi+2];  // txy
          TXiIm[9*iGrSet+2]=TgrIm[18*iGrSet+ 5]*normal[3*iXi+0]+TgrIm[18*iGrSet+ 2]*normal[3*iXi+2];  // txz
          TXiIm[9*iGrSet+3]=TgrIm[18*iGrSet+ 6]*normal[3*iXi+0]+TgrIm[18*iGrSet+11]*normal[3*iXi+2];  // tyx
          TXiIm[9*iGrSet+4]=TgrIm[18*iGrSet+ 9]*normal[3*iXi+0]+TgrIm[18*iGrSet+10]*normal[3*iXi+2];  // tyy
          TXiIm[9*iGrSet+5]=TgrIm[18*iGrSet+11]*normal[3*iXi+0]+TgrIm[18*iGrSet+ 8]*normal[3*iXi+2];  // tyz
          TXiIm[9*iGrSet+6]=TgrIm[18*iGrSet+12]*normal[3*iXi+0]+TgrIm[18*iGrSet+17]*normal[3*iXi+2];  // tzx
          TXiIm[9*iGrSet+7]=TgrIm[18*iGrSet+15]*normal[3*iXi+0]+TgrIm[18*iGrSet+16]*normal[3*iXi+2];  // tzy
          TXiIm[9*iGrSet+8]=TgrIm[18*iGrSet+17]*normal[3*iXi+0]+TgrIm[18*iGrSet+14]*normal[3*iXi+2];  // tzz
        }
        if (calcTg0)
        {
          TXi0Re[9*iGrSet+0]=Tgr0Re[18*iGrSet+ 0]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+ 5]*normal[3*iXi+2];  // txx
          TXi0Re[9*iGrSet+1]=Tgr0Re[18*iGrSet+ 3]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+ 4]*normal[3*iXi+2];  // txy
          TXi0Re[9*iGrSet+2]=Tgr0Re[18*iGrSet+ 5]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+ 2]*normal[3*iXi+2];  // txz
          TXi0Re[9*iGrSet+3]=Tgr0Re[18*iGrSet+ 6]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+11]*normal[3*iXi+2];  // tyx
          TXi0Re[9*iGrSet+4]=Tgr0Re[18*iGrSet+ 9]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+10]*normal[3*iXi+2];  // tyy
          TXi0Re[9*iGrSet+5]=Tgr0Re[18*iGrSet+11]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+ 8]*normal[3*iXi+2];  // tyz
          TXi0Re[9*iGrSet+6]=Tgr0Re[18*iGrSet+12]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+17]*normal[3*iXi+2];  // tzx
          TXi0Re[9*iGrSet+7]=Tgr0Re[18*iGrSet+15]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+16]*normal[3*iXi+2];  // tzy
          TXi0Re[9*iGrSet+8]=Tgr0Re[18*iGrSet+17]*normal[3*iXi+0]+Tgr0Re[18*iGrSet+14]*normal[3*iXi+2];  // tzz
          if (tg0Cmplx)
          {
            TXi0Im[9*iGrSet+0]=Tgr0Im[18*iGrSet+ 0]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+ 5]*normal[3*iXi+2];  // txx
            TXi0Im[9*iGrSet+1]=Tgr0Im[18*iGrSet+ 3]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+ 4]*normal[3*iXi+2];  // txy
            TXi0Im[9*iGrSet+2]=Tgr0Im[18*iGrSet+ 5]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+ 2]*normal[3*iXi+2];  // txz
            TXi0Im[9*iGrSet+3]=Tgr0Im[18*iGrSet+ 6]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+11]*normal[3*iXi+2];  // tyx
            TXi0Im[9*iGrSet+4]=Tgr0Im[18*iGrSet+ 9]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+10]*normal[3*iXi+2];  // tyy
            TXi0Im[9*iGrSet+5]=Tgr0Im[18*iGrSet+11]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+ 8]*normal[3*iXi+2];  // tyz
            TXi0Im[9*iGrSet+6]=Tgr0Im[18*iGrSet+12]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+17]*normal[3*iXi+2];  // tzx
            TXi0Im[9*iGrSet+7]=Tgr0Im[18*iGrSet+15]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+16]*normal[3*iXi+2];  // tzy
            TXi0Im[9*iGrSet+8]=Tgr0Im[18*iGrSet+17]*normal[3*iXi+0]+Tgr0Im[18*iGrSet+14]*normal[3*iXi+2];  // tzz
          }
          else
          {
            TXi0Im[9*iGrSet+0]=0.0;
            TXi0Im[9*iGrSet+1]=0.0;
            TXi0Im[9*iGrSet+2]=0.0;
            TXi0Im[9*iGrSet+3]=0.0;
            TXi0Im[9*iGrSet+4]=0.0;
            TXi0Im[9*iGrSet+5]=0.0;
            TXi0Im[9*iGrSet+6]=0.0;
            TXi0Im[9*iGrSet+7]=0.0;
            TXi0Im[9*iGrSet+8]=0.0;
          }
        }
      }
    }
  }
}
