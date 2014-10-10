function Ft = timflexmat(Fomg,omg,delt,nFlex)
%TIMFLEXMAT   Hybrid frequency-time domain flexibility matrices.
%
%  Ft = TIMFLEXMAT(Fomg,omg,delt,nFlex) computes the time domain flexibility
%  influence matrices from the flexibility in the frequency domain.
%  
%  Fomg  Frequency  domain flexibility (nDof * nDof * nFreq).
%  omg   Frequency sampling (1 * nFreq) or (nFreq * 1).
%  delt  Time step.
%  nFlex Number of time domain flexibility matrices.
%  Ft    Time domain influence matrices (nDof * nDof * nFlex).
%

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

nFreq=length(omg);
nDof=size(Fomg,1);
Ft=zeros(nDof,nDof,nFlex);

ind0=find(~omg);
ind1=find(omg);

shape=zeros(nFreq,1);

%%%% INTEGRATION OF F1
shape(ind0)=delt/2;
shape(ind1)=(1-exp(i*delt*omg(ind1))+i*delt*omg(ind1))./(delt*omg(ind1).^2);



for iDof=1:nDof
  for jDof=1:nDof
    Ft(iDof,jDof,1)=1/pi*(trapz(omg.',real(squeeze(Fomg(iDof,jDof,:)).*shape)));
  end
end

%%%% INTEGRATION OF F2->FN
shape(ind0)=delt;
shape(ind1)=(2-2*cos(delt*omg(ind1))) ./ (delt*omg(ind1).^2);
if (nFlex>1)
  for iDof=1:nDof
    for jDof=1:nDof
      Ft(iDof,jDof,2:end)=1/pi*real(filon(omg,[i*delt*(1:(nFlex-1))],(shape.*squeeze(Fomg(iDof,jDof,:))).','exp','linear'));
    end
  end
end
