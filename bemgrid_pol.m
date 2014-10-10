function rec=bemgrid_pol(r,theta)
%BEMGRID_POL  Receiver grid in polar coordinates.
%
%   rec=BEMGRID_POL(r,theta) constructs the matrix with receiver point
%   coordinates of a grid in polar coordinates for use with the BEMXFER 
%   command. The polar coordinate system is defined for 2D problems
%   in the (x,z)-plane, where a positive theta angle corresponds to
%   a rotation from the z-axis to the x-axis.
%   
%   r      r-coordinates of the grid points (nr * 1) or (1 * nr).
%   theta  theta-coordinates of the grid points (ntheta * 1) or (1 * ntheta).
%   rec    Coordinates of the receiver points (nRec * 3).

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

nrRec=numel(r);
ntRec=numel(theta);

nRec=nrRec*ntRec;
rec=zeros(nRec,3);

iRec=1;
for itRec=1:ntRec
  for irRec=1:nrRec
    rec(iRec,1)=r(irRec)*sin(theta(itRec));
    rec(iRec,2)=0;
    rec(iRec,3)=r(irRec)*cos(theta(itRec));
    iRec=iRec+1;
  end
end
