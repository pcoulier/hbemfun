function rec=bemgrid_cyl(r,theta,z)
%BEMGRID_CYL  Receiver grid in cylindrical coordinates.
%
%   rec=BEMGRID_CYL(r,theta,z) constructs the matrix with receiver point
%   coordinates of a grid in cylindrical coordinates for use with the BEMXFER 
%   command.
%   
%   r      r-coordinates of the grid points (nrRec * 1) or (1 * nrRec).
%   theta  theta-coordinates of the grid points (ntRec * 1) or (1 * ntRec).
%   z      z-coordinates of the grid points (nzRec * 1) or (1 * nzRec).
%   rec    Coordinates of the receiver points (nRec * 3).

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

nrRec=numel(r);
ntRec=numel(theta);
nzRec=numel(z);

nRec=nrRec*ntRec*nzRec;
rec=zeros(nRec,3);

iRec=1;
for izRec=1:nzRec
  for itRec=1:ntRec
    for irRec=1:nrRec
      rec(iRec,1)=r(irRec)*cos(theta(itRec));
      rec(iRec,2)=r(irRec)*sin(theta(itRec));
      rec(iRec,3)=z(izRec);
      iRec=iRec+1;
    end
  end
end
