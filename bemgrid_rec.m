function rec=bemgrid_rec(x,y,z)
%BEMGRID_REC  Receiver grid in Cartesian coordinates.
%
%   rec=BEMGRID_REC(x,y,z) constructs the matrix with receiver point
%   coordinates of a rectangular grid for use with the BEMXFER command.
%   
%   x    x-coordinates of the grid points (nx * 1) or (1 * nx).
%   y    y-coordinates of the grid points (ny * 1) or (1 * ny).
%   z    z-coordinates of the grid points (nz * 1) or (1 * nz).
%   rec  Coordinates of the receiver points (nRec * 3).

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

nXrec=numel(x);
nYrec=numel(y);
nZrec=numel(z);

nRec=nXrec*nYrec*nZrec;
rec=zeros(nRec,3);

iRec=1;
for iZrec=1:nZrec
  for iYrec=1:nYrec
    for iXrec=1:nXrec
      rec(iRec,1)=x(iXrec);
      rec(iRec,2)=y(iYrec);
      rec(iRec,3)=z(iZrec);
      iRec=iRec+1;
    end
  end
end
