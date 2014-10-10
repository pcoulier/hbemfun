function [phi]=bemrigid(col,x0)
%BEMRIGID   Rigid body modes.
%
%   [phi]=bemrigid(col,x0) computes the six threedimensional rigid body
%   modes with respect to the center of rotation x0. The first three modes 
%   correspond to a unit translation in the x,y and z-direction, respectively. 
%   The last three modes are positive unit rotations around the x, y and z-axis 
%   respectively. 
%
%   col   Coordinate array (nCol * 3). Usually corresponds to the collocation
%         point coordinates.
%   x0    Center of rotation (3 * 1) or (1 * 3). Default: [0 0 0].
%   phi   Rigid body modes  (nDof * 6).

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

if nargin<2, x0=[0 0 0]; end

nCol=size(col,1);
nDof=3*nCol;

phi=zeros(nDof,6);

for iCol=1:nCol
  x=col(iCol,1)-x0(1);
  y=col(iCol,2)-x0(2);
  z=col(iCol,3)-x0(3);
  phi(3*(iCol-1)+1:3*(iCol-1)+3,:)=[1  0  0   0  z -y;
                                    0  1  0  -z  0  x;
                                    0  0  1   y -x  0];
end

