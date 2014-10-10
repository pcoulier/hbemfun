function [ds,cs,h]=bemgridplot_pol(varargin)
%BEMGRIDPLOT_POL  Plot results on a Polar receiver grid.
%
%   BEMGRIDPLOT_POL(r,theta,y,u,c) creates a color plot of a scalar wave field 
%   c(r,t,y) on a deformed mesh.  The mesh displacements are given by the vector 
%   field u that has been computed with the BEMXFER command for a receiver grid 
%   defined with the BEMGRID_POL command.  The grid is defined in the polar 
%   coordinate system and the deformation is defined in Cartesian coordinates.
%
%   r       Receiver grid r-coordinates (nr * 1).
%   theta   Receiver grid theta-coordinates (nt * 1).
%   y       Receiver grid y-coordinates (ny * 1).
%   u       Receiver grid displacements (nRecDof * 1). This vector contains
%           the displacementsof the grid points in Cartesian coordinates
%   c       Scalar wave field (nRec * 1). If the vector c is not specified,
%           the scalar wave field is derived from the displacement field
%           using the optional 'ScalarField' parameter.
%
%   BEMGRIDPLOT_POL(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%
%   'DefScale'     Deformation scale.  Default: 'auto'.
%   'ColorScale'   2-element vector specifying the color scale.  Default: 'auto'.
%   'ScalarField'  Scalar wave field used if input argument c is not specified.
%                  either 'norm', 'ur', 'ut', 'ux', 'uy', 'uz'. Default: 'norm'.
%
%   Additional parameters are redirected to the SURF function.
%
%   [ds,cs,h] = BEMGRIDPLOT_POL(...) returns the deformation scale ds, the color
%   scale cs, and a handle h to the SURF object(s).

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

% DETERMINE INPUT ARGUMENTS
rRec=varargin{1};
tRec=varargin{2};
yRec=varargin{3};
iarg=1;
while iarg<=nargin && ~ischar(varargin{iarg}), iarg=iarg+1; end
narg=iarg-1;
switch narg
case {4,5}
  uRec=varargin{4};
  if narg==5, cRec=varargin{5}; else cRec=[]; end
otherwise
  error('Incorrect number of input arguments.');
end
paramlist=varargin(narg+1:end);
[scalarField,paramlist]=cutparam('ScalarField','norm',paramlist);

if ~isreal(uRec)
  warning('Ignoring imaginary part of input argument uRec');
  uRec=real(uRec);
end
nrRec=length(rRec);
ntRec=length(tRec);
nyRec=length(yRec);
nRec=nrRec*ntRec*nyRec;

if ~isempty(uRec)
  if size(uRec,1)==nRec
    ux=0;
    uy=reshape(uRec(1:1:nRec),[nrRec ntRec nyRec]);
    uz=0;
  elseif size(uRec,1)==2*nRec
    ux=reshape(uRec(1:2:2*nRec),[nrRec ntRec nyRec]);
    uy=0;
    uz=reshape(uRec(2:2:2*nRec),[nrRec ntRec nyRec]);
  elseif size(uRec,1)==3*nRec
    ux=reshape(uRec(1:3:3*nRec),[nrRec ntRec nyRec]);
    uy=reshape(uRec(2:3:3*nRec),[nrRec ntRec nyRec]);
    uz=reshape(uRec(3:3:3*nRec),[nrRec ntRec nyRec]);
  else
    error('First dimension of uRec should correspond with the number of degrees of freedom');
  end
else
  ux=0;
  uy=0;
  uz=0;
end

theta=repmat(tRec(:).',[nrRec 1 nzRec]);
ur=ux.*sin(theta)+uz.*cos(theta);
ut=ux.*cos(theta)-uz.*sin(theta);

if ~isempty(cRec)
  c=permute(reshape(cRec(1:nRec),[ntRec nrRec nyRec]),[3 2 1]);
else
  switch lower(scalarField)
  case 'norm'
    c=sqrt(ur.^2+ut.^2+uy.^2);
  case 'ur'
    c=ur;
  case 'ut'
    c=ut;
  case 'ux'
    c=ux;
  case 'uy'
    c=uy;
  case 'uz'
    c=uz;
  otherwise
    error('Unknown value for key option ''scalarField''');
  end
end

[ds,cs,h]=waveplot_cyl(rRec,pi/2-tRec,yRec,ur,ut,uy,c,paramlist{:});
set(gca,'YDir','normal');
set(gca,'ZDir','normal');

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout<1, clear('ds'); end
if nargout<2, clear('cs'); end
if nargout<3, clear('h'); end

%===============================================================================
function [ds,cs,h]=waveplot_cyl(varargin) %% MODIFIED...

%WAVEPLOT_CYL   Plot a wave field in cylindrical coordinates.
%
%   WAVEPLOT_CYL(r,theta,z,ur,ut,uz) creates a color plot of a scalar wave field
%   c(r,theta,z) on a deformed mesh.  The mesh displacements are given by the
%   vector field (ur,ut,uz).  The mesh and the deformation are defined in
%   cylindrical coordinates.
%
%   WAVEPLOT_CYL(r,theta,z,u,c) is an alternative syntax for a 2-D wave field
%   u = (ur,uz) or 3-D wave field u = (ur,ut,uz).
%
%   r       Vertex r-coordinates (n1 * n2 * n3) or (n1 * 1).
%   theta   Vertex theta-coordinates (n1 * n2 * n3) or (n2 * 1).
%   z       Vertex z-coordinates (n1 * n2 * n3) or (n3 * 1).
%   ur      Mesh displacements in r-direction (n1 * n2 * n3).
%   ut      Mesh displacements in theta-direction (n1 * n2 * n3).
%   uz      Mesh displacements in z-direction (n1 * n2 * n3).
%   u       Mesh displacements (2 * n1 * n2 * n3) or (3 * n1 * n2 * n3).
%   c       Scalar wave field (n1 * n2 * n3).  Default: SQRT(ur^2+ut^2+uz^2).
%
%   If one or more of the dimensions of ur, ut, uz, or c is equal to 1 instead
%   of n1, n2, or n3, this function attempts to apply REPMAT to ur, ut, uz, or c
%   in order to generate an (n1 * n2 * n3) matrix.
%
%   WAVEPLOT_CYL(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%
%   'DefScale'     Deformation scale.  Default: 'auto'.
%   'ColorScale'   2-element vector specifying the color scale.  Default: 'auto'.
%
%   Additional parameters are redirected to the SURF function.
%
%   [ds,cs,h] = WAVEPLOT_CYL(...) returns the deformation scale ds, the color
%   scale cs, and a handle h to the SURF object(s).

% Mattias Schevenels
% March 2009

% DETERMINE INPUT ARGUMENTS
r=varargin{1};
theta=varargin{2};
z=varargin{3};
iarg=1;
while iarg<=nargin && ~ischar(varargin{iarg}), iarg=iarg+1; end
narg=iarg-1;
switch narg
case {4,5}
  u=varargin{4};
  if isempty(u) || (numel(u)==1 && u==0), u=[0 0 0]; end
  u=squeeze(u);
  if size(u,1)==1, u=u.'; end
  udim=[size(u) 1];
  switch udim(1)
  case 2
    ur=reshape(u(1,:),udim(2:end));
    ut=0;
    uz=reshape(u(2,:),udim(2:end));
  case 3
    ur=reshape(u(1,:),udim(2:end));
    ut=reshape(u(2,:),udim(2:end));
    uz=reshape(u(3,:),udim(2:end));
  otherwise
    error('The first non-singleton dimension of the wave field u must be 2 or 3.');
  end
  if narg==5, c=varargin{5}; else c=sqrt(ur.^2+ut.^2+uz.^2); end
case {6,7}
  ur=varargin{4};
  ut=varargin{5};
  uz=varargin{6};
  if narg==7, c=varargin{7}; else c=sqrt(ur.^2+ut.^2+uz.^2); end
otherwise
  error('Incorrect number of input arguments.');
end
if isempty(r), r=0; end
if isempty(theta), theta=0; end
if isempty(z), z=0; end
if isempty(ur), ur=0; end
if isempty(ut), ut=0; end
if isempty(uz), uz=0; end
if isempty(c), c=nan; end
paramlist=varargin(narg+1:end);

% TRANSFORM R,THETA,Z TO 3D MATRICES AND DETERMINE MESH SIZE
r=squeeze(r);
theta=squeeze(theta);
z=squeeze(z);
rvec=(ndims(r)==2) && (min(size(r))==1);
tvec=(ndims(theta)==2) && (min(size(theta))==1);
zvec=(ndims(z)==2) && (min(size(z))==1);
if ~(all([rvec tvec zvec]) || all(~[rvec tvec zvec]))
  error('''r'', ''theta'', and ''z'' must be either 3 vectors or 3 equally sized 3D matrices.');
end
if rvec
  n1=length(r);
  n2=length(theta);
  n3=length(z);
  r=repmat(reshape(r,[n1 1 1]),[1 n2 n3]);
  theta=repmat(reshape(theta,[1 n2 1]),[n1 1 n3]);
  z=repmat(reshape(z,[1 1 n3]),[n1 n2 1]);
else
  n1=size(r,1);
  n2=size(r,2);
  n3=size(r,3);
end
if length(find([n1 n2 n3]>1))<2, error('Unable to plot 1-D objects.'); end

% CHECK SIZE OF R,THETA,Z
rdim=size(r);
tdim=size(theta);
zdim=size(z);
if length(rdim)~=length(tdim) || any(rdim~=tdim), error('Incompatible input arguments: check size of ''r'' and ''theta''.'); end
if length(rdim)~=length(zdim) || any(rdim~=zdim), error('Incompatible input arguments: check size of ''r'' and ''z''.'); end

% REPMAT U,V,W,C TO OBTAIN 3D MATRICES
try ur=xrepmat(ur,[n1 n2 n3]); catch error('Incompatible input arguments: check size of ''ur''.'); end
try ut=xrepmat(ut,[n1 n2 n3]); catch error('Incompatible input arguments: check size of ''ut''.'); end
try uz=xrepmat(uz,[n1 n2 n3]); catch error('Incompatible input arguments: check size of ''uz''.'); end
try c=xrepmat(c,[n1 n2 n3]); catch error('Incompatible input arguments: check size of ''c''.'); end

% CHANGE COORDINATE SYSTEM
C=cos(theta(:));
S=sin(theta(:));
[x,y]=pol2cart(theta,r);
uxrec=ur(:).*C-ut(:).*S;
uyrec=ur(:).*S+ut(:).*C;
uzrec=uz(:);
uxrec=reshape(uxrec,[n1 n2 n3]);
uyrec=reshape(uyrec,[n1 n2 n3]);
uzrec=reshape(uzrec,[n1 n2 n3]);

% CALL WAVEPLOT_REC
[ds,cs,h]=waveplot_rec(x,z,y,uxrec,uzrec,uyrec,c,paramlist{:});

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout<1, clear('ds'); end
if nargout<2, clear('cs'); end
if nargout<3, clear('h'); end

%-------------------------------------------------------------------------------

% RESHAPE/REPMAT MATRIX TO SPECIFIED SIZE
function x=xrepmat(x,tdim);
% x       Matrix to reshape.
% tdim    Target dimensions.
ntdim=length(tdim);
x=squeeze(x);
if ndims(x)==2 && min(size(x))==1, x=x(:); end
xdim=size(x);
xdim=xdim(xdim>1);
nxdim=length(xdim);
order=zeros(ntdim,1);
for ixdim=1:nxdim
  itdim=find(tdim==xdim(ixdim),1);
  if isempty(itdim), error('Dimensions do not match.'); end
  tdim(itdim)=nan;
  order(itdim)=ixdim;
end
rep=tdim;
rep(order~=0)=1;
order(order==0)=nxdim+1:ntdim;
x=permute(x,order);
x=repmat(x,rep);


%===============================================================================
function [ds,cs,h]=waveplot_rec(varargin)

%WAVEPLOT_REC   Plot a wave field in Cartesian coordinates.
%
%   WAVEPLOT_REC(x,y,z,ux,uy,uz,c) creates a color plot of a scalar wave field
%   c(x,y,z) on a deformed mesh.  The mesh displacements are given by the vector
%   field (ux,uy,uz).  The mesh and the deformation are defined in Cartesian
%   coordinates.
%
%   WAVEPLOT_REC(x,y,z,u,c) is an alternative syntax for a 2-D wave field
%   u = (ux,uz) or 3-D wave field u = (ux,uy,uz).
%
%   x       Vertex x-coordinates (n1 * n2 * n3) or (n1 * 1).
%   y       Vertex y-coordinates (n1 * n2 * n3) or (n2 * 1).
%   z       Vertex z-coordinates (n1 * n2 * n3) or (n3 * 1).
%   ux      Mesh displacements in x-direction (n1 * n2 * n3).
%   uy      Mesh displacements in y-direction (n1 * n2 * n3).
%   uz      Mesh displacements in z-direction (n1 * n2 * n3).
%   u       Mesh displacements (2 * n1 * n2 * n3) or (3 * n1 * n2 * n3).
%   c       Scalar wave field (n1 * n2 * n3).  Default: SQRT(ux^2+uy^2+uz^2).
%
%   If one or more of the dimensions of ux, uy, uz, or c is equal to 1 instead
%   of n1, n2, or n3, this function attempts to apply REPMAT to ux, uy, uz, or c
%   in order to generate an (n1 * n2 * n3) matrix.
%
%   WAVEPLOT_REC(...,ParamName,ParamValue) sets the value of the specified
%   parameters.  The following parameters can be specified:
%
%   'DefScale'     Deformation scale.  Default: 'auto'.
%   'ColorScale'   2-element vector specifying the color scale.  Default: 'auto'.
%
%   Additional parameters are redirected to the SURF function.
%
%   [ds,cs,h] = WAVEPLOT_REC(...) returns the deformation scale ds, the color
%   scale cs, and a handle h to the SURF object(s).

% Mattias Schevenels
% March 2009

% DETERMINE INPUT ARGUMENTS
x=varargin{1};
y=varargin{2};
z=varargin{3};
iarg=1;
while iarg<=nargin && ~ischar(varargin{iarg}), iarg=iarg+1; end
narg=iarg-1;
switch narg
case {4,5}
  u=varargin{4};
  if isempty(u) || (numel(u)==1 && u==0), u=[0 0 0]; end
  u=squeeze(u);
  if size(u,1)==1, u=u.'; end
  udim=[size(u) 1];
  switch udim(1)
  case 2
    ux=reshape(u(1,:),udim(2:end));
    uy=0;
    uz=reshape(u(2,:),udim(2:end));
  case 3
    ux=reshape(u(1,:),udim(2:end));
    uy=reshape(u(2,:),udim(2:end));
    uz=reshape(u(3,:),udim(2:end));
  otherwise
    error('The first non-singleton dimension of the wave field u must be 2 or 3.');
  end
  if narg==5, c=varargin{5}; else c=sqrt(ux.^2+uy.^2+uz.^2); end
case {6,7}
  ux=varargin{4};
  uy=varargin{5};
  uz=varargin{6};
  if narg==7, c=varargin{7}; else c=sqrt(ux.^2+uy.^2+uz.^2); end
otherwise
  error('Incorrect number of input arguments.');
end
if isempty(x), x=0; end
if isempty(y), y=0; end
if isempty(z), z=0; end
if isempty(ux), ux=0; end
if isempty(uy), uy=0; end
if isempty(uz), uz=0; end
if isempty(c), c=nan; end
paramlist=varargin(narg+1:end);
[ds,paramlist]=cutparam('DefScale','auto',paramlist);
[cs,paramlist]=cutparam('ColorScale','auto',paramlist);
FaceColorSpecified=any(strcmpi(paramlist,'FaceColor'));

% REMOVE SPURIOUS IMAGINARY PART FROM UX,UY,UZ IF PRESENT
urmax=max([max(abs(real(ux(:)))),max(abs(real(uy(:)))),max(abs(real(uz(:))))]);
uimax=max([max(abs(imag(ux(:)))),max(abs(imag(uy(:)))),max(abs(imag(uz(:))))]);
if uimax/urmax>1e-8
  error('Unable to plot a complex-valued wave field.');
else
  ux=real(ux);
  uy=real(uy);
  uz=real(uz);
end

% TRANSFORM X,Y,Z TO 3D MATRICES AND DETERMINE MESH SIZE
x=squeeze(x);
y=squeeze(y);
z=squeeze(z);
xvec=(ndims(x)==2) && (min(size(x))==1);
yvec=(ndims(y)==2) && (min(size(y))==1);
zvec=(ndims(z)==2) && (min(size(z))==1);
if ~(all([xvec yvec zvec]) || all(~[xvec yvec zvec]))
  error('''x'', ''y'', and ''z'' must be either 3 vectors or 3 equally sized 3D matrices.');
end
if xvec
  n1=length(x);
  n2=length(y);
  n3=length(z);
  x=repmat(reshape(x,[n1 1 1]),[1 n2 n3]);
  y=repmat(reshape(y,[1 n2 1]),[n1 1 n3]);
  z=repmat(reshape(z,[1 1 n3]),[n1 n2 1]);
else
  n1=size(x,1);
  n2=size(x,2);
  n3=size(x,3);
end
if length(find([n1 n2 n3]>1))<2, error('Unable to plot 1-D objects.'); end

% CHECK SIZE OF X,Y,Z
xdim=size(x);
ydim=size(y);
zdim=size(z);
if length(xdim)~=length(ydim) || any(xdim~=ydim), error('Incompatible input arguments: check size of ''x'' and ''y''.'); end
if length(xdim)~=length(zdim) || any(xdim~=zdim), error('Incompatible input arguments: check size of ''x'' and ''z''.'); end

% REPMAT UX,UY,UZ,C TO OBTAIN 3D MATRICES
try ux=xrepmat(ux,[n1 n2 n3]); catch error('Incompatible input arguments: check size of ''ux''.'); end
try uy=xrepmat(uy,[n1 n2 n3]); catch error('Incompatible input arguments: check size of ''uy''.'); end
try uz=xrepmat(uz,[n1 n2 n3]); catch error('Incompatible input arguments: check size of ''uz''.'); end
try c=xrepmat(c,[n1 n2 n3]); catch error('Incompatible input arguments: check size of ''c''.'); end

% DETERMINE AUTO SCALES USING SOME OBSCURE HEURISTIC RULE
if strcmpi(ds,'auto') || strcmpi(cs,'auto')
  Lx=max(x(:))-min(x(:));
  Ly=max(y(:))-min(y(:));
  Lz=max(z(:))-min(z(:));
  L=sqrt(Lx.^2+Ly.^2+Lz.^2);
  if strcmpi(ds,'auto')
    U=sqrt(ux(:).^2+uy(:).^2+uz(:).^2);
    [dum,imax]=max(U);
    d=sqrt((x(:)-x(imax)).^2+(y(:)-y(imax)).^2+(z(:)-z(imax)).^2);
    R=L/8;
    U=U(d>R);
    Umax=max(U);
    if Umax~=0
      ds=0.02*L/Umax;
    else
      ds=1;
    end
  end
  if strcmpi(cs,'auto')
    C=abs(c(:));
    [dum,imax]=max(C);
    d=sqrt((x(:)-x(imax)).^2+(y(:)-y(imax)).^2+(z(:)-z(imax)).^2);
    R=L/8;
    C=C(d>R);
    Cmax=max(C);
    if min(c(:))<0
      Cmin=-Cmax;
    else
      Cmin=0;
    end
    if Cmax~=0
      cs=[Cmin Cmax];
    else
      cs=nan;
    end
  end
end

% PREPARE FIGURE
nextplot=get(gca,'NextPlot');
if ~strcmpi(nextplot,'add')
  if strcmpi(nextplot,'replace')
    cla('reset');
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'YDir','reverse');
    set(gca,'ZDir','reverse');
    if all(x(~isnan(x))==0)
      view(-90,0);
      box('on');
    elseif all(y(~isnan(y))==0)
      view(0,0);
      box('on');
    elseif all(z(~isnan(z))==0)
      view(0,90);
      box('on');
    else
      view(142.5,30);
      box('off');
    end
  else
    cla;
  end
end
hold('on');

% PLOT SURFACES
if any([n1 n2 n3]==1)
  nn=[n1 n2 n3];
  mm=nn(nn>1);
  h=surf(reshape(x+ds*ux,mm),reshape(y+ds*uy,mm),reshape(z+ds*uz,mm),reshape(c,mm),paramlist{:});
else
  vertices=[];
  faces=[];
  cdata=[];
  h=zeros(6,1);
  for iface=1:6
    switch iface
    case 1, x1=x(1,:,:);   y1=y(1,:,:);   z1=z(1,:,:);   ux1=ux(1,:,:);   uy1=uy(1,:,:);   uz1=uz(1,:,:);   c1=c(1,:,:);
    case 2, x1=x(end,:,:); y1=y(end,:,:); z1=z(end,:,:); ux1=ux(end,:,:); uy1=uy(end,:,:); uz1=uz(end,:,:); c1=c(end,:,:);
    case 3, x1=x(:,1,:);   y1=y(:,1,:);   z1=z(:,1,:);   ux1=ux(:,1,:);   uy1=uy(:,1,:);   uz1=uz(:,1,:);   c1=c(:,1,:);
    case 4, x1=x(:,end,:); y1=y(:,end,:); z1=z(:,end,:); ux1=ux(:,end,:); uy1=uy(:,end,:); uz1=uz(:,end,:); c1=c(:,end,:);
    case 5, x1=x(:,:,1);   y1=y(:,:,1);   z1=z(:,:,1);   ux1=ux(:,:,1);   uy1=uy(:,:,1);   uz1=uz(:,:,1);   c1=c(:,:,1);
    case 6, x1=x(:,:,end); y1=y(:,:,end); z1=z(:,:,end); ux1=ux(:,:,end); uy1=uy(:,:,end); uz1=uz(:,:,end); c1=c(:,:,end);
    end
    nn=size(x1);
    mm=nn(nn>1);
    h(iface)=surf(reshape(x1+ds*ux1,mm),reshape(y1+ds*uy1,mm),reshape(z1+ds*uz1,mm),reshape(c1,mm),paramlist{:});
    hold('on');
  end
end
if ~FaceColorSpecified
  shading('interp');
end
if ~isnan(cs)
  caxis(cs);
end

% RESET NEXTPLOT STATE
set(gca,'NextPlot',nextplot);

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout<1, clear('ds'); end
if nargout<2, clear('cs'); end
if nargout<3, clear('h'); end

%-------------------------------------------------------------------------------

% CUT PARAMETER FROM LIST
function [value,paramlist]=cutparam(name,default,paramlist);
value=default;
for iarg=length(paramlist)-1:-1:1
  if strcmpi(name,paramlist{iarg})
    value=paramlist{iarg+1};
    paramlist=paramlist([1:iarg-1 iarg+2:end]);
    break
  end
end
