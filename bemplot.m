function scale=bemplot(varargin)
%BEMPLOT   Boundary element plotting.
%
%   s = BEMPLOT(nod,elt,typ,u,c) makes a color plot of a scalar field on the 
%   deformed boudary element mesh.
%
%   nod    Nodes.
%   elt    Elements.
%   typ    Element types.
%   u      Displacements for all degrees of freedom (nDof * 1). 
%          Default: zero displacements.
%   c      Color values at the boundary element collocation 
%          points (nCol * 1). Default: no coloring.
%   s      Deformation scale.
% 
%   The color scaling is determined by the range of c, or by the current
%   setting of CAXIS.  The scaled color values are used as indices into
%   the current COLORMAP. The shading model is set by the SHADING
%   command. Both 2D and 3D meshes are supported. For 2D plots, a 
%   hatched diagram plot is used instread of coloring.
%
%   The arguments u and c are optional:
%   BEMPLOT(nod,elt,typ) plots the undeformed mesh.
%
%   BEMPLOT(...,'KeyName','KeyValue',...) sets the value of the specified 
%   property. Multiple property values can be set with a single statement. 
%   Accepted properties are:
%
%   PARAMETER     DEFAULT   DESCRIPTION
%
%   'PlotNode'              Plot nodes. Default 'off' for 3D plots, 'on' for
%                           2D plots.
%   'PlotCol'     'on'      Plot collocation points.
%   'PlotUndef'   'off'     Plot the undeformed mesh.
%   'NumberNode'  'off'     Display node numbers.
%   'NumberCol'   'off'     Display collocation point indices.
%   'DefScale'    'auto'    Deformation scale for deformation plots.
%   'nDiv'         4        Number of element divisions for plotting.
%                           A small number results in low rendering times,
%                           a higher number results in smoother plots for
%                           curved/deformed elements and color plots.
%
%   In addition, the following properties can be set for 3D plots:
%
%   PARAMETER     DEFAULT           DESCRIPTION
%
%   'EltEdge'     'on'              Plot the element boundaries. Default 'on'.
%   'EltBack'     'off'             Color the element side with negative normal. 
%   'EltColor'    [0.90 0.90 0.90]  Element color RGB value if no coloring is used.
%   'BackColor'   [0.39 0.47 0.64]  Element back color RGB value.
%   'AutoColor'   'off'             Use the displacement norm for coloring.
%
%   The following properties can be set for 2D plots:
%
%   PARAMETER     DEFAULT           DESCRIPTION
%
%   'HatchColor'  [0.39 0.47 0.64]  Hatch color RGB value for diagram plots.
%
%
%   All unknown parameters are redirected to the PATCH function
%   that is used to plot the element surfaces for 3D plots or to the
%   LINE function  used to plot the elements in 2D plots.

% Stijn Francois
% October 2007

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

% PROCESS INPUT ARGUMENTS
if nargin<3, error('minimum 3 input arguments required.'); end
iarg=1;
while iarg<=nargin && ~isstr(varargin{iarg}), iarg=iarg+1; end
nod=varargin{1};
elt=varargin{2};
typ=varargin{3};
if iarg<5, u=[]; else, u=varargin{4}; end
if iarg<6, c=[]; else, c=varargin{5}; end
paramlist=varargin(iarg:end);
if isempty(elt), elt = zeros(0,3); end % allow for an empty elt array 

probDim=bemdimension(elt,typ);

[nDiv,paramlist]=cutparam('ndiv',4,paramlist);
[eltbound,paramlist]=cutparam('EltEdge','on',paramlist);
[PlotCol,paramlist]=cutparam('PlotCol','off',paramlist);
if probDim==2, [PlotNode,paramlist]=cutparam('PlotNode','on',paramlist); end 
if probDim==3, [PlotNode,paramlist]=cutparam('PlotNode','off',paramlist); end
[NumberCol,paramlist]=cutparam('NumberCol','off',paramlist);
[NumberNode,paramlist]=cutparam('NumberNode','off',paramlist);
[scale,paramlist]=cutparam('DefScale','auto',paramlist);
[backplot,paramlist]=cutparam('EltBack','off',paramlist);
[clearplot,paramlist]=cutparam('ClearPlot','off',paramlist);
[plotundef,paramlist]=cutparam('PlotUndef','off',paramlist);
[AutoColor,paramlist]=cutparam('AutoColor','off',paramlist);
[EltColor,paramlist]=cutparam('EltColor',[0.90 0.90 0.90],paramlist);
[BackColor,paramlist]=cutparam('BackColor',[0.39 0.47 0.64],paramlist);
[HatchColor,paramlist]=cutparam('BackColor',[0.39 0.47 0.64],paramlist);

% MESH PROPERTIES
nCol=size(bemcollpoints(nod,elt,typ),1);
nNod=size(nod,1);
[col,colType,colID]=bemcollpoints(nod,elt,typ);
for iCol=1:nCol
  if colType(iCol)==1
    [dumarg,propind]=id2prop(elt(:,1),colID(iCol));
  else
    [dumarg,propind]=id2prop(nod(:,1),colID(iCol));
  end
end

if ~isreal(u)
  warning('Imaginary parts of complex input argument ''u'' are ignored');
  u=real(u);
end
if ~isreal(c)
  warning('Imaginary parts of complex input argument ''c'' are ignored');
  c=real(c);
end
if ~isempty(c), c=reshape(c,nCol,1); end

% ELEMENT TYPES (all element typeIDs used in elt)
typeID=unique(elt(:,2));
nTypeID=length(typeID);

% COMPUTE PROBLEM DIMENSION
xlen=abs(max(nod(:,2))-min(nod(:,2)));
ylen=abs(max(nod(:,3))-min(nod(:,3)));
zlen=abs(max(nod(:,4))-min(nod(:,4)));
lenscale= sqrt(xlen^2+ylen^2+zlen^2);

% AUTO DEFORMATION SCALING
if ischar(scale)
  if probDim==3
    xdef=max(abs(u(1:3:end)));
    ydef=max(abs(u(2:3:end)));
    zdef=max(abs(u(3:3:end)));
  elseif probDim==2
    xdef=max(abs(u(1:2:end)));
    ydef=0;
    zdef=max(abs(u(2:2:end)));
  end
  defsum=xdef+ydef+zdef;
  if ~(defsum==0)
    scale=0.05*lenscale/defsum;
  else
    scale=1;
  end
end

if strcmpi(AutoColor,'on') && ~isempty(u) && probDim==3
  c=sqrt(u(1:3:end).^2+u(2:3:end).^2+u(3:3:end).^2);
end

if isempty(u), scale=NaN; end 

% PREPARE FIGURE
if ~ishold
  clf;
  if probDim==2, view(0,0); end
end

% 3D PLOTTING ROUTINE
if probDim==3
  % LOOP OVER ELEMENT TYPES
  for iTyp=1:nTypeID
    plotElt=elt(find(elt(:,2)==typeID(iTyp)),:);
    nPlotElt=size(plotElt,1);

    % ELEMENT PROPERTIES
    [parent,nEltNod,nEltCol,TypeN,TypeM,NodXi]=bemeltdef(typeID(iTyp),typ);
    if nEltCol==1
      if parent==0         % 1D element
        colXi=[0 0];
      elseif parent==1     % 2D trianglular element
        colXi=[1/3 1/3];
      elseif parent==2     % 2D quadrilateral element
        colXi=[0 0];
      end
    else
      colXi=NodXi;
    end

    [xi,face,edge]=parentgrid(nDiv,parent);
    nXi =size(xi,1);
    nfac=size(face,1);
    N=transpose(bemshape(TypeN,xi));
    Nr=transpose(bemshape(TypeN,colXi));
    M=transpose(bemshape(TypeM,xi));
    Mr=transpose(bemshape(TypeM,colXi));
    Mn=transpose(bemshape(TypeM,NodXi));
    
    ptcoord=zeros(nPlotElt*nXi,3);
    undefcoord=zeros(nPlotElt*nXi,3);
    totface=zeros(nPlotElt*nfac,size(face,2));
    colCoord=zeros(nPlotElt*nEltCol,3);
    NodCoord=zeros(nPlotElt*nEltNod,3);

    colval=zeros(nPlotElt*nXi,1);

    % BACK COLORING
    if strcmpi(backplot,'on')
      normal = bemnormal(nod,plotElt,typ,xi);
      backcoord=zeros(nPlotElt*nXi,3);
    end

    % LOOP OVER ELEMENTS
    for iPlotElt=1:nPlotElt
      [nodcoord,nodind]=id2prop(nod,plotElt(iPlotElt,3:nEltNod+2));

      ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:)=N*nodcoord;
      if strcmpi(plotundef,'on')
        undefcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:)=...
                             ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:);
      end
      colCoord((iPlotElt-1)*nEltCol+1:iPlotElt*nEltCol,:)=Nr*nodcoord;
      NodCoord((iPlotElt-1)*nEltNod+1:iPlotElt*nEltNod,:)=nodcoord;

      nEltNod=length(nodind);
      if (nEltCol==1)
        for iCol=1:nCol
          if (colType(iCol)==1 & colID(iCol)==plotElt(iPlotElt,1))
            eltColIndex=iCol;
          end
        end
      else
        eltColIndex=zeros(1,nEltNod);
        for iNod=1:nEltNod
          for iCol=1:nCol
            if (colType(iCol)==2 & colID(iCol)==nod(nodind(iNod),1))
              eltColIndex(iNod)=iCol;
            end
          end
        end
        
      end
      if ~isempty(u)
        ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:) = ...
        ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:) ...
          +  scale*M*[u(3*eltColIndex-2) u(3*eltColIndex-1) u(3*eltColIndex-0)];
        colCoord((iPlotElt-1)*nEltCol+1:iPlotElt*nEltCol,:)= ...
        colCoord((iPlotElt-1)*nEltCol+1:iPlotElt*nEltCol,:) ...
          + scale*Mr*[u(3*eltColIndex-2) u(3*eltColIndex-1) u(3*eltColIndex-0)];
        NodCoord((iPlotElt-1)*nEltNod+1:iPlotElt*nEltNod,:)= ...
        NodCoord((iPlotElt-1)*nEltNod+1:iPlotElt*nEltNod,:) ...
          + scale*Mn*[u(3*eltColIndex-2) u(3*eltColIndex-1) u(3*eltColIndex-0)];
      end
      if strcmpi(backplot,'on')
        backcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:) = ...
          ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:) ...
                                         - 1e-4*squeeze(normal(:,:,iPlotElt)).';
      end
      totface((iPlotElt-1)*nfac+1:(iPlotElt)*nfac,:)=face+(iPlotElt-1)*nXi;
      if ~isempty(c)
        colval((iPlotElt-1)*nXi+1:iPlotElt*nXi)=M*c(eltColIndex);
      end
    end

    % PLOT PATCHES
    if ~isempty(c)
      patch('Vertices',ptcoord,...
            'Faces',totface,...
            'CData',colval,...
            'LineStyle','none',...
            'FaceColor','interp',...
             paramlist{:});
    else % uncoulored plot
      if ((EltColor>=0)&(EltColor<=1))
      patch('Vertices',ptcoord,...
            'Faces',totface,...
            'FaceColor',EltColor,...
            'Linestyle','none',...
            'Facelighting','flat',...
            'AmbientStrength',0.90,...
            'DiffuseStrength',0.40);
      end
    end
    if strcmpi(backplot,'on')
      patch('Vertices',backcoord,...
            'Faces',totface,...
            'FaceColor',BackColor,...
            'Linestyle','none',...
            'Facelighting','flat',...
            'AmbientStrength',0.80,...
            'DiffuseStrength',1.00);
    end

    % PLOT ELEMENT BOUNDARIES
    if strcmpi(eltbound,'on')
      for iPlotElt=1:nPlotElt
        coordutil=ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:);
        line(coordutil(edge,1),...
             coordutil(edge,2),...
             coordutil(edge,3),...
            'Color',[0.25 0.25 0.25],...
            'linewidth',0.25);
        if strcmpi(backplot,'on')
          coordutil=backcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:);
          line(coordutil(edge,1),...
               coordutil(edge,2),...
               coordutil(edge,3),...
              'Color',[0.25 0.25 0.25],...
              'linewidth',0.25);
        end
      end
    end

    % PLOT UNDEFORMED MESH
    if strcmpi(plotundef,'on')
      for iPlotElt=1:nPlotElt
        coordutil=undefcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:);
        line(coordutil(edge,1),...
             coordutil(edge,2),...
             coordutil(edge,3),...
            'Color',[0.4 0.4 0.4],...
            'linewidth',0.25);
      end
    end

    % PLOT NODES
    if strcmpi(PlotNode,'on')
      line(NodCoord(:,1),...
           NodCoord(:,2),...
           NodCoord(:,3),...
           'Marker','o',...
           'Linestyle','none',...
           'MarkerSize',4,...
           'MarkerFaceColor','none',...
           'MarkerEdgeColor','k');
    end

    % PLOT COLLOCATION POINTS
    if strcmpi(PlotCol,'on')
      line(colCoord(:,1),...
           colCoord(:,2),...
           colCoord(:,3),...
           'Marker','o',...
           'Linestyle','none',...
           'MarkerSize',4,...
           'MarkerFaceColor','k',...
           'MarkerEdgeColor','none');
    end
    if strcmpi(NumberCol,'on')
      text(colCoord(:,1)+lenscale/50,...
           colCoord(:,2)+lenscale/50,...
           colCoord(:,3)+lenscale/50,...
           num2str([1:nCol]'),...
           'VerticalAlignment','Bottom',...
           'HorizontalAlignment','left',...
           'Clipping','off');
    end
    
    if strcmpi(NumberNode,'on')
      text(nod(:,2)+lenscale/50,...
           nod(:,3)+lenscale/50,...
           nod(:,4)+lenscale/50,...
           num2str(nod(:,1)),...
           'VerticalAlignment','Bottom',...
           'HorizontalAlignment','left',...
           'Clipping','off');
    end
  end

%===============================================================================
% 2D PLOTTING ROUTINE
elseif (probDim==2)
  if ~isempty(c)
    maxconval=max(abs(c));
    if (maxconval>0)
      cScale=lenscale/maxconval/5;
    else
      cScale=lenscale/50;
    end
  end

  % LOOP OVER ELEMENT TYPES
  for iTyp=1:nTypeID
    plotElt=elt(find(elt(:,2)==typeID(iTyp)),:);
    nPlotElt=size(plotElt,1);

    % ELEMENT PROPERTIES
    [parent,nEltNod,nEltCol,TypeN,TypeM,NodXi]=bemeltdef(typeID(iTyp),typ);
    if nEltCol==1
      colXi=0;
    else
      colXi=NodXi(:,1);
    end

    % PARENT GRID
    xi=linspace(-1,1,nDiv+1)';
    nXi=nDiv+1;
    normal = bemnormal(nod,plotElt,typ,[xi zeros(nDiv+1,1)]);

    N=transpose(bemshape(TypeN,xi));
    Nr=transpose(bemshape(TypeN,colXi));
    M=transpose(bemshape(TypeM,xi));
    Mr=transpose(bemshape(TypeM,colXi));
    Mn=transpose(bemshape(TypeM,NodXi(:,1)));

    ptcoord=zeros(nPlotElt*nXi,3);
    undefcoord=zeros(nPlotElt*nXi,3);
    colCoord=zeros(nPlotElt*nEltCol,3);
    NodCoord=zeros(nPlotElt*nEltNod,3);

    colval=zeros(nPlotElt*nXi,1);
    
    % LOOP OVER ELEMENTS
    for iPlotElt=1:nPlotElt
      [nodcoord,nodind]=id2prop(nod,plotElt(iPlotElt,3:nEltNod+2));
      nEltNod=length(nodind);
      if (nEltCol==1)
        for iCol=1:nCol
          if (colType(iCol)==1 & colID(iCol)==plotElt(iPlotElt,1))
            eltColIndex=iCol;
          end
        end
      else
        eltColIndex=zeros(1,nEltNod);
        for iNod=1:nEltNod
          for iCol=1:nCol
            if (colType(iCol)==2 & colID(iCol)==nodind(iNod))
              eltColIndex(iNod)=iCol;
            end
          end
        end
      end

      ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:)=N*nodcoord;
      undefcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:)...
                                    =ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:);

      colCoord((iPlotElt-1)*nEltCol+1:iPlotElt*nEltCol,:)=Nr*nodcoord;
      NodCoord((iPlotElt-1)*nEltNod+1:iPlotElt*nEltNod,:)=nodcoord;

      if ~isempty(u)
        ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:) = ...
        ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,:) ...
         + scale*M*[u(2*eltColIndex-1) 0*u(2*eltColIndex-0) u(2*eltColIndex-0)];
        colCoord((iPlotElt-1)*nEltCol+1:iPlotElt*nEltCol,:)= ...
        colCoord((iPlotElt-1)*nEltCol+1:iPlotElt*nEltCol,:)...
         +scale*Mr*[u(2*eltColIndex-1) 0*u(2*eltColIndex-0) u(2*eltColIndex-0)];
        NodCoord((iPlotElt-1)*nEltNod+1:iPlotElt*nEltNod,:)= ...
        NodCoord((iPlotElt-1)*nEltNod+1:iPlotElt*nEltNod,:)...
        + scale*Mn*[u(2*eltColIndex-1) 0*u(2*eltColIndex-0) u(2*eltColIndex-0)];
      end

      line(ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,1),...
           ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,2),...
           ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,3),...
           'Color','k',...
           'linewidth',1.50,...
           'clipping','off',...
           paramlist{:});

      if ~isempty(c)
        colutil=M*c(eltColIndex);
        normx=normal(1,:,iPlotElt)';
        normz=normal(3,:,iPlotElt)';
        line(ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,1)+cScale*colutil.*normx,...
             ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,2),...
             ptcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,3)+cScale*colutil.*normz,...
            'Color',HatchColor,...
            'linewidth',0.25,...
            'clipping','off');
        for iXi=1:nXi
          line([ptcoord((iPlotElt-1)*nXi+iXi,1);ptcoord((iPlotElt-1)*nXi+iXi,1)+cScale*colutil(iXi)*normx(iXi)],...
               [ptcoord((iPlotElt-1)*nXi+iXi,2);ptcoord((iPlotElt-1)*nXi+iXi,2)],...
               [ptcoord((iPlotElt-1)*nXi+iXi,3);ptcoord((iPlotElt-1)*nXi+iXi,3)+cScale*colutil(iXi)*normz(iXi)],...
               'Color',HatchColor,...
               'linewidth',0.25,...
               'clipping','off');
        end
      end
    end
    
    % PLOT UNDEFORMED MESH
    if strcmpi(plotundef,'on')
      for iPlotElt=1:nPlotElt
        line(undefcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,1),...
             undefcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,2),...
             undefcoord((iPlotElt-1)*nXi+1:iPlotElt*nXi,3),...
             'Color',[0.50 0.50 0.50],...
             'linewidth',0.50);
      end
    end

    % PLOT NODES
    if strcmpi(PlotNode,'on')
      line(NodCoord(:,1),...
           NodCoord(:,2),...
           NodCoord(:,3),...
          'Marker','o',...
          'Linestyle','none',...
          'MarkerSize',5,...
          'linewidth',1.00,...
          'MarkerFaceColor','none',...
          'MarkerEdgeColor','k',...
          'clipping','off');
    end

    % PLOT COLLOCATION POINTS
    if strcmpi(PlotCol,'on')
      line(colCoord(:,1),...
           colCoord(:,2),...
           colCoord(:,3),...
          'Marker','o',...
          'Linestyle','none',...
          'MarkerSize',5,...
          'linewidth',0.05,...
          'MarkerFaceColor','k',...
          'MarkerEdgeColor','none',...
          'clipping','off');
    end

    if strcmpi(NumberCol,'on')
      text(colCoord(:,1)+lenscale/50,...
           colCoord(:,2),...
           colCoord(:,3)+lenscale/50,...
           num2str([1:nCol]'),...
           'VerticalAlignment','Bottom',...
           'HorizontalAlignment','left',...
           'Clipping','off');
    end
  end
  
  if strcmpi(NumberNode,'on')
    text(nod(:,2)+lenscale/50,...
         nod(:,3),...
         nod(:,4)+lenscale/50,...
         num2str(nod(:,1)),...
         'VerticalAlignment','Bottom',...
         'HorizontalAlignment','left',...
         'Clipping','off');
  end
end

axis off;
axis equal;

% RETURN OUTPUT ARGUMENTS ONLY IF REQUESTED
if nargout<1, clear('scale'); end

%===============================================================================
function [node,face,edge]=parentgrid(ndiv,parent)
%PARENTGRID Patch Grid of the parent element (3D plots).
%   [NOD,FACE,EDGE]=PARENTGRID(NDIV,PARENT) computes the coordinates NOD,
%   the faces FACE and the edge EDGE of the patch for the parent element
%   in the natural coordinate system (xi1,xi2). A distinction is made
%   between PARENT=1 for a triangular parent element and PARENT=2
%   for a quadrilateral parent element. The resulting matrices can be passed
%   to a patch plot.

np = ndiv+1;
% LINE ELEMENT
if (parent==0)
  xi=0:1/ndiv:1;
  node=xi;
  face=[];
  edge=[];

% TRIANGULAR PARENT ELEMENT
elseif (parent==1)
  xi1=0:1/ndiv:1;
  xi2=1:-1/ndiv:0;
  % NODES
  node=zeros((np)*(np+1)/2,2);
  iutil=1;
  for ixi2=1:np
    for ixi1=1:ixi2
      node(iutil,1) =xi1(ixi1);
      node(iutil,2) =xi2(ixi2);
      iutil=1+iutil;
    end
  end
  % FACES
  face=zeros(ndiv^2,3);
  uprow=1;
  eltutil=1;
  for idiv=1:ndiv
    nuprow = length(uprow);
    downrow = uprow(end)+1:uprow(end)+nuprow+1;
     % upright triangles
      for ipt=1:nuprow
        face(eltutil,:) = [downrow(ipt) downrow(ipt+1) uprow(ipt)];
        eltutil = eltutil+1;
      end
     % downright triangles
      for ipt=1:nuprow-1
        face(eltutil,:) = [downrow(ipt+1) uprow(ipt+1) uprow(ipt)];
        eltutil = eltutil+1;
      end
    uprow = downrow;
  end
  % EDGES
  edge=zeros(1,3*ndiv+1);
  edge(1)=1;
  iedge=2;
  util=1;
  for idiv=1:ndiv
    edge(iedge)=edge(iedge-1)+util;
    util=util+1;
    iedge=iedge+1;
  end
  edge(iedge:iedge+ndiv-1)=edge(iedge-1)+1:edge(iedge-1)+ndiv;
  iedge=iedge+ndiv;
  for idiv=1:ndiv
    edge(iedge)=edge(iedge-1)-util;
    util=util-1;
    iedge=iedge+1;
  end

% RECTANGULAR PARENT ELEMENT
elseif (parent==2)
  xi1=-1:2/ndiv:1;
  xi2=-1:2/ndiv:1;
  % NODES
  node=zeros(np^2,2);
  node(:,1)=repmat(xi1.',np,1);
  node(:,2)=reshape(repmat(xi2,np,1),np^2,1);
  % FACES
  face=zeros(ndiv^2,4);
  for idiv=1:ndiv
    for jdiv=1:ndiv
      face(idiv+ndiv*(jdiv-1),1)= (ndiv+1)*(jdiv-1)+idiv;
      face(idiv+ndiv*(jdiv-1),2)= (ndiv+1)*(jdiv-1)+idiv+1;
      face(idiv+ndiv*(jdiv-1),3)= (ndiv+1)*(jdiv)  +idiv+1;
      face(idiv+ndiv*(jdiv-1),4)= (ndiv+1)*(jdiv)  +idiv;
    end
  end
  % EDGES
  edge=[1:np-1 np:np:np^2 np^2-1:-1:np^2-np+1 np^2-2*np+1:-np:1];
else
  error('Parent element undefined')
end

%===============================================================================
function [prop,propind]=id2prop(item,itemID)
%ID2PROP item properties from ID.
%   PROP=ID2PROP(ITEM,ITEMID)

nID=numel(itemID);
propind=[];
for iID=1:nID
  ind=find(item(:,1)==itemID(iID),1);
  propind=[propind,ind];
end
prop=item(propind,2:end);

%===============================================================================
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
