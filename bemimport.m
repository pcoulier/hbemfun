function [nod,elt] = bemimport(importformat,varargin)
%BEMIMPORT   Boundary element mesh importing.
%
%   [nod,elt] = BEMIMPORT('ansys',nodfile,eltfile,typ,typID)
%   imports a boundary element mesh from Ansys. The elements used
%   in ANSYS should have the same nodal connection as the element type
%   "typID".
%
%   nodfile  File with node definitions, generated in Ansys with the
%            Ansys-command NWRITE.
%   eltfile  File with element definitions, generated in Ansys with the
%            Ansys-command EWRITE.
%   typ      Element types.
%   typID    Element type ID.
%   nod      Nodes.
%   elt      Elements.

% CHECK BEMFUN LICENSE
bemfunlicense('VerifyOnce');

switch lower(importformat)
case 'ansys'
  nodalfile = varargin{1};
  eltfile   = varargin{2};
  typ       = varargin{3};
  typID     = varargin{4};

  [Parent,nNod,nCol,TypeN,TypeM,NodDef,EltDim]=bemeltdef(typID,typ);

  % IMPORT NODAL DATA
  tmpfile = [tempname '.tmp'];
  fid = fopen(nodalfile,'rt');
  fid2 = fopen(tmpfile,'wt');
  fseek(fid,0,'eof');
  endpos=ftell(fid);
  fseek(fid,0,'bof');

  % ADD SPACES IN FRONT OF MINUS SIGNS
  % IN THE NODAL FILE
  prevchar  = ' ';
  while ftell(fid) < endpos
    currchar = fscanf(fid,'%c',1);
    if (currchar == '-') & ~(prevchar=='E')
      currchar = ' -';
    end
    fprintf(fid2,'%c',currchar);
    prevchar=currchar;
  end
  fclose(fid);
  fclose(fid2);

  % READ NODAL DATA
  fid2 = fopen(tmpfile,'rt');
  nodaldata = textscan(fid2,'%n %n %n %n','emptyValue',0);
  
  if (EltDim==2)
    nod=[nodaldata{1} nodaldata{2} nodaldata{3} nodaldata{4}];
  elseif (EltDim==1)
    % in this case, y and z-coordinates are swapped...
    nod=[nodaldata{1} nodaldata{2} nodaldata{4} nodaldata{3}];
  end
  
  fclose(fid2);
  delete(tmpfile);

  % READ ELEMENT DATA
  fid=fopen(eltfile,'rt');
  eltdata=textscan(fid,'%n %n %n %n %n %n %n %n %n %n %n %n %n %n','emptyValue',0);
  fclose(fid); 
  
  nElt=size(eltdata{1},1);
  elt = [eltdata{14} repmat([typID],nElt,1) eltdata{1:nNod}];
  
otherwise
    error('Import format unknown');
end
