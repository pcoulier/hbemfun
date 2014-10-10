%BEMCOLLPOINTS   Boundary element collocation points.
%
%   [col,colTyp,ID] = BEMCOLLPOINTS(nod,elt,typ) returns the boundary
%   element collocation points.
%
%   nod    Nodes.
%   elt    Elements.
%   typ    Element types.
%   col    Collocation point coordinates (nCol * 3).
%   colTyp Collocation type (nCol * 1). 1 for a centroid collocation, 2 for a
%          nodal collocation point.
%   ID     The corresponding element number (if colTyp=1) or node number
%          (if colTyp=2).
