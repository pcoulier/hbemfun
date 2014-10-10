%BEMELTDEF   Boundary element properties.
%
%   [map,nNod,nCol,nID,mID,nodDef,eltDim,axiSym,nGauss,nEltDiv]
%                                                         =BEMELTDEF(typID,typ)
%   returns various element properties.
%
%   typID   Element type ID.
%   typ     Element types.
%   map     Parent element mapping (for 2D elements only). 1 for a triangular
%           mapping, 2 for a rectangular mapping. 0 for 1D elements
%   nNod    Number of element nodes.
%   nCol    Number of element collocation points.
%   nID     Shape function ID for element Geometry (see BEMSHAPE).
%   mID     Boundary element interpolation function ID (see BEMSHAPE).
%   nodDef  Natural coordinates of the nodes (Nnod * 2).
%   eltDim  Element dimension. 1 for 1D element, 2 for 2D element.
%   axiSym  1 for an axisymmetric geometry, 0 otherwise.
%   nGauss  Number of Gaussian points.
%   nEltDiv Number of element divisions for integration.
