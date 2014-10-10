%BEMNORMAL   Element normals.
%
%   normal = BEMNORMAL(nod,elt,typ,xi) computes the normal in the specified
%   points for all elements of the boundary element mesh.
%
%   nod    Nodes.
%   elt    Elements.
%   typ    Element types.
%   Xi     Natural coordinates of the sampling points (nXi * 2).
%          For 1D elements, the second column is not used.
%   normal Element normals (3 * nXi * nElt).
