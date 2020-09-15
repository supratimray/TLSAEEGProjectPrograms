% MD: 06/05/2016
function newMatrix = removeDimIfSingleton(oldMatrix,dim)
    z = size(oldMatrix);
    if nargin<2; dim = 1; end
    if z(dim) == 1
        z(dim) = [];
        newMatrix = reshape(oldMatrix,[z 1]);
    else
        newMatrix = oldMatrix;
    end
end