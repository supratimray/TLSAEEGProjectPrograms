% cellArray must be 1xN cell; each vector of the cell must be a matrix of size Mx1 
% This function returns an MxN matrix

function newMatrix = concatenateCellArrayToMatrix(cellArray)
    cols = size(cellArray,2);
    cellElementRows = cellfun(@length,cellArray);
    numRowsElement = unique(cellElementRows);
    discordantElementCol = [];
    if length(numRowsElement)>1 
        numRowsElement = mode(cellElementRows);
        discordantElementCol = find(cellElementRows ~= numRowsElement);        
        disp(['Discrepency in no. of elements in cells at ',num2str(discordantElementCol),' column(s).']);
    end
    
    newMatrix = zeros(numRowsElement,cols);    
    for iCol = 1:cols
        clear vector
        vector = cellArray{1,iCol};
        if ~isempty(discordantElementCol); if ismember(iCol,discordantElementCol); vector = resizeVector(vector,numRowsElement); end; end;
        newMatrix(:,iCol) = vector;
    end
end

function newVector = resizeVector(oldVector,Size)
    if isrow(oldVector); oldVector = oldVector'; end
    if length(oldVector)<Size
        newVector = cat(1,oldVector,(oldVector(end)+zeros((Size-length(oldVector)),1)));
    else
        newVector = oldVector(1:Size);
    end
end
