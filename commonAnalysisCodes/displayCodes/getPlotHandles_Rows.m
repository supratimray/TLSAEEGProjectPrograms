function [hPlot,hPlotPos] = getPlotHandles_Rows(plotPos,plotParts,plotGap)    
    if ~exist('plotGap','var') || isempty(plotGap); plotGap = 0.01; end;
    plPartsR = fliplr(plotParts);
    plotCount = length(plotParts);
    for iPlot = 1:length(plotParts)
        plotHeight(iPlot) = plotPos(4)*(plPartsR(iPlot)/sum(plPartsR))-plotGap; %#ok<*AGROW>
        if iPlot>1; plotStart(iPlot) = plotStart(iPlot-1)+plotHeight(iPlot-1)+2*plotGap; else plotStart(iPlot) = plotPos(2); end
        hPlot(plotCount) = subplot('position',[plotPos(1) plotStart(iPlot) plotPos(3) plotHeight(iPlot)]);
        hPlotPos{plotCount} = [plotPos(1) plotStart(iPlot) plotPos(3) plotHeight(iPlot)];
        plotCount = plotCount-1;
    end
end