function [hPlot,hPlotPos] = getPlotHandles_Columns(plotPos,plotParts,plotGap)    
    if ~exist('plotGap','var') || isempty(plotGap); plotGap = 0.01; end;
    plPartsR = plotParts;
%     plotCount = length(plotParts);
    for iPlot = 1:length(plotParts)
        plotWidth(iPlot) = (plotPos(3)-plotGap*(length(plotParts)-1))*(plPartsR(iPlot)/sum(plPartsR)); %#ok<*AGROW>
        if iPlot>1; plotStart(iPlot) = plotStart(iPlot-1)+plotWidth(iPlot-1)+plotGap; else plotStart(iPlot) = plotPos(1); end
        hPlot(iPlot) = subplot('position',[plotStart(iPlot) plotPos(2) plotWidth(iPlot) plotPos(4)]);
        hPlotPos{iPlot} = [plotStart(iPlot) plotPos(2) plotWidth(iPlot) plotPos(4)];
%         plotCount = iPlot-1;
    end
end