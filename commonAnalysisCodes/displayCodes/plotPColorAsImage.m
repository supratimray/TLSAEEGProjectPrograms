function [imData,xlims,ylims,clims] = plotPColorAsImage(spHandle,x,y,C,xlims,ylims,clims,colormapName)

    randNum = randi(100000);
    axesUnits = spHandle.Units;
    axesPosition = spHandle.Position;
    axesPosition(3) = axesPosition(3)*2;
    figPC = figure(randNum);
    figPC.PaperUnits = axesUnits;
    figPC.PaperPosition = axesPosition;
    figPC.Visible = 'off';
    pcHandle = gca; 
    pcolor(pcHandle,x,y,C);
    shading(pcHandle,'interp'); 
    if ~exist('colormapName','var') || isempty(colormapName); colormapName = 'magma'; end;
    colormap(colormapName)
    
    if ~exist('xlims','var') || isempty(xlims); xlims = xlim; end
    if ~exist('ylims','var') || isempty(ylims); ylims = ylim; end
    if ~exist('clims','var') || isempty(clims); clims = caxis; end
    
    xlim(pcHandle,xlims);
    ylim(pcHandle,ylims);
    caxis(pcHandle,clims);
    
    pcHandle.TickDir = 'in';
    pcHandle.YTickLabel = [];
    pcHandle.XTickLabel = [];
    set(pcHandle,'position',[0 0 1 1],'units','normalized','box','off');
    print(figPC,['pcolorToImageTempFig' num2str(randNum)],'-dbmp','-r600','-opengl');
    close(figPC);
    imData = imread(['pcolorToImageTempFig' num2str(randNum) '.bmp']);
    delete(['pcolorToImageTempFig' num2str(randNum) '.bmp']);
end