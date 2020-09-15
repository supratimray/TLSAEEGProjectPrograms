function [imData,xlims,ylims,clims] = pasteAsImage(spHandle,xlims,ylims,clims,colormapName)

    randNum = randi(100000);
    axesUnits = spHandle.Units;
    axesPosition = spHandle.Position;
    axesPosition(3) = axesPosition(3)*2;
    figPC = figure(randNum);
    figPC.PaperUnits = axesUnits;
    figPC.PaperPosition = axesPosition;
    figPC.Visible = 'off';
    pcHandle = copyobj(spHandle,figPC);

    if ~exist('colormapName','var') || isempty(colormapName); colormapName = 'magma'; end;
    colormap(colormapName)
    
    if ~exist('xlims','var') || isempty(xlims); xlims = xlim; end
    if ~exist('ylims','var') || isempty(ylims); ylims = ylim; end
    if ~exist('clims','var') || isempty(clims); clims = caxis; end
    
    axis normal
    xlim(pcHandle,xlims);
    ylim(pcHandle,ylims);
    caxis(pcHandle,clims);
    
    pcHandle.TickDir = 'in';
    pcHandle.YTick = [];
    pcHandle.XTick = [];
    axis off
    set(pcHandle,'position',[0 0 1 1],'units','normalized','box','off');
    print(figPC,['pasteAsImageTempFig' num2str(randNum)],'-dbmp','-r600','-opengl');
    close(figPC);
    imData = imread(['pasteAsImageTempFig' num2str(randNum) '.bmp']);
    delete(['pasteAsImageTempFig' num2str(randNum) '.bmp']);
    
    subplot(spHandle); cla;
    imagesc(xlims,fliplr(ylims),imData); set(spHandle,'ydir','normal','tickdir','out');
    xlim(xlims); ylim(ylims);
       
end