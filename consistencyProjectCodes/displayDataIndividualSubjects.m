function displayDataIndividualSubjects(pos,dTFPowerDBAllSubjects,dPowerVsFreqAllSubjects,timeValsTF,freqValsTF,freqVals)

xLims = [-0.5 1.3]; cLims = [-1.5 3]; yLims = [0 100]; dPSDLims = [-3 6];

numGroups = size(dTFPowerDBAllSubjects,1);
numSubjects = length(pos);
hPlotsTF = getPlotHandles(numSubjects,numGroups,[0.1 0.05 0.4 0.9],0.01,0.01,1);
hPlotsDPSD = getPlotHandles(numSubjects,1,[0.55 0.05 0.15 0.9],0.005,0.005,1);

for i=1:numSubjects
    for j=1:numGroups
        pcolor(hPlotsTF(i,j),timeValsTF,freqValsTF,squeeze(dTFPowerDBAllSubjects(j,pos(i),:,:))'); 
        shading(hPlotsTF(i,j),'interp');
        if (i==numSubjects) && (j==1)
            set(hPlotsTF(i,j),'Tickdir','out','XTick',[0 0.8],'YTick',[0 50 100],'XTicklabel',[0 0.8],'YTicklabel',[0 50 100]);
            xlabel(hPlotsTF(i,j),'Time (s)'); ylabel(hPlotsTF(i,j),'Frequency (Hz)'); 
        elseif (i==numSubjects) && (j==numGroups)
            set(hPlotsTF(i,j),'Tickdir','out','XTick',[0 0.8],'YTick',[0 50 100],'XTicklabel',[0 0.8],'YTicklabel',[]);
            xlabel(hPlotsTF(i,j),'Time (s)');
        else
            set(hPlotsTF(i,j),'Tickdir','out','XTick',[0 0.8],'YTick',[0 50 100],'XTicklabel',[],'YTicklabel',[]);
        end
        caxis(hPlotsTF(i,j),cLims);
        axis(hPlotsTF(i,j),[xLims yLims]);
        
        plot(hPlotsDPSD(i),freqVals,10*squeeze(dPowerVsFreqAllSubjects(j,pos(i),:)));
        hold(hPlotsDPSD(i),'on');
        if (i==numSubjects)
            set(hPlotsDPSD(i),'Tickdir','out','XTick',[0 50 100],'YTick',[-3 0 3 6],'XTicklabel',[0 50 100],'YTicklabel',[-3 0 3 6]);
            xlabel(hPlotsDPSD(i),'Frequency (Hz)');
        else
            set(hPlotsDPSD(i),'Tickdir','out','XTick',[0 50 100],'YTick',[-3 0 3 6],'XTicklabel',[],'YTicklabel',[]);
        end
        axis(hPlotsDPSD(i),[yLims dPSDLims]);
    end
end
colormap magma;