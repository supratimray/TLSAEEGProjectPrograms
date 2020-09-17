function displayDataIndividualSubjects(pos,dTFPowerDBAllSubjects,dPowerVsFreqAllSubjects,timeValsTF,freqValsTF,freqVals)

xLims = [-0.5 1.3]; cLims = [-1 3]; yLims = [0 100];

numGroups = size(dTFPowerDBAllSubjects,1);
numSubjects = length(pos);
hPlotsTF = getPlotHandles(numSubjects,numGroups,[0.05 0.05 0.6 0.9],0.005,0.005,1);
hPlotsDPSD = getPlotHandles(numSubjects,1,[0.75 0.05 0.15 0.9],0.005,0.005,1);

for i=1:numSubjects
    for j=1:numGroups
        pcolor(hPlotsTF(i,j),timeValsTF,freqValsTF,squeeze(dTFPowerDBAllSubjects(j,pos(i),:,:))'); 
        shading(hPlotsTF(i,j),'interp');
        caxis(hPlotsTF(i,j),cLims);
        axis(hPlotsTF(i,j),[xLims yLims]);
        
        plot(hPlotsDPSD(i),freqVals,10*squeeze(dPowerVsFreqAllSubjects(j,pos(i),:)));
        hold(hPlotsDPSD(i),'on');
    end
end