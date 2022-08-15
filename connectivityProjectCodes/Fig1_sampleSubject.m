function Fig1_sampleSubject(dataForDisplay,sample_SubjIndx,chanlocs,elecGroups,freqRanges,figpath)
% usage: Fig1_sampleSubject(dataForDisplayAllGroups{2,2},subjIndx,freqRanges{2}) (agegroup2, SG)
%%%%%% defining electrode clusters %%%%%
elecCluster = defineElecClusters();
colorsCluster = {'#EF4444','#EF4444','#009F75','#009F75','#FAA31B','#FAA31B','#88C6ED','#88C6ED'};
ClusterInfo.elecCluster = elecCluster;
ClusterInfo.colorsCluster = colorsCluster;

freqVals = dataForDisplay.freqVals;
cLims = [0 1]; % for unipolar setting
noseDir = '+X';

%%%% defining figure layout %%%%%%
figR = figure('numbertitle', 'off','name','Figure 1: Example subject');
figR.PaperType = 'a4';
figR.PaperUnits = 'centimeters';
figR.PaperSize = [12.3 24.7]; % Nature specifications
figR.PaperOrientation = 'Landscape';
figR.PaperPosition = [0 0 figR.PaperSize];
figR.Color = [1 1 1]; % White background

plotTopo{1} = getPlotHandles(2,1,[0.02 0.1 0.13 0.8],0.01,0.01,1);
plotTopo{2} = getPlotHandles(2,1,[0.24 0.1 0.15 0.8],0.01,0.01,1);
plotTopo{3} = getPlotHandles(2,1,[0.40 0.1 0.19 0.8],0.01,0.01,1);
plotTopo{4} = getPlotHandles(2,3,[0.65 0.1 0.3 0.8],0.01,0.01,1);

%%%%%% making power topoplots %%%%%
% absolute stimulus period power
subplot(plotTopo{1}(1)); 
powTopos = dataForDisplay.stPowerTopoAllSubjects{sample_SubjIndx};
topoplot_murty(powTopos,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir',noseDir,'emarker',{'.','k',12,1},'emarkercolors',powTopos);

cbr = colorbar;
cbr.Position = cbr.Position + [0.065 0 0 0.2];
cbr.FontSize = 13;
cbr.Ticks = [-2 0 1];
cbr.TickLabels = {'-2', '0', '1'};
caxis([-2 1]);

% change in power across stimulus and baseline
subplot(plotTopo{1}(2)); 
powTopos = dataForDisplay.diffPowerTopoAllSubjects{sample_SubjIndx};
topoplot_murty(powTopos,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir',noseDir,'emarker',{'.','k',12,1},'emarkercolors',powTopos);

cbr = colorbar;
cbr.Position = cbr.Position + [0.065 0 0 0.2];
cbr.FontSize = 13;
cbr.Ticks = [-4 0 4];
cbr.TickLabels = {'-4', '0', '4'};
caxis([-4 4]);

%%%%%%%% PSD with frequency  %%%%
for elecClusterSide = 1:2 % limited to left and right side
    subplot(plotTopo{2}(elecClusterSide));
    BLpsd = dataForDisplay.logBLPowerVsFreqAllSubjects{sample_SubjIndx,elecClusterSide};
    STpsd = dataForDisplay.logSTPowerVsFreqAllSubjects{sample_SubjIndx,elecClusterSide};
    hold on;
    line([freqRanges(1) freqRanges(1)],[-20 20],'LineStyle','--','Color','k');
    line([freqRanges(2) freqRanges(2)],[-20 20],'LineStyle','--','Color','k');
    set(plotTopo{2}(elecClusterSide),'TickDir','out','TickLength',[0.02, 0.01]);
    plot(freqVals,10*(STpsd-BLpsd),'LineWidth',1.5,'Color','k');
    line([0 100],[0 0],'LineStyle',':','Color','k');
    ax = gca;
    ax.YColor = 'k';
    ylim([-4 6]);
    if(elecClusterSide==1)
        set(plotTopo{2}(elecClusterSide),'XTick',[0 50 100],'YTick',[-2:2:6],'XTicklabel',[],'YTicklabel',[],'FontSize',11);
    else % elecClusterSide==2
        set(plotTopo{2}(elecClusterSide),'XTick',[0 50 100],'YTick',[-2:2:6],'XTicklabel',[0 50 100],'YTicklabel',[-2:2:6],'FontSize',11);
        ylabel(plotTopo{2}(elecClusterSide),'Change in power (dB)','FontSize',11); %,'Color','b');
        xlabel(plotTopo{2}(elecClusterSide),'Frequency (Hz)','FontSize',11);
    end
end

%%%%%%%%%%%%%%%%%% Conn. topos  %%%%%%%%%%%%%%%%%%
topoData = cell(2,3);
for elecClusterSide = 1:2
        axS = subplot(plotTopo{3}(elecClusterSide));
        ConnData = dataForDisplay.connFreqBandsAllSubjects;
        ppcConn = ConnData{sample_SubjIndx,elecClusterSide};
        outConn = squeeze(nanmean(ppcConn,2));
        topoplot(outConn,chanlocs,'plotrad',0.5,'emarker',{'.','k',10,1}); hold on;         
        topoplot([],chanlocs,'hcolor','none','plotrad',0.5,'plotchans',elecGroups{elecClusterSide},'emarker',{'o','w',6,2}); % highlighting reference electrodes
        for k = 1:size(ppcConn,2) % across "electrodes"
            topoData{elecClusterSide,k} = squeeze(ppcConn(1,k,:));
        end
        if(elecClusterSide==1) 
            cbr = colorbar(axS);
            cbr.Position = cbr.Position + [0.065 0 0 0];
            cbr.FontSize = 13;
            cbr.Ticks = [0 1];
            cbr.TickLabels = {'0', '1'};
        end
        caxis(cLims);
end

% Conn. with distance & color coding
for elecClusterSide = 1:2
    for electrodeIndx = 1:3
        subplot(plotTopo{4}(elecClusterSide,electrodeIndx));
        plotConnData_sample(topoData{elecClusterSide,electrodeIndx},elecGroups{elecClusterSide}(electrodeIndx),ClusterInfo,chanlocs);
        set(plotTopo{4}(elecClusterSide,electrodeIndx),'XTick',[-1 0 1],'YTick',[cLims(1):(diff(cLims)/2):cLims(2)],'XTicklabel',[],'YTicklabel',[],'FontSize',13);
        set(plotTopo{4}(elecClusterSide,electrodeIndx),'Xdir','reverse');
        ylim(cLims);
        set(gca,'TickDir','out','TickLength',[0.03, 0.01]);
        if(elecClusterSide==2 && electrodeIndx==1)
            set(plotTopo{4}(elecClusterSide,electrodeIndx),'XTick',[-1 0 1],'YTick',[cLims(1):(diff(cLims)/2):cLims(2)],'XTicklabel',[-1 0 1],'YTicklabel',[cLims(1):(diff(cLims)/2):cLims(2)],'FontSize',13);
            xlabel(plotTopo{4}(elecClusterSide,electrodeIndx),'cos( \Delta \theta)','FontSize',13,'Interpreter','tex'); 
            ytickformat('%.1f');
        else
            set(plotTopo{4}(elecClusterSide,electrodeIndx),'XTick',[-1 0 1],'YTick',[cLims(1):(diff(cLims)/2):cLims(2)],'XTicklabel',[],'YTicklabel',[],'FontSize',13);
        end
    end
end

print(figR,'-painters',fullfile(figpath,'Fig1_sample_ppc_SG'),'-dtiff','-r300');
print(figR,'-painters',fullfile(figpath,'Fig1_sample_ppc_SG'),'-dsvg','-r300');

% NOTE: following code generates electrodes with color coding, as per clusters defined
% hf = figure;
% makeRefelecTopo(chanlocs,elecCluster,colorsCluster);
% print(hf,'-painters',fullfile(figpath,'topo_categories'),'-dsvg','-r300');
end

function plotConnData_sample(Gconn,elecRef,ClusterInfo,chanlocs)
loc = getElecLocAngles(chanlocs);
tot = 1:64;
% Fitting Conn vs. dist (Binned version of this is in Murty's dual gamma paper)
dist = sqrt((loc.azi(elecRef)-loc.azi(tot)).^2+(loc.ele(elecRef)-loc.ele(tot)).^2);
cos_dist = cos((dist/180)*pi);
for k = 1:length(ClusterInfo.elecCluster)
    if(mod(k,2))
        scatter(cos_dist(ClusterInfo.elecCluster{k}),Gconn(ClusterInfo.elecCluster{k}),20,'o','MarkerEdgeColor',ClusterInfo.colorsCluster{k},'MarkerFaceColor',ClusterInfo.colorsCluster{k}); hold on;
    else
        scatter(cos_dist(ClusterInfo.elecCluster{k}),Gconn(ClusterInfo.elecCluster{k}),20,'^','MarkerEdgeColor',ClusterInfo.colorsCluster{k},'MarkerFaceColor','w'); hold on;
    end
end 
end
function makeRefelecTopo(chanlocs,elecCluster,colorsCluster) % used only once to generate electrodes with color code
    for elecC=1:8
        if(mod(elecC,2))
            topoplot([],chanlocs,'hcolor','none','plotchans',elecCluster{elecC},'emarker',{'.',colorsCluster{elecC},12,1});
        else
            topoplot([],chanlocs,'hcolor','none','plotchans',elecCluster{elecC},'emarker',{'^',colorsCluster{elecC},4,1});
        end
    end
end
function loc = getElecLocAngles(chanlocs)
azi = zeros(1,length(chanlocs)); ele = azi;
for e = 1:length(chanlocs)
    azi(e) = chanlocs(e).sph_theta;
    ele(e) = chanlocs(e).sph_phi;
end
loc.azi = azi;
loc.ele = ele;
end
function elecCluster = defineElecClusters()
% Lcentral
elecCluster{1}= [8 9 13 18 19 32+[11 15 16 20]];
% Rcentral
elecCluster{2}= [10 11 15 20 21 32+[12 17 18 22]];
% Ltemporal
elecCluster{3}= [12 17 32+[9 10 19]];
% Rtemporal 
elecCluster{4}= [16 22 32+[13 14 23]];
% Lvisual 
elecCluster{5}= [23 24 28 29 32+[24 25 28 29]];
% Rvisual 
elecCluster{6}= [26 27 31 32 32+[26 27 31 32]];
% Lfrontal 
elecCluster{7}= [1 3 4 32+[1 2 5 6]];
% Rfrontal 
elecCluster{8}= [2 6 7 32+[3 4 7 8]]; 
end