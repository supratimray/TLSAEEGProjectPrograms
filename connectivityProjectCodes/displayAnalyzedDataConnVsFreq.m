% This program displays the data saved in the analyzedData folder

% subjectNameList is a cell array, with each cell containing a list of
% subjects. For example, subjectNameList{1} could contain all MCI/AD subjects
% while subjectNameList{2} could contain all their controls.

function displayAnalyzedDataConnVsFreq(folderSourceString,subjectNameLists,strList,projectName,refType,protocolType,stRange,freqRanges,freqRangeNames,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData,useMedianFlagData,useBLConnData)
if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('freqRanges','var')
    freqRanges{1} = [8 12]; freqRangeNames{1} = 'Alpha';
    freqRanges{2} = [20 34]; freqRangeNames{2} = 'Slow gamma';
    freqRanges{3} = [36 66]; freqRangeNames{3} = 'Fast gamma';
end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end
if ~exist('useCleanData','var');    useCleanData=0;                     end
if ~exist('useMedianFlagData','var');  useMedianFlagData = 1;           end

connMethod = 'ppc';
numGroups = length(subjectNameLists);
numFreqRanges = length(freqRanges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataForDisplayAllGroups = cell(numGroups,numFreqRanges);

for i=1:numFreqRanges % analysis done separately for each frequency
    disp(['Working on Freq: ' num2str(freqRanges{i})]);
    dataForDisplayAllGroups{1,i} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1},projectName,refType,protocolType,stRange,freqRanges(i),connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData,useBLConnData);
    dataForDisplayAllGroups{2,i} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{2},projectName,refType,protocolType,stRange,freqRanges(i),connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData,useBLConnData);
end

%%%%%%%%%%%%%%%%% Setting parameters for Display script %%%%%%%%%%%%%%%%%%%
elecGroupsCell = getElectrodeList('actiCap64',refType,0,1); % has electrode divisions
elecClusterSide = 1; % 1 for left, 2 for right, 3 for back
refElectrodes = cell2mat(elecGroupsCell{elecClusterSide+1});
chanlocs = getMontageDetails(refType);
[electrodeGroupList,groupNameList] = getElectrodeGroupsConn(refElectrodes,chanlocs);
numElectrodeGroups = length(electrodeGroupList);

%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10; displaySettings.tickLengthMedium = [0.025 0];

%%%%%%%%%%%%%%%%%%%%%%%%%% Display Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotTopo{1} = getPlotHandles(numFreqRanges,numGroups,[0.05 0.1 0.20 0.8],0.01,0.01,1); % power topoplots
plotTopo{2} = getPlotHandles(numFreqRanges,numGroups,[0.30 0.1 0.20 0.8],0.01,0.01,1); % conn topoplots
plotTopo{3} = getPlotHandles(numFreqRanges,1,[0.55 0.1 0.075 0.8],0.01,0.01,1); % conn-distance profile
plotTopo{4} = getPlotHandles(2,numElectrodeGroups/2,[0.65 0.1 0.325 0.8],0,0.02,1);

for freqBand = 1:numFreqRanges

    % 1. change in power topoplots
    for iGrp = 1:numGroups
        axS = subplot(plotTopo{1}(freqBand,iGrp));
        selpowTopos = dataForDisplayAllGroups{iGrp,freqBand}.diffPowerTopoAllSubjects;
        if(useMedianFlagData)
            avgpowTopos = median(cell2mat(selpowTopos'),2,'omitnan'); % Plots in Murty's study were with this setting
        else
            avgpowTopos = mean(cell2mat(selpowTopos'),2,'omitnan');
        end
        topoplot_murty(avgpowTopos,chanlocs,'electrodes','off','style','blank','drawaxis','off','emarker',{'.','k',12,1},'emarkercolors',avgpowTopos);
        clim([-1 1]);
        topoplot([],chanlocs,'hcolor','none','plotchans',refElectrodes,'emarker',{'o','k',10,1});
        if(iGrp == 1)
            ht = text(-0.65,0,freqRangeNames{freqBand},'HorizontalAlignment','center');
            set(ht,'Rotation',90);
        end
        if(freqBand == 1)
            title(strList{1},'FontSize',displaySettings.fontSizeLarge,'FontWeight','normal');
            if(iGrp == numGroups)
                cbr = colorbar(axS);
                cbr.Position = cbr.Position + [0.07 -0.02 0.0015 0.07];
                cbr.FontSize = displaySettings.fontSizeLarge;
                title(strList{2},'FontSize',displaySettings.fontSizeLarge,'FontWeight','normal');
            end
        end
    end

    % 2. conn topoplots
    tbinnedMeanConn = cell(numGroups,1);

    for iGrp = 1:numGroups
        axS = subplot(plotTopo{2}(freqBand,iGrp));
        connData = dataForDisplayAllGroups{iGrp,freqBand}.connFreqBandsAllSubjects(:,elecClusterSide);

        if(useMedianFlagData)
            outConn = squeeze(median(cell2mat(connData'),2,'omitnan'));
        else
            outConn = squeeze(mean(cell2mat(connData'),2,'omitnan'));
        end

        tbinnedMeanConn{iGrp} = getMeanConn(connData,electrodeGroupList);
        
        topoplot(outConn,chanlocs,'emarker',{'.','k',10,1}); hold on;
        topoplot_santosh(outConn,chanlocs,'style','contour','plotrad',0.5,'headrad',0,'ccolor','r','numcontour',3); % for drawing 25% contour
        topoplot([],chanlocs,'hcolor','none','plotrad',0.5,'headrad',0,'plotchans',refElectrodes,'emarker',{'o','k',8,1}); % highlighting reference electrodes
        clim([0 1]);
        if(freqBand == 1)
            if(iGrp==1)
                title(strList{1},'FontSize',displaySettings.fontSizeLarge,'FontWeight','normal');
            else
                title(strList{2},'FontSize',displaySettings.fontSizeLarge,'FontWeight','normal');
                cbr = colorbar(axS);
                cbr.Position = cbr.Position + [0.07 -0.02 0.0015 0.07];
                cbr.FontSize = displaySettings.fontSizeLarge;
            end
        end
    end

    displaySettings.fontSizeLarge = 10;
    displaySettings.tickLengthMedium = [0.025 0];
    displaySettings.colorNames(1,:) = [0.8 0 0.8];      % Purple
    displaySettings.colorNames(2,:) = [0.25 0.41 0.88]; % Cyan
    % 3. conn-distance profile
    displayAndcompareData(plotTopo{3}(freqBand),tbinnedMeanConn,1:numElectrodeGroups,displaySettings,[0 1],1,useMedianFlagData,1);
end

% ConnVsFreq
meanConnVsFreq = cell(numGroups,1);

for iGrp=1:numGroups
    connData = dataForDisplayAllGroups{iGrp,1}.connAllSubjects(:,elecClusterSide);
    meanConnVsFreq{iGrp} = getMeanConn(connData,electrodeGroupList);
end

freqVals = dataForDisplayAllGroups{1,1}.freqVals;
for i=1:numElectrodeGroups
    data{1} = squeeze(meanConnVsFreq{1}(:,i,:));
    data{2} = squeeze(meanConnVsFreq{2}(:,i,:));
    displayAndcompareData(plotTopo{4}(i),data,freqVals,displaySettings,[0 1],1,useMedianFlagData,1);
    title(plotTopo{4}(i),groupNameList{i});
end

end

function chanlocs = getMontageDetails(refType)
capLayout = 'actiCap64';
clear cL bL chanlocs iElec electrodeList noseDir
switch refType
    case 'unipolar'
        cL = load([capLayout '.mat']);
        chanlocs = cL.chanlocs;
    case 'laplacian'
        cL = load([capLayout '.mat']);
        chanlocs = cL.chanlocs;
    case 'bipolar'
        cL = load(['bipolarChanlocs' capLayout '.mat']);
        chanlocs = cL.eloc;
end
end
function meanData = getMeanConn(connData,electrodeGroupList)
numSubjects = length(connData);
[numRefElectrodes,numElectrodeGroups] = size(electrodeGroupList);

sizeConnData = size(squeeze(connData{1}));
if length(sizeConnData)==2 %2-D array
    outData = zeros(numSubjects,numRefElectrodes,numElectrodeGroups);
else
    outData = zeros(numSubjects,numRefElectrodes,numElectrodeGroups,sizeConnData(3));
end

for i=1:numSubjects
    tmpConn = squeeze(connData{i});

    for j=1:numRefElectrodes
        for k=1:numElectrodeGroups
            if length(sizeConnData)==2
                outData(i,j,k) = mean(tmpConn(j,electrodeGroupList{j,k}),'omitnan');
            else
                outData(i,j,k,:) = mean(tmpConn(j,electrodeGroupList{j,k},:),2,'omitnan');
            end
        end
    end
end
meanData = squeeze(mean(outData,2));
end
function [electrodeGroupList,groupNameList,binnedCenters] = getElectrodeGroupsConn(refElectrodes,montageChanlocs)

%%%%% Discretize connectivty into bins depending on distance from seed %%%%
binWidth = 0.25;
binEdges = -1:binWidth:1;
nbins = length(binEdges)-1;
binnedCenters = binEdges(1:end-1)+(binWidth/2);
loc = getElecLocAngles(montageChanlocs);

numElectrodes = length(refElectrodes);
electrodeGroupList = cell(numElectrodes,nbins);
groupNameList = cell(1,nbins);

for e=1:numElectrodes
    dist = sqrt((angl_dist(loc.azi(refElectrodes(e)),loc.azi,'a')).^2+(angl_dist(loc.ele(refElectrodes(e)),loc.ele,'e')).^2);
    fitx = cos((dist/180)*pi);
    binned_fitx = discretize(fitx,binEdges);
    
    for b = 1:nbins % number of bins
        electrodeGroupList{e,nbins-b+1} = find(binned_fitx == b);
        if e==1
            groupNameList{nbins-b+1} = [num2str(binEdges(b)) '<x<' num2str(binEdges(b+1))];
        end
    end
end
end
function loc = getElecLocAngles(chanlocs)
azi = zeros(1,length(chanlocs)); ele = zeros(1,length(chanlocs));
for e = 1:length(chanlocs)
    azi(e) = chanlocs(e).sph_theta;
    ele(e) = chanlocs(e).sph_phi;
end
loc.azi = azi;
loc.ele = ele;
end
function out_theta = angl_dist(in_theta_ref,in_theta,val)
if(strcmp(val,'a')) % azimuth (addressing Cz issue)
    if(in_theta_ref > 90)
        in_theta(14) = 90;
    elseif(in_theta_ref < -90)
        in_theta(14) = -90;
    end
end
in_theta = abs(in_theta-in_theta_ref);
out_theta = zeros(1,length(in_theta));
for i=1:length(in_theta)
    if(in_theta(i) > 180)
        out_theta(i) = 360 - in_theta(i); % to get shortest angular distance
    else
        out_theta(i) = in_theta(i);
    end
end
if(strcmp(val,'a')) % azimuth (addressing Cz issue)
    out_theta(14) = 0;
end
end
function displayAndcompareData(hPlot,data,xs,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,nonMatchedFlag)

if ~exist('displaySignificanceFlag','var'); displaySignificanceFlag=0;  end
if ~exist('useMedianFlag','var');           useMedianFlag=1;            end
if ~exist('nonMatchedFlag','var');          nonMatchedFlag=1;           end

if useMedianFlag
    getLoc = @(g)(squeeze(median(g,1,'omitnan')));
else
    getLoc = @(g)(squeeze(mean(g,1,'omitnan')));
end

numGroups = length(data);

axes(hPlot);
for i=1:numGroups
    clear bootStat mData sData
    mData = getLoc(data{i});
    if useMedianFlag
        bootStat = bootstrp(1000,getLoc,data{i});
        sData = std(bootStat);
    else
        sData = std(data{i},[],1)/sqrt(size(data{i},1));
    end

    patch([xs';flipud(xs')],[mData'-sData';flipud(mData'+sData')],displaySettings.colorNames(i,:),'linestyle','none','FaceAlpha',0.4);
    hold on;
    plot(xs,mData,'color',displaySettings.colorNames(i,:),'linewidth',1);
end

set(gca,'fontsize',displaySettings.fontSizeLarge);
set(gca,'TickDir','out','TickLength',displaySettings.tickLengthMedium);

if exist('yLims','var') && ~isempty(yLims)
    ylim(yLims);
else
    yLims = ylim;
end

if displaySignificanceFlag % Do significance Testing

    allData = [];
    allIDs = [];
    for j=1:numGroups
        allData = cat(1,allData,data{j});
        allIDs = cat(1,allIDs,j+zeros(size(data{j},1),1));
    end

    for i=1:length(xs)
        if useMedianFlag
            p=kruskalwallis(allData(:,i),allIDs,'off');
        else
            if nonMatchedFlag
                [~,p]=ttest2(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
            else
                [~,p]=ttest(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
            end
        end
        % Get patch coordinates
        yVals = yLims(1)+[0 0 diff(yLims)/20 diff(yLims)/20];

        clear xMidPos xBegPos xEndPos
        xMidPos = xs(i);
        if i==1
            xBegPos = xMidPos;
        else
            xBegPos = xMidPos-(xs(i)-xs(i-1))/2;
        end
        if i==length(xs)
            xEndPos = xMidPos;
        else
            xEndPos = xMidPos+(xs(i+1)-xs(i))/2;
        end
        clear xVals; xVals = [xBegPos xEndPos xEndPos xBegPos]';

        if (p<0.05)
            patch(xVals,yVals,'c','linestyle','none');
        end
        if (p<0.01)
            patch(xVals,yVals,'k','linestyle','none');
        end
    end
end
end