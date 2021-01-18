% This program displays the data saved in the analyzedData folder

% subjectNameList is a cell array, with each cell containing a list of
% subjects. For example, subjectNameList{1} could contain all MCI/AD subjects
% while subjectNameList{2} could contain all their controls.

% Option added to send a 2D array of subjectList in which each entry is a
% list of subjects.

function displayAnalyzedData(folderSourceString,subjectNameLists,strList,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,useMedianFlag,spatialFrequenciesToRemove,useCleanData)

if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('gamma1Range','var');     gamma1Range = [20 34];              end
if ~exist('gamma2Range','var');     gamma2Range = [36 66];              end
if ~exist('alphaRange','var');      alphaRange = [8 12];                end
if ~exist('useMedianFlag','var');   useMedianFlag = 1;                  end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end
if ~exist('useCleanData','var');    useCleanData=0;                     end

numGroups = length(subjectNameLists);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeLims = [-0.5 1.2]; freqLims = [0 100]; BLPSDLims = [-3 3];

if strcmpi(protocolType,'SF_ORI')
    cLims = [-2 2];
    barLims = [-1 2];
elseif strcmpi(protocolType,'TFCP')
    cLims = [-4 12];
    barLims = [4 10];
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10; displaySettings.tickLengthMedium = [0.025 0];

if numGroups==2 % Murty's color scheme
    colormap magma;
    colorNames = hot(8); colorNames([1:3,end-2:end],:) = [];
    displaySettings.colorNames = colorNames;
else
    displaySettings.colorNames = hot(numGroups);
end

noseDir = '+X';
chanlocs = getMontageDetails(refType);

barWidth = 0.4;
barPosList = (barWidth/2)*(-(numGroups-1):2:(numGroups-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataForDisplayAllGroups = cell(1,numGroups);
nonMatchedFlag = ischar(subjectNameLists{1}{1});

for i=1:numGroups
    disp(['Getting data for group: ' strList{i}]);
    subjectNameList = subjectNameLists{i};
    if nonMatchedFlag
        dataForDisplayAllGroups{i} = combineAnalyzedData(folderSourceString,subjectNameList,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,spatialFrequenciesToRemove,useCleanData);
    else
        numLists = length(subjectNameList);
        dataForDisplayAllGroupsTMP = cell(1,numLists);
        for j=1:numLists
            dataForDisplayAllGroupsTMP{j} = combineAnalyzedData(folderSourceString,subjectNameList{j},projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,spatialFrequenciesToRemove,useCleanData);
        end
        dataForDisplayAllGroups{i} = combineDataForDisplay(dataForDisplayAllGroupsTMP,0);
    end
end

numRanges = length(dataForDisplayAllGroups{1}.rangeNames);

strListDisplay = cell(1,numGroups);
for i=1:numGroups
    strListDisplay{i} = [strList{i} '(' num2str(size(dataForDisplayAllGroups{i}.erpData,1)) ')'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hERP = subplot('Position',[0.05 0.8 0.15 0.15]);
hDeltaPSD = subplot('Position',[0.275 0.8 0.125 0.15]);
hBLPSD = subplot('Position',[0.45 0.8 0.125 0.15]);
hBar = subplot('Position',[0.625 0.8 0.125 0.15]);

hDeltaTF = getPlotHandles(numGroups,1,[0.05 0.1 0.15 0.55],0,0.05); linkaxes(hDeltaTF);
hTopo = getPlotHandles(numGroups,numRanges,[0.275 0.1 0.475 0.55],0.05,0.05);

dataERP = cell(1,numGroups);
dataDeltaPSD = cell(1,numGroups);
dataBLPSD = cell(1,numGroups);
dataBar = cell(1,numGroups);

for i=1:numGroups
    dataForDisplay = dataForDisplayAllGroups{i};
    
    dataERP{i} = dataForDisplay.erpData;
    dataDeltaPSD{i} = 10*(dataForDisplay.logSTPowerVsFreqAllSubjects - dataForDisplay.logBLPowerVsFreqAllSubjects);
    dataBLPSD{i} = dataForDisplay.logBLPowerVsFreqAllSubjects;
    dataBar{i} = dataForDisplay.powerDBAllSubjects;
    
    % Time-Frequency plots
    if useMedianFlag
        mTF = squeeze(median(dataForDisplay.dTFPowerDBAllSubjects,3))';
    else
        mTF = squeeze(mean(dataForDisplay.dTFPowerDBAllSubjects,3))';
    end
    pcolor(hDeltaTF(i),dataForDisplay.timeValsTF,dataForDisplay.freqValsTF,mTF);
    shading(hDeltaTF(i),'interp'); 
    caxis(hDeltaTF(i),cLims); 
    xlim(hDeltaTF(i),timeLims);
    ylim(hDeltaTF(i),freqLims);
    
    set(hDeltaTF(i),'box','off');
    set(hDeltaTF(i),'fontsize',displaySettings.fontSizeLarge);
    set(hDeltaTF(i),'TickDir','out','TickLength',displaySettings.tickLengthMedium);
  
    makeBox(hDeltaTF(i),timeLims,gamma1Range,'w',1,'-','H');
    makeBox(hDeltaTF(i),timeLims,gamma2Range,'w',1,'--','H');
    makeBox(hDeltaTF(i),timeLims,alphaRange,'w',1,':','H');
    
    if i==numGroups
        set(hDeltaTF(i),'xTick',[0 0.8],'xtickLabel',[0 0.8]);
        xlabel(hDeltaTF(i),'Time (s)');
    else
        set(hDeltaTF(i),'xTick',[0 0.8],'xtickLabel',[]);
    end
    
    set(hDeltaTF(i),'yTick',[0 50 100],'YTicklabel',[0 50 100]);
    ylabel(hDeltaTF(i),'Frequency (Hz)');

    % Topoplots
    for j=1:numRanges
        axes(hTopo(i,j)); %#ok<LAXES>
        if useMedianFlag
            x = squeeze(nanmedian(dataForDisplay.powerDBTopoAllSubjects(:,j,:),3));
        else
            x = squeeze(nanmean(dataForDisplay.powerDBTopoAllSubjects(:,j,:),3));
        end

        x(isnan(x)) = 999;
        topoplot_murty(x,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir',noseDir,'emarkercolors',x); 
        caxis(cLims);
        if i==1
            title(dataForDisplay.rangeNames{j},'fontsize',displaySettings.fontSizeLarge);
        end
    end
    
    % Set Colorbar
    tmpPos = get(hDeltaTF(i),'Position');
    cbPos = [tmpPos(1)+tmpPos(3)+0.015 tmpPos(2) 0.015 tmpPos(4)];
    hCB = colorbar('Position',cbPos,'Limits',cLims);
    
    set(hCB,'ydir','normal','box','off');
    set(hCB,'xtick',[],'ytick',[cLims(1) 0 cLims(2)],'yticklabel',[cLims(1) 0 cLims(2)]);
    set(hCB,'yaxislocation','right');
    %ylabel(hCB,'dB');
    set(hCB,'fontsize',displaySettings.fontSizeLarge,'TickDir','out','TickLength',displaySettings.tickLengthMedium(1));
    
    title(hDeltaTF(i),strListDisplay{i},'color',displaySettings.colorNames(i,:));
end

% ERPs
displayAndcompareData(hERP,dataERP,dataForDisplay.timeVals,displaySettings,[],0,useMedianFlag);
xlim(hERP,timeLims);
set(hERP,'xTick',[0 0.8],'XTicklabel',[0 0.8]);
xlabel(hERP,'Time (s)');
ylabel(hERP,'Voltage (\muV)');
title(hERP,'ERP');

% DeltaPSD
smoothSigma=[]; 
displayAndcompareData(hDeltaPSD,dataDeltaPSD,dataForDisplay.freqVals,displaySettings,cLims,1,useMedianFlag,smoothSigma,nonMatchedFlag);
axis(hDeltaPSD,[freqLims cLims]);
makeBox(hDeltaPSD,gamma1Range,cLims,'k',1,'-','V');
makeBox(hDeltaPSD,gamma2Range,cLims,'k',1,'--','V'); 
makeBox(hDeltaPSD,alphaRange,cLims,'k',1,':','V');

set(hDeltaPSD,'xTick',[0 50 100],'XTicklabel',[0 50 100]);
xlabel(hDeltaPSD,'Frequency (Hz)');
ylabel(hDeltaPSD,'\DeltaPower (dB)');
title(hDeltaPSD,'Change in power');

% BaselinePSD
displayAndcompareData(hBLPSD,dataBLPSD,dataForDisplay.freqVals,displaySettings,BLPSDLims,1,useMedianFlag);
xlim(hBLPSD,freqLims);
set(hBLPSD,'xTick',[0 50 100],'XTicklabel',[0 50 100]);
xlabel(hBLPSD,'Frequency (Hz)');
ylabel(hBLPSD,'log(Power(\muV^2))');
title(hBLPSD,'Baseline PSDs');
axes(hBLPSD);

for i=1:numGroups
    text(0.3,1-0.15*i,strListDisplay{i},'color',displaySettings.colorNames(i,:),'units','normalized','fontsize',displaySettings.fontSizeLarge);
end

% BarPlots
axes(hBar);
for i=1:numRanges
    dataTMP = []; groupID = [];
    for j=1:numGroups
        d = dataBar{j}(:,i);
        dataTMP = cat(1,dataTMP,d(:));
        groupID = cat(1,groupID,j+zeros(length(d),1));
        
        if useMedianFlag
            mD = median(d);
            sD = std(bootstrp(10000,@median,d));
        else
            mD = mean(d);
            sD = std(d)/sqrt(length(d));
        end
        errorbar(i+barPosList(j),mD,sD,'linestyle','none','marker','none','color',displaySettings.colorNames(j,:)); hold on;
        bar(i+barPosList(j),mD,barWidth,'facecolor',displaySettings.colorNames(j,:),'facealpha',0.85);
        
        % Plot individual data points
        scatter(i+barPosList(j) + zeros(1,length(d)),d,30,[0.5 0.5 0.5]);
    end

    if useMedianFlag
        [pD,tblD] = kruskalwallis(dataTMP',groupID','off');
        disp([dataForDisplay.rangeNames{i} '; KW test: X2(' num2str(tblD{4,3}) ')=' num2str(round(tblD{2,5},2)) ', p=' num2str(pD)]);
    else
        if nonMatchedFlag
            [~,pD,~,stats] = ttest2(dataBar{1}(:,i),dataBar{2}(:,i));
        else
            [~,pD,~,stats] = ttest(dataBar{1}(:,i),dataBar{2}(:,i));
        end
        disp([dataForDisplay.rangeNames{i} '; t-test: tstat(' num2str(stats.df) ')=' num2str(round(stats.tstat,2)) ', p=' num2str(pD)]);
    end
end

axis([0 numRanges+1 barLims]);
ylabel('\DeltaPower (dB)');
set(gca,'xTick',1:numRanges,'xTicklabel',dataForDisplay.rangeNames);
set(gca,'fontsize',displaySettings.fontSizeLarge,'TickDir','out','TickLength',displaySettings.tickLengthMedium);

end
function displayAndcompareData(hPlot,data,xs,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,smoothSigma,nonMatchedFlag)

if ~exist('displaySignificanceFlag','var'); displaySignificanceFlag=0;  end
if ~exist('useMedianFlag','var');           useMedianFlag=1;            end
if ~exist('smoothSigma','var');             smoothSigma=[];             end
if ~exist('nonMatchedFlag','var');          nonMatchedFlag=1;           end

if useMedianFlag
    getLoc = @(g)(squeeze(median(g,1)));
else
    getLoc = @(g)(squeeze(mean(g,1)));
end

numGroups = length(data);

if ~isempty(smoothSigma)
    windowLen = 5*smoothSigma;
    window = exp(-0.5*(((1:windowLen) - (windowLen+1)/2)/smoothSigma).^2);
    window = window/sum(window); %sqrt(sum(window.^2));
    
    for i=1:numGroups
        data{i} = convn(data{i},window,'same');
    end
end

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
           patch(xVals,yVals,'k','linestyle','none');
       end
       if (p<0.01)
           patch(xVals,yVals,'g','linestyle','none');
       end
   end
end
end
function chanlocs = getMontageDetails(refType)

capLayout = 'actiCap64';
clear cL bL chanlocs iElec electrodeList noseDir
switch refType
    case 'unipolar'
        cL = load([capLayout '.mat']);
        chanlocs = cL.chanlocs;
    case 'bipolar'
        cL = load(['bipolarChanlocs' capLayout '.mat']);
        chanlocs = cL.eloc;
end
end
function Y = combineDataForDisplay(X,useMedianFlag)
N = length(X);

Y = X{1};

Y.erpData=[];
Y.logBLPowerVsFreqAllSubjects=[];
Y.logSTPowerVsFreqAllSubjects=[];
Y.powerDBAllSubjects=[];
Y.dTFPowerDBAllSubjects=[];
Y.powerDBTopoAllSubjects=[];

for i=1:N
    if useMedianFlag
        Y.erpData(i,:) = median(X{i}.erpData,1);
        Y.logBLPowerVsFreqAllSubjects(i,:) = median(X{i}.logBLPowerVsFreqAllSubjects,1);
        Y.logSTPowerVsFreqAllSubjects(i,:) = median(X{i}.logSTPowerVsFreqAllSubjects,1);
        Y.powerDBAllSubjects(i,:) = median(X{i}.powerDBAllSubjects,1);
        Y.dTFPowerDBAllSubjects(:,:,i) = squeeze(median(X{i}.dTFPowerDBAllSubjects,3));
        Y.powerDBTopoAllSubjects(:,:,i) = squeeze(nanmedian(X{i}.powerDBTopoAllSubjects,3));
    else
        Y.erpData(i,:) = mean(X{i}.erpData,1);
        Y.logBLPowerVsFreqAllSubjects(i,:) = mean(X{i}.logBLPowerVsFreqAllSubjects,1);
        Y.logSTPowerVsFreqAllSubjects(i,:) = mean(X{i}.logSTPowerVsFreqAllSubjects,1);
        Y.powerDBAllSubjects(i,:) = mean(X{i}.powerDBAllSubjects,1);
        Y.dTFPowerDBAllSubjects(:,:,i) = squeeze(mean(X{i}.dTFPowerDBAllSubjects,3));
        Y.powerDBTopoAllSubjects(:,:,i) = squeeze(nanmean(X{i}.powerDBTopoAllSubjects,3));
    end
end
end
