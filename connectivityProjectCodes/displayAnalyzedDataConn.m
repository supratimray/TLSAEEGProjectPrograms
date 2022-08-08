% This program displays the data saved in the analyzedData folder

% subjectNameList is a cell array, with each cell containing a list of
% subjects. For example, subjectNameList{1} could contain all MCI/AD subjects
% while subjectNameList{2} could contain all their controls.

function [slope_all,pvals] = displayAnalyzedDataConn(folderSourceString,figR,figpath,subjectNameLists,methodOptions,strList,projectName,refType,protocolType,stRange,freqRanges,freqRangeNames,removeMicroSaccadesFlag,useMedianFlagBarPlot,spatialFrequenciesToRemove,useCleanData,useMedianFlagData)
if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('freqRanges','var')
    freqRanges{1} = [8 12]; freqRangeNames{1} = 'Alpha';
    freqRanges{2} = [20 34]; freqRangeNames{2} = 'Slow gamma';
    freqRanges{3} = [36 66]; freqRangeNames{3} = 'Fast gamma';
end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('useMedianFlagBarPlot','var');  useMedianFlagBarPlot = 1;     end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end
if ~exist('useCleanData','var');    useCleanData=0;                     end
if ~exist('useMedianFlagData','var');  useMedianFlagData = 1;           end

numControls = methodOptions.CaseVsControl.unmatched.numControls;
numGroups = length(subjectNameLists);
numFreqRanges = length(freqRanges);
pow_label = methodOptions.pow_label;
combinedMatching = methodOptions.combinedMatching;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataForDisplayAllGroups = cell(numGroups,numFreqRanges);
subj_division = cell(1,3); % useful only in "matched" case for both conditions to count the subjects from the two electrode sides
noGoodControl = false(numFreqRanges,length(subjectNameLists{2})); % for marking subjects with slow gamma change in power < 0 dB

pow_elec = cell(1,numFreqRanges);
slopes_elec = cell(1,numFreqRanges);
for i=1:numFreqRanges % analysis done separately for each frequency
    pow_elec{i} = [];  % for regression with connectivity & age
    slopes_elec{i} = []; % for regression with power & age
    if strcmp(methodOptions.comparisonType,'MidVsOld')
        if strcmp(methodOptions.powcontrolType,'unmatched')
            subjectNameListTMP{1} = subjectNameLists{1};
            subjectNameListTMP{2} = subjectNameLists{2};
        elseif strcmp(methodOptions.powcontrolType,'matched')
            subjectNameListTMP = getPowerMatchedSubjectLists(folderSourceString,subjectNameLists,projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.sideToShow,pow_label,combinedMatching,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
            if(methodOptions.combineOppSide)
                if(combinedMatching)
                    subjectNameListTMP = getPowerMatchedSubjectLists(folderSourceString,subjectNameLists,projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.sideToShow,pow_label,combinedMatching,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                    csubjectNameListTMP = subjectNameListTMP(3:4);
                else
                    csubjectNameListTMP = getPowerMatchedSubjectLists(folderSourceString,subjectNameLists,projectName,refType,protocolType,stRange,freqRanges(i),getOppSide(methodOptions.sideToShow),pow_label,combinedMatching,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                end
                subj_division{i}(1) = length(subjectNameListTMP{1}); subj_division{i}(2) = length(csubjectNameListTMP{1});
                subjectNameListTMP{1} = [subjectNameListTMP{1}  csubjectNameListTMP{1}];
                subjectNameListTMP{2} = [subjectNameListTMP{2}  csubjectNameListTMP{2}]; 
            end   
        end
        dataForDisplayAllGroups{1,i} = combineAnalyzedDataConn(folderSourceString,subjectNameListTMP{1},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
        dataForDisplayAllGroups{2,i} = combineAnalyzedDataConn(folderSourceString,subjectNameListTMP{2},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
        selElecSide = @(x)(x(methodOptions.sideToShow));
        if(strcmp(methodOptions.pow_label,'st'))
            pow_elec{i} = cellfun(selElecSide,dataForDisplayAllGroups{1,i}.stPowerAllSubjects);
            pow_elec{i} = [pow_elec{i}; cellfun(selElecSide,dataForDisplayAllGroups{2,i}.stPowerAllSubjects)];
        elseif(strcmp(methodOptions.pow_label,'bl'))
            pow_elec{i} = cellfun(selElecSide,dataForDisplayAllGroups{1,i}.blPowerAllSubjects);
            pow_elec{i} = [pow_elec{i}; cellfun(selElecSide,dataForDisplayAllGroups{2,i}.blPowerAllSubjects)];
        else
            pow_elec{i} = cellfun(selElecSide,dataForDisplayAllGroups{1,i}.diffPowerAllSubjects);
            pow_elec{i} = [pow_elec{i}; cellfun(selElecSide,dataForDisplayAllGroups{2,i}.diffPowerAllSubjects)]; % for regression
        end
    elseif strcmp(methodOptions.comparisonType,'CaseVsControl')
        numCases = length(subjectNameLists{2});
        dataForDisplayAllGroupsTMP = cell(numGroups,numCases);
        subj_division{i} = [numCases numCases];
        for iCase=1:numCases
            dataForDisplayAllGroupsTMP{2,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{2}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
            if strcmp(methodOptions.powcontrolType,'unmatched') % For each case, simply average all controls - note that this is like the matched case in the ADGammaProject
                tmp = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                sortedSubjs = sortSubjsOnExpDate(tmp,dataForDisplayAllGroupsTMP{2,iCase}); % sorting Controls based on Case Exp. Date
                if(methodOptions.rejectLowPowSubjs)
                    if i==1 % For alpha, flip the power vals so that we reject the ones with high power (no alpha suppression), but only if the power values are diff
                        if strcmp(methodOptions.pow_label,'diff')
                            tmp2 = flipValues(tmp.diffPowerAllSubjects(sortedSubjs));
                        else
                            tmp2 = tmp.diffPowerAllSubjects(sortedSubjs);
                        end
                        validSortedSubjs = sortedSubjs(findValidIndx(tmp.connFreqBandsAllSubjects(sortedSubjs,:),methodOptions.sideToShow,tmp2));
                    else
                        validSortedSubjs = sortedSubjs(findValidIndx(tmp.connFreqBandsAllSubjects(sortedSubjs,:),methodOptions.sideToShow,tmp.diffPowerAllSubjects(sortedSubjs)));
                    end
                else
                    validSortedSubjs = sortedSubjs(findValidIndx(tmp.connFreqBandsAllSubjects(sortedSubjs,:),methodOptions.sideToShow));
                end
                nValidSortedSubjs = length(validSortedSubjs);
                if(~isempty(numControls))
                    if(numControls > nValidSortedSubjs)
                        numControls = nValidSortedSubjs;
                    end
                    dataForDisplayAllGroupsTMP{1,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase}(validSortedSubjs(1:numControls)),projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                else % when numControls == []
                    if (methodOptions.rejectLowPowSubjs)
                        dataForDisplayAllGroupsTMP{1,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase}(validSortedSubjs),projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                    else
                        dataForDisplayAllGroupsTMP{1,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                    end
                end
            elseif strcmp(methodOptions.powcontrolType,'matched')
                tmp = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                if(strcmp(methodOptions.pow_label,'diff')) % choosing the variable used to select closest control
                    powerCaseAllSides = cell2mat(dataForDisplayAllGroupsTMP{2,iCase}.diffPowerAllSubjects);
                    powerAllControlsAllSides = cell2mat(tmp.diffPowerAllSubjects');
                elseif(strcmp(methodOptions.pow_label,'st'))
                    powerCaseAllSides = cell2mat(dataForDisplayAllGroupsTMP{2,iCase}.stPowerAllSubjects);
                    powerAllControlsAllSides = cell2mat(tmp.stPowerAllSubjects');
                elseif(strcmp(methodOptions.pow_label,'bl'))
                    powerCaseAllSides = cell2mat(dataForDisplayAllGroupsTMP{2,iCase}.blPowerAllSubjects);
                    powerAllControlsAllSides = cell2mat(tmp.blPowerAllSubjects');
                end
                absDiffPower = abs(powerAllControlsAllSides(methodOptions.sideToShow,:) - powerCaseAllSides(methodOptions.sideToShow));
                if(methodOptions.rejectLowPowSubjs) % ignoring controls with less power for SG/FG and more power for alpha
                    if i==1 && strcmp(methodOptions.pow_label,'diff')
                        controlNegative = powerAllControlsAllSides(methodOptions.sideToShow,:)>0;
                    else
                        controlNegative = powerAllControlsAllSides(methodOptions.sideToShow,:)<0;
                    end
                    absDiffPower(controlNegative) = NaN;
                end
                validSubjs = findValidIndx(tmp.connFreqBandsAllSubjects,methodOptions.sideToShow);
                absDiffPowerValid = absDiffPower(validSubjs); % the  control should be on higher end over case ALWAYS!!!
                findBestControl = validSubjs(absDiffPowerValid==min(absDiffPowerValid));
                if(methodOptions.combineOppSide)
                    absDiffPower = abs(powerAllControlsAllSides(getOppSide(methodOptions.sideToShow),:) - powerCaseAllSides(getOppSide(methodOptions.sideToShow)));
                    if(methodOptions.rejectLowPowSubjs) % ignoring controls with less (more) power for SG/FG (alpha)
                        if (i==1) && strcmp(methodOptions.pow_label,'diff')
                            controlNegative = powerAllControlsAllSides(getOppSide(methodOptions.sideToShow),:)>0;
                        else
                            controlNegative = powerAllControlsAllSides(getOppSide(methodOptions.sideToShow),:)<0;
                        end
                        absDiffPower(controlNegative) = NaN;
                    end
                    validSubjs = findValidIndx(tmp.connFreqBandsAllSubjects,getOppSide(methodOptions.sideToShow));
                    absDiffPowerValid = absDiffPower(validSubjs);
                    findBestControl = cat(2,findBestControl,validSubjs(absDiffPowerValid==min(absDiffPowerValid)));
                end
                dataForDisplayAllGroupsTMP{1,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase}(findBestControl),projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
            end
            if(isempty(findBestControl))
                noGoodControl(i,iCase) = true;
            end
        end
        dataForDisplayAllGroupsTMP = dataForDisplayAllGroupsTMP(:,~noGoodControl(i,:));
        dataForDisplayAllGroups{1,i} = combinedData(dataForDisplayAllGroupsTMP(1,:),methodOptions.powcontrolType,methodOptions.combineOppSide);
        dataForDisplayAllGroups{2,i} = combinedData(dataForDisplayAllGroupsTMP(2,:),methodOptions.powcontrolType,methodOptions.combineOppSide);
        selElecSide = @(x)(x(methodOptions.sideToShow));
        if(strcmp(methodOptions.pow_label,'st'))
            pow_elec{i} = cellfun(selElecSide,dataForDisplayAllGroups{2,i}.stPowerAllSubjects);
        elseif(strcmp(methodOptions.pow_label,'bl'))
            pow_elec{i} = cellfun(selElecSide,dataForDisplayAllGroups{2,i}.blPowerAllSubjects);
        else
            pow_elec{i} = cellfun(selElecSide,dataForDisplayAllGroups{2,i}.diffPowerAllSubjects);
        end
    end
    
    selectLeftRight = @(x)(x{1}(1:2)); % selecting left and right electrode groups
    if(strcmp(methodOptions.comparisonType,'CaseVsControl') && methodOptions.rejectLowPowSubjs) % discarding cases if needed
        numCases = size(dataForDisplayAllGroups{2,i}.connAllSubjects,1);
        pow_bothSides = cellfun(selectLeftRight,dataForDisplayAllGroups{2,i}.diffPowerAllSubjects,'UniformOutput',false);
        avgPow = mean(cell2mat(pow_bothSides'),1); 
        if i==1 && strcmp(methodOptions.pow_label,'diff')
            goodCases = avgPow(~noGoodControl(i,:)) < 0; % criteria for detecting cases with alpha pow < 0 dB
        else
            goodCases = avgPow(~noGoodControl(i,:)) > 0; % criteria for detecting cases with SG/FG pow > 0 dB
        end
        % Santosh to modify
        out(1,:) = selectCases(dataForDisplayAllGroups(1,:),[goodCases goodCases]); % repition is to account for both the electrode sides
        out(2,:) = selectCases(dataForDisplayAllGroups(2,:),goodCases);
        dataForDisplayAllGroups = out;
        subj_division{i} = [nnz(goodCases) nnz(goodCases)];
        fprintf('Out of %d cases, %d cases are selected\n',numCases,nnz(goodCases));
    end
end

%%%%%%%%%%%%%%%%%%%% Setting parameters for Display script %%%%%%%%%%%%%%%%%%%%
elecGroupsCell = getElectrodeList('actiCap64',refType,0,1); % has electrode divisions
elecGroups = cell(1,length(elecGroupsCell));
for iG = 1:length(elecGroupsCell)-1
    elecGroups{iG} = cell2mat(elecGroupsCell{iG+1});
end
elecGroups{length(elecGroupsCell)} = cell2mat(elecGroupsCell{1});
elecClusterSide = methodOptions.sideToShow; % for ease
chanlocs = getMontageDetails(refType);
colorMrkr = {[0 0 0],[0.5 0.5 0.5]};
case_clrs = {[0.9294    0.6941    0.1255], [0.8510    0.3255    0.0941]};
if(strcmp(refType,'laplacian')) % this conditional is not being used now
    cLims = [0 0.5]; % connectivity value limits
    distXlim = [0.5 1]; % xlim for connectivity vs. inter-electrode distance profile
    barLims = [0 0.2]; % ylim for connectivity spread bar plots
    plotBinWidth = 0.1; % bin width for connectivity vs. inter-electrode distance profile
    meanAnglLimits = [0.75 1]; % angle limits in cosine space for spread measure
else
    distXlim = [-1 1]; % xlim for connectivity vs. inter-electrode distance profile
    cLims = [0 1]; % connectivity value limits
    barLims = [0.5 1.2]; % ylim for connectivity spread bar plots
    plotBinWidth = 0.25; % bin width for connectivity vs. inter-electrode distance profile
    meanAnglLimits = [-0.5 0.5]; % angle limits in cosine space for spread measure
end

% modifying ylim for connectivity spread bar plots based on conn. measure
if(strcmp(methodOptions.connMethod,'coh'))
    barLims = [0.2 1]; 
elseif(strcmp(methodOptions.connMethod,'plv'))
    barLims = [0.15 0.5];
elseif(strcmp(methodOptions.connMethod,'ppc'))
    barLims = [0 0.6];
end
xtickSet = [-1:0.5:1];
ElecGLabels = {'Left', 'Right', 'Back'};
bN = 1000; % number of bootstrap iterations for median SEM
freqVals = dataForDisplayAllGroups{1,1}.freqVals;
alphaErrorbar = 0.65; % transparency level for errorbars in connectivity with inter-electrode distance plot
connProfile_pvalThres = 0.01;

% Santosh to add this function
if strcmp(methodOptions.comparisonType,'MidVsOld') && strcmp(methodOptions.powcontrolType,'unmatched')
    Fig1_sampleSubject(dataForDisplayAllGroups{2,2},113,chanlocs,elecGroups,freqRanges{2},figpath);
end

set(groot,'CurrentFigure',figR); % to change the current figure handle back to the main "figR"

%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10; displaySettings.tickLengthMedium = [0.025 0];
%%%%%%%%%%%%%%%%%%%%%%%%%% Display Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(methodOptions.powcontrolType,'unmatched'))
    plotTopo{1} = getPlotHandles(3,2,[0.05 0.1 0.20 0.8],0.01,0.01,1); % power topoplots
    plotTopo{2} = getPlotHandles(3,2,[0.30 0.1 0.22 0.8],0.01,0.01,1); % conn topoplots
    plotTopo{3} = getPlotHandles(3,1,[0.58 0.1 0.10 0.8],0.01,0.01,1); % conn-distance profile
    plotTopo{4} = getPlotHandles(3,3,[0.72 0.1 0.26 0.8],0.01,0.01,1); % bar plot - one side
    for freqBand = 1:3
        % 1. change in power topoplots
        for iGrp = 1:2
            axS = subplot(plotTopo{1}(freqBand,iGrp));
            selpowTopos = dataForDisplayAllGroups{iGrp,freqBand}.diffPowerTopoAllSubjects;
            if(useMedianFlagData)
                avgpowTopos = nanmedian(cell2mat(selpowTopos'),2); % Plots in Murty's study were with this setting
            else
                avgpowTopos = nanmean(cell2mat(selpowTopos'),2);
            end
            topoplot_murty(avgpowTopos,chanlocs,'electrodes','off','style','blank','drawaxis','off','emarker',{'.','k',12,1},'emarkercolors',avgpowTopos);
            caxis([-1 1]);
            topoplot([],chanlocs,'hcolor','none','plotchans',elecGroups{elecClusterSide},'emarker',{'o','k',10,1});
            if(iGrp == 1)
                ht = text(-0.65,0,freqRangeNames{freqBand},'HorizontalAlignment','center');
                set(ht,'Rotation',90);
            end
            if(freqBand == 1)
                title(strList{1},'FontSize',displaySettings.fontSizeLarge,'FontWeight','normal');
                if(iGrp == 2)
                    cbr = colorbar(axS);
                    cbr.Position = cbr.Position + [0.07 -0.02 0.0015 0.07];
                    cbr.FontSize = displaySettings.fontSizeLarge;
                    title(strList{2},'FontSize',displaySettings.fontSizeLarge,'FontWeight','normal');
                end
            end
        end
    end
else % 'matched' case
    plotTopo{1} = getPlotHandles(3,1,[0.05 0.1 0.20 0.8],0.01,0.01,1); % del. psd - one side
    plotTopo{2} = getPlotHandles(3,2,[0.30 0.1 0.22 0.8],0.01,0.01,1); % conn topoplots
    plotTopo{3} = getPlotHandles(3,1,[0.58 0.1 0.10 0.8],0.01,0.01,1); % conn-distance profile
    plotTopo{4} = getPlotHandles(3,3,[0.72 0.1 0.26 0.8],0.01,0.01,1); % bar plot - one side
    for freqBand = 1:3
        % 1. change in PSD traces
        subplot(plotTopo{1}(freqBand)); hold on;
        line([0 100],[0 0],'Color','k');
        for iGrp = 1:2
            BLpsd_allSubjs = cell2mat(dataForDisplayAllGroups{iGrp,freqBand}.logBLPowerVsFreqAllSubjects(:,elecClusterSide));
            STpsd_allSubjs = cell2mat(dataForDisplayAllGroups{iGrp,freqBand}.logSTPowerVsFreqAllSubjects(:,elecClusterSide));
            if(methodOptions.combineOppSide)
                cBLpsd_allSubjs = cell2mat(dataForDisplayAllGroups{iGrp,freqBand}.logBLPowerVsFreqAllSubjects(:,getOppSide(elecClusterSide)));
                cSTpsd_allSubjs = cell2mat(dataForDisplayAllGroups{iGrp,freqBand}.logSTPowerVsFreqAllSubjects(:,getOppSide(elecClusterSide)));
                    BLpsd = nanmean([BLpsd_allSubjs; cBLpsd_allSubjs],1);
                    STpsd = nanmean([STpsd_allSubjs; cSTpsd_allSubjs],1);
            else
                    BLpsd = nanmean(BLpsd_allSubjs,1);
                    STpsd = nanmean(STpsd_allSubjs,1);
            end
            mDELpsd = 10*(STpsd-BLpsd);
            Nsubs(iGrp) = size(dataForDisplayAllGroups{iGrp,freqBand}.logBLPowerVsFreqAllSubjects,1);
            if(iGrp == 1 && methodOptions.combineOppSide) % because of appending both sides in the beginning
                Nsubs(iGrp) = Nsubs(iGrp)/2;
            end
            msemDELpsd = std(bootstrp(bN,@median,10*(STpsd_allSubjs-BLpsd_allSubjs)));
            ylim([-2 2]);
            patch([freqVals';flipud(freqVals')],[mDELpsd'-msemDELpsd';flipud(mDELpsd'+msemDELpsd')],case_clrs{iGrp},'linestyle','none','FaceAlpha',0.5);
            if(freqBand==3)
                set(plotTopo{1}(freqBand),'XTick',[0 50 100],'YTick',[-2 -1 0 1 2],'XTicklabel',[0 50 100],'YTicklabel',[-2 -1 0 1 2],'FontSize',displaySettings.fontSizeLarge);
                ylabel(plotTopo{1}(freqBand),'Change in power (dB)');
                xlabel(plotTopo{1}(freqBand),'Frequency (Hz)');
            else
                set(plotTopo{1}(freqBand),'XTick',[0 50 100],'YTick',[-2 -1 0 1 2],'XTicklabel',[],'YTicklabel',[],'FontSize',displaySettings.fontSizeLarge);
            end
        end
        line([freqRanges{freqBand}(1) freqRanges{freqBand}(1)],[-20 20],'LineStyle','--','Color','k');
        line([freqRanges{freqBand}(2) freqRanges{freqBand}(2)],[-20 20],'LineStyle','--','Color','k');
        if(freqBand == 3)
            ylabel({[freqRangeNames{freqBand}],'Change in power'});
        else
            ylabel({[freqRangeNames{freqBand}],'',''});
        end
        legend(plotTopo{1}(freqBand),[strList{1} ' (' int2str(Nsubs(1)) ')'],[strList{2} ' (' int2str(Nsubs(2)) ')'],'Location','southeast');
    end
end

% Replace 3 by numFreqRanges
slope_all = cell(3,3);
pvals = zeros(3,3); % freqBands x elecSides
Nsubs = zeros(2,3); % Grps x freqBands
for freqBand = 1:3
    % 2. conn topoplots
    binnedMeanConn = cell(numGroups,1);
    tbinnedMeanConn = cell(numGroups,1);
    nSubs = zeros(2,3);
    for iGrp = 1:2
        axS = subplot(plotTopo{2}(freqBand,iGrp));
        ConnData = dataForDisplayAllGroups{iGrp,freqBand}.connFreqBandsAllSubjects;
        if(methodOptions.combineOppSide)
            elecClusterSides = [elecClusterSide getOppSide(elecClusterSide)];
            if(strcmp(methodOptions.powcontrolType,'unmatched'))
                ppcConn{1} = ConnData(:,elecClusterSides(1));
                ppcConn{2} = mirrorTopoplotData(ConnData(:,elecClusterSides(2)));
            else % 'matched'
                if(strcmp(methodOptions.comparisonType,'MidVsOld'))
                    ppcConn{1} = ConnData(1:subj_division{freqBand}(1),elecClusterSides(1));
                    ppcConn{2} = mirrorTopoplotData(ConnData((subj_division{freqBand}(1)+1):(subj_division{freqBand}(1)+subj_division{freqBand}(2)),elecClusterSides(2)));
                else
                    if(iGrp==1) % control
                        ppcConn{1} = ConnData(1:subj_division{freqBand}(1),elecClusterSides(1));
                        ppcConn{2} = mirrorTopoplotData(ConnData((subj_division{freqBand}(1)+1):(subj_division{freqBand}(1)+subj_division{freqBand}(2)),elecClusterSides(2)));
                    else % cases (Grp==2)
                        ppcConn{1} = ConnData(:,elecClusterSides(1));
                        ppcConn{2} = mirrorTopoplotData(ConnData(:,elecClusterSides(2)));
                    end
                end
            end
        if(useMedianFlagData)
            outConn = nanmedian([cell2mat(ppcConn{1}') cell2mat(ppcConn{2}')],2);
        else
            outConn = nanmean([cell2mat(ppcConn{1}') cell2mat(ppcConn{2}')],2);
        end
        tbinnedMeanConn{iGrp}{1} = getMeanConn(ppcConn{1},elecGroups{elecClusterSide},chanlocs,meanAnglLimits);
        tbinnedMeanConn{iGrp}{2} = getMeanConn(ppcConn{2},elecGroups{elecClusterSide},chanlocs,meanAnglLimits);
        if(strcmp(methodOptions.powcontrolType,'unmatched') || (strcmp(methodOptions.powcontrolType,'matched') && combinedMatching))
            for k = 1:length(tbinnedMeanConn{iGrp}{1})
                binnedMeanConn{iGrp}{k} = nanmean([tbinnedMeanConn{iGrp}{1}{k}; tbinnedMeanConn{iGrp}{2}{k}],1);
            end
        else % matched && (~combinedMatching)
            binnedMeanConn{iGrp} = [tbinnedMeanConn{iGrp}{1} tbinnedMeanConn{iGrp}{2}]; % because 'matched' subjects can't be averaged (different number of subjects for each side)
        end
        tbinnedMeanConn{iGrp}{3} = binnedMeanConn{iGrp}; % combined "meanConn" (across left & right)
        else % ~(methodOptions.combineOppSide)
            ppcConn = ConnData(:,elecClusterSide);
            if(useMedianFlagData)
                outConn = nanmedian(cell2mat(ppcConn'),2);
            else
                outConn = nanmean(cell2mat(ppcConn'),2);
            end
            binnedMeanConn{iGrp} = getMeanConn(ppcConn,elecGroups{elecClusterSide},chanlocs,meanAnglLimits);
        end
        Nsubs(iGrp,freqBand) = length(binnedMeanConn{iGrp});
        
        topoplot(outConn,chanlocs,'emarker',{'.','k',10,1}); hold on;
        topoplot_santosh(outConn,chanlocs,'style','contour','plotrad',0.5,'headrad',0,'ccolor','r','numcontour',3); % for drawing 25% contour
        topoplot([],chanlocs,'hcolor','none','plotrad',0.5,'headrad',0,'plotchans',elecGroups{elecClusterSide},'emarker',{'o','k',8,1}); % highlighting reference electrodes
        caxis(cLims);
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

        % 3. conn-distance profile
        % "Conn. with distance & color coding"        
        if(methodOptions.combineOppSide)
            [h{iGrp},tbinned_Y{iGrp}{1}] = plotIndivConnData(ppcConn{1},plotBinWidth,elecGroups{elecClusterSide},chanlocs,colorMrkr{iGrp},useMedianFlagData,0);
            [h{iGrp},tbinned_Y{iGrp}{2}] = plotIndivConnData(ppcConn{2},plotBinWidth,elecGroups{elecClusterSide},chanlocs,colorMrkr{iGrp},useMedianFlagData,0);
            if(strcmp(methodOptions.powcontrolType,'unmatched') || (strcmp(methodOptions.powcontrolType,'matched') && combinedMatching))
                for k = 1:length(tbinnedMeanConn{iGrp}{1})
                    binned_Y{iGrp}(k,:) = nanmean([tbinned_Y{iGrp}{1}(k,:); tbinned_Y{iGrp}{2}(k,:)],1);
                end
            else % matched && (~combinedMatching)
                binned_Y{iGrp} = [tbinned_Y{iGrp}{1}; tbinned_Y{iGrp}{2}]; 
            end
        else
            [h{iGrp},binned_Y{iGrp}] = plotIndivConnData(ppcConn,plotBinWidth,elecGroups{elecClusterSide},chanlocs,colorMrkr{iGrp},useMedianFlagData,0);
        end
        mSEM = @(data)std(bootstrp(1000,@nanmedian,data)); std_binned_binY = mSEM(binned_Y{iGrp});
        if(useMedianFlagData)
            avg_binned_binY = nanmedian(binned_Y{iGrp},1); 
        else
            avg_binned_binY = nanmean(binned_Y{iGrp},1);
        end
        binEdges = -1:plotBinWidth:1; binned_binX = binEdges(1:end-1)+(plotBinWidth/2);
        subplot(plotTopo{3}(freqBand));
        h{iGrp} = errorbar(binned_binX,avg_binned_binY,std_binned_binY,'-o','Color',colorMrkr{iGrp},'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility','on');
        set([h{iGrp}.Bar, h{iGrp}.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h{iGrp}.Line.ColorData(1:3); 255*alphaErrorbar]);
        set(h{iGrp}.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h{iGrp}.Cap.EdgeColorData(1:3); 255*alphaErrorbar]);
        hold on;
    end
    if(freqBand==3)
        set(plotTopo{3}(freqBand),'XTick',xtickSet,'YTick',[cLims(1):(diff(cLims)/2):cLims(2)],'XTicklabel',xtickSet,'YTicklabel',[cLims(1):(diff(cLims)/2):cLims(2)],'TickLength',[0.02 0.025],'FontSize',displaySettings.fontSizeLarge);
        xlabel(plotTopo{3}(freqBand),'cos( \Delta \theta)','FontSize',13,'Interpreter','tex');
        ytickformat('%.2f');
    else
        if(freqBand == 1)
            title([upper(methodOptions.connMethod) ' connectivity'],'FontSize',displaySettings.fontSizeLarge,'FontWeight','normal');
        end
        set(plotTopo{3}(freqBand),'XTick',xtickSet,'YTick',[cLims(1):(diff(cLims)/2):cLims(2)],'XTicklabel',[],'YTicklabel',[],'TickLength',[0.02 0.025],'FontSize',displaySettings.fontSizeLarge);
    end
    
    if(strcmp(methodOptions.powcontrolType,'unmatched'))
        if(freqBand == 3)
            lgdSpread = legend('show');
            lgdSpread.String = concatStrnum(strList,Nsubs(:,freqBand));
            lgdSpread.FontSize = 8;        
        end
    else % matched
        lgdSpread = legend('show');
        lgdSpread.String = concatStrnum(strList,Nsubs(:,freqBand));
        lgdSpread.FontSize = 8;   
    end
        
    set(plotTopo{3}(freqBand),'Xdir','reverse');
    ylim(cLims);
    xlim(distXlim);
    
    % plotting significance label above each bin (stat. test)
    binEdges = -1:plotBinWidth:1;
    binned_binX = binEdges(1:end-1)+(plotBinWidth/2);
    for binIndx = 1:size(binned_Y{1},2)
        d = [binned_Y{1}(:,binIndx); binned_Y{2}(:,binIndx)];
        group = [zeros(size(binned_Y{1},1),1); ones(size(binned_Y{2},1),1)];
        valid = ~isnan(d);
        pval_bin(binIndx) = kruskalwallis(d(valid),group(valid),'off');
        if(pval_bin(binIndx) < connProfile_pvalThres)
            baseLevel = max(nanmedian(binned_Y{1}(:,binIndx)),nanmedian(binned_Y{2}(:,binIndx)));
            scatter(binned_binX(binIndx),baseLevel+0.1,30,'filled','p','k','HandleVisibility','off');
        end
    end
    
    % 4. Bar plots
    slope_elecAvg = cell(1,3);
    for elecSide=1:3 
        for iGrp = 1:2
            slope_elecAvg{elecSide}{iGrp} = nanmean(cell2mat(tbinnedMeanConn{iGrp}{elecSide}'),2); % averaging across individual electrodes
            nSubs(iGrp,elecSide) = length(slope_elecAvg{elecSide}{iGrp});
        end
    end
    
    % storing slope values for regression
    if strcmp(methodOptions.comparisonType,'MidVsOld')
        for iGrp = 1:2
            slopes_elec{freqBand} = [slopes_elec{freqBand}; slope_elecAvg{3}{iGrp}]; % for regression
        end
    else
        slopes_elec{freqBand} =  slope_elecAvg{3}{2}; % saving slopes for only cases
    end
    
    % saving for a scatter of healthy and MCI, AD
    if(strcmp(methodOptions.powcontrolType,'unmatched')) % only then all the subjects would be considered
        if strcmp(methodOptions.comparisonType,'MidVsOld')
            save([methodOptions.comparisonType '_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow} '.mat'],'pow_elec','slopes_elec');
        else
        save([methodOptions.comparisonType '_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{methodOptions.sideToShow} '_' methodOptions.caseType '.mat'],'pow_elec','slopes_elec'); 
        end
    end
    
    barWidth = diff(barLims);
    yticking = [barLims(1):(barWidth/2):barLims(2)];
    pval_y = yticking(end)-(yticking(end)/10);
    slope_all(freqBand,:) = slope_elecAvg;
    
    %%%%%% barplots %%%%%%
    for elecSide = 1:3 % Left, Right & Combined
    subplot(plotTopo{4}(freqBand,elecSide)); % Left side
    slope_extracted = slope_elecAvg{elecSide};
    if(strcmp(methodOptions.comparisonType,'MidVsOld'))
        swarmchart([ones(nSubs(1,elecSide),1); 2*ones(nSubs(2,elecSide),1)],[slope_extracted{1}; slope_extracted{2}],20,'k','filled','MarkerFaceAlpha',0.5); hold on;
        if useMedianFlagBarPlot
            bar([1 2],[nanmedian(slope_extracted{1}) nanmedian(slope_extracted{2})],'FaceAlpha',0.2); hold on;
            SEMdata = @(data,bN)(getSEMedian(rmmissing(data),bN));
            errorbar([1 2],[nanmedian(slope_extracted{1}) nanmedian(slope_extracted{2})],[SEMdata(slope_extracted{1},bN) SEMdata(slope_extracted{2},bN)] ,'o','Color','#CB4779','MarkerFaceColor','w','LineWidth',1.5);
            pvals(freqBand,elecSide) = compare_unequalSets(slope_extracted{1},slope_extracted{2},4);
        else
            bar([1 2],[nanmean(slope_extracted{1}) nanmean(slope_extracted{2})],'FaceAlpha',0.2); hold on;
            SEMdata = @(data)std(data,[],1)/sqrt(size(data,1));
            errorbar([1 2],[nanmean(slope_extracted{1}) nanmean(slope_extracted{2})],[SEMdata(slope_extracted{1},bN) SEMdata(slope_extracted{2},bN)] ,'o','Color','#CB4779','MarkerFaceColor','w','LineWidth',1.5);
            [~,pvals(freqBand,elecSide),~,~] = ttest2(slope_extracted{1},slope_extracted{2},'tail','right');
        end
    else % strcmp(methodOptions.comparisonType,'CaseVsControl')
        if useMedianFlagBarPlot
            bar([1 2],[nanmedian(slope_extracted{1}) nanmedian(slope_extracted{2})],'m','FaceAlpha',0.2); hold on;
            pvals(freqBand,elecSide) = signrank(slope_extracted{1},slope_extracted{2},'tail','right');
        else
            bar([1 2],[nanmean(slope_extracted{1}) nanmean(slope_extracted{2})],'m','FaceAlpha',0.2); hold on;
            [~,pvals(freqBand,elecSide),~,~] = ttest(slope_extracted{1},slope_extracted{2},'tail','right');
        end
        for d=1:length(slope_extracted{1})
            plot([1 2], [slope_extracted{1}(d) slope_extracted{2}(d)],'k','Marker','o','MarkerFaceColor','w'); hold on;
        end
        set(plotTopo{4}(freqBand,elecSide),'box','off');
    end
    text(0.5,pval_y,['p = ' num2str(pvals(freqBand,elecSide),'%.3f')],'Clipping','on');
    ylim(barLims);
    
    if(freqBand==3 && elecSide==1)
        set(plotTopo{4}(freqBand,elecSide),'XTick',[1 2],'YTick',yticking,'XTicklabel',strList,'YTicklabel',yticking,'TickLength',[0.02, 0.025],'XTickLabelRotation',20);
        ytickformat('%.2f');
    else
        set(plotTopo{4}(freqBand,elecSide),'XTick',[1 2],'YTick',yticking,'XTicklabel',[],'YTicklabel',[],'TickLength',[0.02 0.025]);
    end
    if(elecSide~=1)
        set(plotTopo{4}(freqBand,elecSide),'YColor','none','Color','none');
    end
    set(plotTopo{4}(freqBand,elecSide),'FontSize',displaySettings.fontSizeLarge);
    end
end

% generating connectivity vs. inter-electrode distance profile for each
% case and corresponding control subject
if(strcmp(methodOptions.comparisonType,'CaseVsControl'))
for freqBand = 1:3
    figR2 = plotIndivSubjConnProfiles(dataForDisplayAllGroups(:,freqBand),elecClusterSide,methodOptions,elecGroups{elecClusterSide},chanlocs,useMedianFlagData); % dataForDisplayAllGroups(:,freqBand), elecGroups{elecClusterSide}
    print(figR2,'-painters',fullfile(figpath,['Fig5_EachSubj_' methodOptions.powcontrolType '_' freqRangeNames{freqBand} '_' methodOptions.connMethod '_' ElecGLabels{elecClusterSide}]),'-dtiff','-r300');
end

else  % for 'MidVsOld' condition
    % Scatter plot of "change in pow" vs. "slope"
    if(strcmp(methodOptions.powcontrolType,'unmatched')) % the following is same for 'unmatched' & 'matched'. Hence, running it for only 'unmatched'
    hf = figure; % slopes with band power
    if(strcmp(pow_label,'diff'))
        t_xpos = [-2 1 1.5];
    else % 'st' condition settings
        t_xpos = [10 -1 -5];
    end
    for freqBand = 1:3
        subplot(3,1,freqBand);
        scatter(pow_elec{freqBand},slopes_elec{freqBand},30,'filled');
        mdl = fitlm(pow_elec{freqBand},slopes_elec{freqBand});
        bCoeffs = table2array(mdl.Coefficients(:,1));
        text(t_xpos(freqBand),0.45,['slope = ' num2str(bCoeffs(1)) '  + ' num2str(bCoeffs(2)) ' *pow'],'FontSize',6,'BackgroundColor','w');
        text(t_xpos(freqBand),0.4,['pval = ' num2str(table2array(mdl.Coefficients(2,4)),'%.3e')],'FontSize',8);
        if(strcmp(methodOptions.pow_label,'diff'))
            xlabel('Diff. power in dB');
        else
            xlabel('Absolute power in log10()');
        end
        ylabel('Mean connectivity');
        ylim([0 0.5]);
        set(gca,'FontSize',6);
    end
    print(hf,'-painters',fullfile(pwd,['SuppFig1_SlopeVsPow_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{elecClusterSide}]),'-dtiff','-r300');
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
    case 'laplacian'
        cL = load([capLayout '.mat']);
        chanlocs = cL.chanlocs;
    case 'bipolar'
        cL = load(['bipolarChanlocs' capLayout '.mat']);
        chanlocs = cL.eloc;
end
end
function matchedSubjectNameLists = getPowerMatchedSubjectLists(folderSourceString,subjectNameLists,projectName,refType,protocolType,stRange,freqRange,sideToShow,pow_label,combinedMatching,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData)
powerMatchingBinWidth = 0.5; % dB
dataForDisplay{1} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1},projectName,refType,protocolType,stRange,freqRange,[],removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
dataForDisplay{2} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{2},projectName,refType,protocolType,stRange,freqRange,[],removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);

if(strcmp(pow_label,'diff'))
    tmp = cell2mat(dataForDisplay{1}.diffPowerAllSubjects'); dataToMatch{1} = tmp(sideToShow,:); cdataToMatch{1} = tmp(getOppSide(sideToShow),:); 
    tmp = cell2mat(dataForDisplay{2}.diffPowerAllSubjects'); dataToMatch{2} = tmp(sideToShow,:); cdataToMatch{2} = tmp(getOppSide(sideToShow),:);
elseif(strcmp(pow_label,'st'))
    tmp = cell2mat(dataForDisplay{1}.stPowerAllSubjects'); dataToMatch{1} = tmp(sideToShow,:); cdataToMatch{1} = tmp(getOppSide(sideToShow),:); 
    tmp = cell2mat(dataForDisplay{2}.stPowerAllSubjects'); dataToMatch{2} = tmp(sideToShow,:); cdataToMatch{2} = tmp(getOppSide(sideToShow),:);
elseif(strcmp(pow_label,'bl')) 
    tmp = cell2mat(dataForDisplay{1}.blPowerAllSubjects'); dataToMatch{1} = tmp(sideToShow,:); cdataToMatch{1} = tmp(getOppSide(sideToShow),:); 
    tmp = cell2mat(dataForDisplay{2}.blPowerAllSubjects'); dataToMatch{2} = tmp(sideToShow,:); cdataToMatch{2} = tmp(getOppSide(sideToShow),:);
end

minVal = min(min(dataToMatch{1}),min(dataToMatch{2}));
maxVal = max(max(dataToMatch{1}),max(dataToMatch{2}));
if(combinedMatching)
    cminVal = min(min(cdataToMatch{1}),min(cdataToMatch{2}));
    cmaxVal = max(max(cdataToMatch{1}),max(cdataToMatch{2}));

    minVal = min(minVal,cminVal);
    maxVal = max(maxVal,cmaxVal);
end
powerBins = minVal:powerMatchingBinWidth:maxVal;
numPowerBins = length(powerBins);
d = powerBins(2)-powerBins(1);
powerBins = [powerBins powerBins(numPowerBins)+d];

matchedSubjectNameLists{1} = [];
matchedSubjectNameLists{2} = [];
for i=1:numPowerBins
    pos1 = intersect(find(dataToMatch{1}>=powerBins(i)),find(dataToMatch{1}<powerBins(i+1)));
    pos2 = intersect(find(dataToMatch{2}>=powerBins(i)),find(dataToMatch{2}<powerBins(i+1)));
    [equalPos1,equalPos2] = getEqualNumOfIndices(pos1,pos2);
    matchedSubjectNameLists{1} = cat(2,matchedSubjectNameLists{1},subjectNameLists{1}(equalPos1));
    matchedSubjectNameLists{2} = cat(2,matchedSubjectNameLists{2},subjectNameLists{2}(equalPos2));
    if(combinedMatching)
        pos3 = intersect(find(cdataToMatch{1}>=powerBins(i)),find(cdataToMatch{1}<powerBins(i+1)));
        pos4 = intersect(find(cdataToMatch{2}>=powerBins(i)),find(cdataToMatch{2}<powerBins(i+1)));
        [cequalPos1,cequalPos2] = getEqualNumOfIndices(pos3,pos4);
        if(length(cequalPos1) <= length(equalPos1))
            % Santosh to modify this part
 %           N = length(cequalPos1); N_other = length(equalPos1);
 %           equalPos1 = equalPos1(randperm(N_other,N));
 %           equalPos2 = equalPos2(randperm(N_other,N));
        else
            N = length(equalPos1); N_other = length(cequalPos1);
            cequalPos1 = cequalPos1(randperm(N_other,N));
            cequalPos2 = cequalPos2(randperm(N_other,N));
        end
        matchedSubjectNameLists{3} = cat(2,matchedSubjectNameLists{1},subjectNameLists{1}(cequalPos1));
        matchedSubjectNameLists{4} = cat(2,matchedSubjectNameLists{2},subjectNameLists{2}(cequalPos2));
    end
end
end

function [x2,y2] = getEqualNumOfIndices(x1,y1)

N1 = length(x1);
N2 = length(y1);

if (N1==0) || (N2==0) % one of the two is an empty array
    x2=[]; y2=[];
    
elseif N1==N2
    x2=x1; y2=y1;
    
elseif N1<N2
    x2 = x1;
    randVals = randperm(N2);
    y2 = y1(sort(randVals(1:N1)));
    
else %N1>N2
    y2 = y1;
    randVals = randperm(N1);
    x2 = x1(sort(randVals(1:N2)));
end
end
function out = combinedData(input,powcontrolType,combineSides)
if(strcmp(powcontrolType,'matched'))
    if(~combineSides)
        dataCase = averageAllSubjects(input,0);
    else
        dataCase = input;
    end
    out = structCatAllSubjects(dataCase,combineSides);
else
    dataCase = averageAllSubjects(input,0);
    out = structCatAllSubjects(dataCase,0);
end
end
function out = averageAllSubjects(data,useMedianFlag)
nCases = length(data);
out = cell(1,nCases);
for iCase = 1:nCases
    out{iCase}.freqVals = data{1,iCase}.freqVals;
    if(~useMedianFlag)
        for freqBand = 1:3
            out{iCase}.logBLPowerVsFreqAllSubjects{1,freqBand} = nanmean(cell2mat(data{1,iCase}.logBLPowerVsFreqAllSubjects(:,freqBand)),1);
            out{iCase}.logSTPowerVsFreqAllSubjects{1,freqBand} = nanmean(cell2mat(data{1,iCase}.logSTPowerVsFreqAllSubjects(:,freqBand)),1);
            out{iCase}.connAllSubjects{1,freqBand} = avgAcrossSubjs(data{1,iCase}.connAllSubjects(:,freqBand));
            out{iCase}.connFreqBandsAllSubjects{1,freqBand} = nanmean(cell2mat(data{1,iCase}.connFreqBandsAllSubjects(:,freqBand)),1);
        end
        out{iCase}.diffPowerAllSubjects = mean(cell2mat(data{1,iCase}.diffPowerAllSubjects'),2);
        out{iCase}.diffPowerTopoAllSubjects = mean(cell2mat(data{1,iCase}.diffPowerTopoAllSubjects'),2);
        out{iCase}.stPowerTopoAllSubjects = mean(cell2mat(data{1,iCase}.stPowerTopoAllSubjects'),2);
    else
        for freqBand = 1:3
            out{iCase}.logBLPowerVsFreqAllSubjects{1,freqBand} = median(cell2mat(data{1,iCase}.logBLPowerVsFreqAllSubjects(:,freqBand)),1);
            out{iCase}.logSTPowerVsFreqAllSubjects{1,freqBand} = median(cell2mat(data{1,iCase}.logSTPowerVsFreqAllSubjects(:,freqBand)),1);
            out{iCase}.connAllSubjects{1,freqBand} = avgAcrossSubjs(data{1,iCase}.connAllSubjects(:,freqBand));
            out{iCase}.connFreqBandsAllSubjects{1,freqBand} = median(cell2mat(data{1,iCase}.connFreqBandsAllSubjects(:,freqBand)),1);
        end
        out{iCase}.diffPowerAllSubjects = median(cell2mat(data{1,iCase}.diffPowerAllSubjects'),2);
        out{iCase}.diffPowerTopoAllSubjects = median(cell2mat(data{1,iCase}.diffPowerTopoAllSubjects'),2);
        out{iCase}.stPowerTopoAllSubjects = median(cell2mat(data{1,iCase}.stPowerTopoAllSubjects'),2);
    end
end
end
function out = structCatAllSubjects(data,combineSides)
nCases = length(data);
out.freqVals = data{1}.freqVals;
if(~combineSides)
    for iCase = 1:nCases
        out.logBLPowerVsFreqAllSubjects(iCase,:) = data{iCase}.logBLPowerVsFreqAllSubjects;
        out.logSTPowerVsFreqAllSubjects(iCase,:) = data{iCase}.logSTPowerVsFreqAllSubjects;
        out.diffPowerAllSubjects(iCase,:) = {data{iCase}.diffPowerAllSubjects};
        out.diffPowerTopoAllSubjects(iCase,:) = {data{iCase}.diffPowerTopoAllSubjects};
        out.stPowerTopoAllSubjects(iCase,:) = {data{iCase}.stPowerTopoAllSubjects};
        out.connAllSubjects(iCase,:) = data{iCase}.connAllSubjects;
        out.connFreqBandsAllSubjects(iCase,:) = data{iCase}.connFreqBandsAllSubjects;
    end
else
    for iCase = 1:nCases
        out.logBLPowerVsFreqAllSubjects(iCase,:) = data{iCase}.logBLPowerVsFreqAllSubjects(1,:);
        out.logSTPowerVsFreqAllSubjects(iCase,:) = data{iCase}.logSTPowerVsFreqAllSubjects(1,:);
        out.diffPowerAllSubjects(iCase,:) = {data{iCase}.diffPowerAllSubjects(1,:)};
        out.diffPowerTopoAllSubjects(iCase,:) = data{iCase}.diffPowerTopoAllSubjects(1,1); % note the difference
        out.stPowerTopoAllSubjects(iCase,:) = {data{iCase}.stPowerTopoAllSubjects(1,:)};
        out.connAllSubjects(iCase,:) = data{iCase}.connAllSubjects(1,:);
        out.connFreqBandsAllSubjects(iCase,:) = data{iCase}.connFreqBandsAllSubjects(1,:);
    end
    if(size(data{iCase}.logBLPowerVsFreqAllSubjects,1)>1) % when it is controls, not case
        % as cases have only one subject and controls have two subjects from the two electrode sides
    for iCase = 1:nCases
        out.logBLPowerVsFreqAllSubjects(nCases+iCase,:) = data{iCase}.logBLPowerVsFreqAllSubjects(2,:);
        out.logSTPowerVsFreqAllSubjects(nCases+iCase,:) = data{iCase}.logSTPowerVsFreqAllSubjects(2,:);
        out.diffPowerAllSubjects(nCases+iCase,:) = {data{iCase}.diffPowerAllSubjects(2,:)};
        out.diffPowerTopoAllSubjects(nCases+iCase,:) = data{iCase}.diffPowerTopoAllSubjects(2,1);
        out.stPowerTopoAllSubjects(nCases+iCase,:) = {data{iCase}.stPowerTopoAllSubjects(2,:)};
        out.connAllSubjects(nCases+iCase,:) = data{iCase}.connAllSubjects(2,:);
        out.connFreqBandsAllSubjects(nCases+iCase,:) = data{iCase}.connFreqBandsAllSubjects(2,:);
    end
    end
end
end
function figR2 = plotIndivSubjConnProfiles(dataForDisplayAllGroups,elecClusterSide,methodOptions,refElecs,chanlocs,medianFlag) % dataForDisplayAllGroups(:,freqBand), elecGroups{elecClusterSide}
colorMrkr = {[0 0 0],[0.5 0.5 0.5]};
cLims = [0 1]; plotBinWidth = 0.25; distXlim = [-1 1]; fontSize = 10;
numValidCases = length(dataForDisplayAllGroups{2}.connFreqBandsAllSubjects);
ageList = methodOptions.CaseVsControl.ageList;
genderList = methodOptions.CaseVsControl.genderList;
% Plotting each subject (in rows) profiles for respective electrodeSides (in columns)
figR2 = figure('numbertitle', 'off','name','Figure Supp: CaseVsControl IndvSubj Conn profile');
figR2.PaperType = 'a4';
figR2.PaperUnits = 'centimeters';
% figR2.PaperSize = [18.3 28]; % Nature specifications
if(numValidCases>6)
    figR2.PaperSize = [16 36];
else
    figR2.PaperSize = [16 18];
end
figR2.PaperOrientation = 'Portrait';
figR2.PaperPosition = [0 0 figR2.PaperSize];
figR2.Color = [1 1 1]; % White background
plotTopoSup = getPlotHandles(numValidCases,5,[0.1 0.05 0.85 0.9],0.06,0.01,1);
columns = {'Controls','Cases'};
% plotting topoplot (1), conn. profile (2)
for caseIndx = 1:numValidCases
    for plotSide = 1:3
        subplot(plotTopoSup(caseIndx,plotSide));
        if(plotSide<3 && caseIndx==1)
            title(columns{plotSide});
        end
        %% plotSide == 1
        if(plotSide<3)
        iGrp = plotSide;
        toposListAvg = selectFreqBandConnAvg(dataForDisplayAllGroups{iGrp}.connFreqBandsAllSubjects{caseIndx,elecClusterSide});
        outConn = toposListAvg;
        topoplot(outConn,chanlocs,'emarker',{'.','k',10,1}); hold on;
        topoplot([],chanlocs,'hcolor','none','plotrad',0.5,'headrad',0,'plotchans',refElecs,'emarker',{'o','k',8,1}); % highlighting reference electrodes
        caxis(cLims);
        else
            %% plotSide == 2
            clear binned_Y;
            for iGrp = 1:2
                [~,binned_Y{iGrp}] = plotIndivConnData(dataForDisplayAllGroups{iGrp}.connFreqBandsAllSubjects(caseIndx,elecClusterSide),plotBinWidth,refElecs,chanlocs,colorMrkr{iGrp},medianFlag); hold on;
                set(plotTopoSup(caseIndx,plotSide),'Xdir','reverse');
                ylim(cLims);
                xlim(distXlim);
            end
        end
        xlimsSet = [-1:0.5:1];
        if(plotSide==3)
            if(caseIndx == numValidCases)
                set(plotTopoSup(caseIndx,plotSide),'XTick',xlimsSet,'YTick',[cLims(1):(diff(cLims)/2):cLims(2)],'XTicklabel',xlimsSet,'YTicklabel',[cLims(1):(diff(cLims)/2):cLims(2)],'TickLength',[0.02 0.025],'FontSize',fontSize);
                xlabel(plotTopoSup(caseIndx,plotSide),'cos( \Delta \theta)','FontSize',13,'Interpreter','tex');
                ytickformat('%.2f');
                lgdSpread = legend('show');
                lgdSpread.String = {'Controls', 'Cases'};
                lgdSpread.FontSize = 6;
                lgdSpread.Location = 'best';
                lgdSpread.Box = 'off';
                lgdSpread.LineWidth = 0.3;
            else
                set(plotTopoSup(caseIndx,plotSide),'XTick',xlimsSet,'YTick',[cLims(1):(diff(cLims)/2):cLims(2)],'XTicklabel',[],'YTicklabel',[],'TickLength',[0.02 0.025],'FontSize',fontSize);
            end
        end
        if(plotSide == 3)
            ylabel(plotTopoSup(caseIndx,plotSide),[int2str(ageList(caseIndx)) '/' genderList{caseIndx}]);
            hYLabel = get(plotTopoSup(caseIndx,plotSide),'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','middle','Position',[1.4 0.75 -1.5]);
        end
    end
    subplot(plotTopoSup(caseIndx,4)); % del PSDs
    for iGrp=1:2
        BLpsd_allSubjs = (dataForDisplayAllGroups{iGrp}.logBLPowerVsFreqAllSubjects{caseIndx,elecClusterSide});
        STpsd_allSubjs = (dataForDisplayAllGroups{iGrp}.logSTPowerVsFreqAllSubjects{caseIndx,elecClusterSide});
        BLpsd = nanmean(BLpsd_allSubjs,1);
        STpsd = nanmean(STpsd_allSubjs,1);
        mDELpsd = 10*(STpsd-BLpsd);
        ylim([-2 2]);
        plot(dataForDisplayAllGroups{1}.freqVals,mDELpsd,'LineWidth',2,'Color',colorMrkr{iGrp}); hold on;
        caxis([-1 1]);
    end
    subplot(plotTopoSup(caseIndx,5));
    for iGrp = 1:2
        topoD = (dataForDisplayAllGroups{iGrp}.diffPowerTopoAllSubjects{caseIndx});
        topoplot_murty(topoD,chanlocs,'electrodes','off','style','blank','drawaxis','off','emarker',{'.','k',10,1},'emarkercolors',topoD);
    end
end
end

function out = findValidIndx(connFreqBandsAllSubjects,sideToShow,varargin)
rejectIndx = false(size(connFreqBandsAllSubjects,1),1);
elecThres = 40;
rejectPowSubjs = false(size(connFreqBandsAllSubjects,1),1);
for i = 1:size(connFreqBandsAllSubjects,1) % across subjects
    finalConn = squeeze(nanmean(connFreqBandsAllSubjects{i,sideToShow},2)); % across freqBands
    if(length(find(isnan(finalConn))) > elecThres)
        rejectIndx(i) = true;
    end
    if(nargin>2)
        if(varargin{1}{i}(sideToShow) < 0)
            rejectPowSubjs(i) = true;
        end
    end
end
allRejectIndx = rejectIndx | rejectPowSubjs;
out = find(~allRejectIndx);
end
function sortedSubjs = sortSubjsOnExpDate(tmp,dataCase)
formatIn = 'ddmmyy';
ControlDates = datenum(tmp.expDate,formatIn);
CaseDate = datenum(dataCase.expDate,formatIn);
diffDates = ControlDates - CaseDate;
[~,sortedSubjs] = sort(diffDates);
end
function out = concatStrnum(str,nums)
out = cell(1,length(str));
for k =1:length(str)
    out{k} = strcat(str{k},' (',num2str(nums(k)),')');
end
end
function out = selectFreqBandConnAvg(input)
out = squeeze(nanmean(input(1,:,:),2)); % averaging across electrodes in this cluster
end

function [h,mean_binY2] = plotIndivConnData(allConn,binWidth,elecGroup,chanlocs,colorMrkr,medianFlag,plotSwitch)
if ~exist('plotSwitch','var');    plotSwitch = 1;       end
connectDiscrete =  '-o';
DiscreteVisibility = 'on';
binEdges = -1:binWidth:1;
nbins = length(binEdges)-1;
binned_binX = binEdges(1:end-1)+(binWidth/2);
loc = getElecLocAngles(chanlocs);
for i= 1:length(allConn)
    for e = 1%:length(elecGroup)
        Gconn = squeeze(allConn{i}(1,e,:));
        dist = sqrt((angl_dist(loc.azi(elecGroup(e)),loc.azi,'a')).^2+(angl_dist(loc.ele(elecGroup(e)),loc.ele,'e')).^2);
        fitx = cos((dist/180)*pi);
        binned_fitx = discretize(fitx,binEdges);
        binned_binY = cell(1,nbins);
        for b = 1:nbins % number of bins
            binned_binY{b} = Gconn(binned_fitx == b);
        end
        if(medianFlag)
            mean_binY1(e,:) = cellfun(@nanmedian,binned_binY); % across the inter-electrode combinations in each bin
        else
            mean_binY1(e,:) = cellfun(@nanmean,binned_binY);
        end
    end
    if(medianFlag)
        mean_binY2(i,:) = nanmedian(mean_binY1,1); % across the electrodes in each group
    else
        mean_binY2(i,:) = nanmean(mean_binY1,1);
    end
end

if(medianFlag)
    median_binned_binY = nanmedian(mean_binY2,1);
    mSEM = @(data)std(bootstrp(1000,@nanmedian,data));
else
    median_binned_binY = nanmean(mean_binY2,1);
    mSEM = @(data)std(bootstrp(1000,@nanmean,data));
end

std_binned_binY = mSEM(mean_binY2);

if(plotSwitch)
    if(length(std_binned_binY)==1)
        h = errorbar(binned_binX,median_binned_binY,nan(1,8),connectDiscrete,'Color',colorMrkr,'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility',DiscreteVisibility);
    else
        h = errorbar(binned_binX,median_binned_binY,std_binned_binY,connectDiscrete,'Color',colorMrkr,'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility',DiscreteVisibility);
    end
else
    h = [];
end
end

function out = avgAcrossSubjs(topoList)
dims = [length(topoList), size(topoList{1},2), size(topoList{1},3)];
accumTopos = zeros(dims);
for subj = 1:length(topoList)
    accumTopos(subj,:,:) = topoList{subj}(1,:,:);
end
out = squeeze(nanmean(accumTopos,1));
end

function out = getMeanConn(input,refElecs,chanlocs,meanAnglLimits)
out = cell(1,length(input));
loc = getElecLocAngles(chanlocs);
for i= 1:length(input)
    for e = 1:length(refElecs)
        dist = sqrt((angl_dist(loc.azi(refElecs(e)),loc.azi,'a')).^2+(angl_dist(loc.ele(refElecs(e)),loc.ele,'e')).^2);
        angl_sep = cos((dist/180)*pi);
        Gconn = squeeze(input{i}(1,e,:));
        elecsInRange = angl_sep >= meanAnglLimits(1) & angl_sep <= meanAnglLimits(2);
        out{i}(e) = nanmean(Gconn(elecsInRange));
    end
end
end
function out = selectCases(dataForDisplayAllGroups,goodCases)
out = dataForDisplayAllGroups;
    for freqB = 1:3
        out{freqB}.logBLPowerVsFreqAllSubjects = dataForDisplayAllGroups{freqB}.logBLPowerVsFreqAllSubjects(goodCases,:); % {5×3 cell}
        out{freqB}.logSTPowerVsFreqAllSubjects = dataForDisplayAllGroups{freqB}.logSTPowerVsFreqAllSubjects(goodCases,:); %{5×3 cell}
        out{freqB}.diffPowerAllSubjects = dataForDisplayAllGroups{freqB}.diffPowerAllSubjects(goodCases); % {5×1 cell}
        out{freqB}.diffPowerTopoAllSubjects = dataForDisplayAllGroups{freqB}.diffPowerTopoAllSubjects(goodCases); % {5×1 cell}
        out{freqB}.stPowerTopoAllSubjects = dataForDisplayAllGroups{freqB}.stPowerTopoAllSubjects(goodCases); % {5×1 cell}
        out{freqB}.connAllSubjects = dataForDisplayAllGroups{freqB}.connAllSubjects(goodCases,:); % {5×3 cell}
        out{freqB}.connFreqBandsAllSubjects = dataForDisplayAllGroups{freqB}.connFreqBandsAllSubjects(goodCases,:); % {5×3 cell}
    end
end
function loc = getElecLocAngles(chanlocs)
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
function out = getOppSide(side)
    if(side==1)
        out =2;
    elseif(side==2)
        out = 1;
    elseif(side==3)
        out = 3;
    end
end
function y=flipValues(x)
numEntries = length(x);
y=cell(numEntries,1);
for i=1:numEntries
    y{i} = -x{i};
end
end