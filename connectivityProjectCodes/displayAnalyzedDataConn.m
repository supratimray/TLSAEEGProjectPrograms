% This program displays the data saved in the analyzedData folder

% subjectNameList is a cell array, with each cell containing a list of
% subjects. For example, subjectNameList{1} could contain all MCI/AD subjects
% while subjectNameList{2} could contain all their controls.

function displayAnalyzedDataConn(folderSourceString,subjectNameLists,methodOptions,strList,projectName,refType,protocolType,stRange,freqRanges,freqRangeNames,removeMicroSaccadesFlag,useMedianFlag,spatialFrequenciesToRemove,useCleanData,numControls)
if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('freqRanges','var')
    freqRanges{1} = [8 12]; freqRangeNames{1} = 'Alpha'; % alpha
    freqRanges{2} = [20 34]; freqRangeNames{2} = 'Slow gamma'; % slow gamma
    freqRanges{3} = [36 66]; freqRangeNames{3} = 'Fast gamma'; % fast gamma
end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('useMedianFlag','var');   useMedianFlag = 1;                  end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end
if ~exist('useCleanData','var');    useCleanData=0;                     end
numControls = methodOptions.CaseVsControl.unmatched.numControls;
numGroups = length(subjectNameLists);
numFreqRanges = length(freqRanges);
pow_label = methodOptions.pow_label;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataForDisplayAllGroups = cell(numGroups,numFreqRanges);
for i=1:numFreqRanges % Analysis done separately for each frequency
    % First, get subjectList for all groups
    pow_elec{i} = [];  % for regression
    slopes_elec{i} = []; % for regression
    if strcmp(methodOptions.comparisonType,'MidVsOld')
        if strcmp(methodOptions.powcontrolType,'unmatched')
            subjectNameListTMP{1} = subjectNameLists{1};
            subjectNameListTMP{2} = subjectNameLists{2};
        elseif strcmp(methodOptions.powcontrolType,'matched')
            subjectNameListTMP = getPowerMatchedSubjectLists(folderSourceString,subjectNameLists,projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.sideToShow,pow_label,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
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
        for iCase=1:numCases
            dataForDisplayAllGroupsTMP{2,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{2}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
            if strcmp(methodOptions.powcontrolType,'unmatched') % For each case, simply average all controls - note that this is like the matched case in the ADGammaProject
                tmp = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                sortedSubjs = sortSubjsOnExpDate(tmp,dataForDisplayAllGroupsTMP{2,iCase});
                validSortedSubjs = sortedSubjs(findValidIndx(tmp.connFreqBandsAllSubjects(sortedSubjs),methodOptions.sideToShow));
                nValidSortedSubjs = length(validSortedSubjs);
                if(~isempty(numControls))
                    if(numControls > nValidSortedSubjs)
                        numControls = nValidSortedSubjs;
                    end
                    dataForDisplayAllGroupsTMP{1,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase}(validSortedSubjs(1:numControls)),projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
                else % when numControls == []
                    dataForDisplayAllGroupsTMP{1,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase},projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
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
                validSubjs = findValidIndx(tmp.connFreqBandsAllSubjects,methodOptions.sideToShow);
                absDiffPowerValid = absDiffPower(validSubjs);
                findBestControl = validSubjs(absDiffPowerValid==min(absDiffPowerValid));
                dataForDisplayAllGroupsTMP{1,iCase} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1}{iCase}(findBestControl),projectName,refType,protocolType,stRange,freqRanges(i),methodOptions.connMethod,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
            end
        end
        dataForDisplayAllGroups{1,i} = combinedData(dataForDisplayAllGroupsTMP(1,:));
        dataForDisplayAllGroups{2,i} = combinedData(dataForDisplayAllGroupsTMP(2,:));
    end
end


%%%%%%%%%%%%%%%%%%%% Setting parameters for Display script %%%%%%%%%%%%%%%%%%%%
elecGroupsCell = getElectrodeList('actiCap64',refType,0,1); % has electrode divisions
elecGroups = cell(1,length(elecGroupsCell));
for iG = 1:length(elecGroupsCell)-1
    elecGroups{iG} = cell2mat(elecGroupsCell{iG+1});
end
elecGroups{4} = cell2mat(elecGroupsCell{1});
elecClusterSide = methodOptions.sideToShow; % for ease
chanlocs = getMontageDetails(refType);
colorMrkr = {[0 0 0],[0.5 0.5 0.5]};
case_clrs = {[0.9294    0.6941    0.1255], [0.8510    0.3255    0.0941]};
if(strcmp(refType,'laplacian'))
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
    barLims = [0 0.4];
end
xtickSet = [-1:0.5:1];
ElecGLabels = {'Left', 'Right', 'Back'};
bN = 1000; % number of bootstrap iterations for median SEM
freqVals = dataForDisplayAllGroups{1,1}.freqVals;
alphaErrorbar = 0.65; % transparency level for errorbars in connectivity with inter-electrode distance plot
%%%%%%%%%%%%%%%%%%%%%%%%%% Display Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10; displaySettings.tickLengthMedium = [0.025 0];
%%%%%%%%%%%%%%%%%%%%%%%%%% Display Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(methodOptions.powcontrolType,'unmatched'))
    plotTopo{1} = getPlotHandles(3,2,[0.02 0.1 0.24 0.8],0.01,0.01,1); % power topoplots
    plotTopo{2} = getPlotHandles(3,2,[0.30 0.1 0.24 0.8],0.01,0.01,1); % conn topoplots
    plotTopo{3} = getPlotHandles(3,1,[0.6 0.1 0.16 0.8],0.01,0.01,1); % conn-distance profile
    plotTopo{4} = getPlotHandles(3,1,[0.82 0.1 0.16 0.8],0.01,0.01,1); % bar plot - one side
    for freqBand = 1:3
        % 1. change in power topoplots
        for iGrp = 1:2
            axS = subplot(plotTopo{1}(freqBand,iGrp)); % mid age
            selpowTopos = dataForDisplayAllGroups{iGrp,freqBand}.diffPowerTopoAllSubjects;
            avgpowTopos = nanmean(cell2mat(selpowTopos'),2);
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
    plotTopo{1} = getPlotHandles(3,1,[0.02 0.1 0.24 0.8],0.01,0.01,1); % del. psd - one side
    plotTopo{2} = getPlotHandles(3,2,[0.30 0.1 0.24 0.8],0.01,0.01,1); % conn topoplots
    plotTopo{3} = getPlotHandles(3,1,[0.6 0.1 0.16 0.8],0.01,0.01,1); % conn-distance profile
    plotTopo{4} = getPlotHandles(3,1,[0.82 0.1 0.16 0.8],0.01,0.01,1); % bar plot - one side
    for freqBand = 1:3
        % 1. change in PSD traces
        subplot(plotTopo{1}(freqBand)); hold on;
        for iGrp = 1:2
            BLpsd_allSubjs = cell2mat(dataForDisplayAllGroups{iGrp,freqBand}.logBLPowerVsFreqAllSubjects(:,elecClusterSide));
            STpsd_allSubjs = cell2mat(dataForDisplayAllGroups{iGrp,freqBand}.logSTPowerVsFreqAllSubjects(:,elecClusterSide));
            BLpsd = nanmean(BLpsd_allSubjs,1);
            STpsd = nanmean(STpsd_allSubjs,1);
            mDELpsd = 10*(STpsd-BLpsd);
            Nsubs(iGrp) = length(dataForDisplayAllGroups{iGrp,freqBand}.logBLPowerVsFreqAllSubjects);
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

for freqBand = 1:3
    % 2. conn topoplots
    topoData = cell(numGroups,1);
    binnedMeanConn = cell(numGroups,1);
    nSubs = zeros(1,2);
    for iGrp = 1:2
        axS = subplot(plotTopo{2}(freqBand,iGrp));
        ConnData = dataForDisplayAllGroups{iGrp,freqBand}.connFreqBandsAllSubjects;
        Nsubs(iGrp) = size(ConnData,1);
        ppcConn = ConnData(:,elecClusterSide);
        toposList = ppcConn;
        toposListAvg = cellfun(@selectFreqBandConnAvg,ppcConn,'UniformOutput',false);
        outConn = nanmean(cell2mat(toposListAvg'),2);
        topoData{iGrp} = avgAcrossSubjs(toposList); %outConn';
        topoplot(outConn,chanlocs,'emarker',{'.','k',10,1}); hold on;
        topoplot_santosh(outConn,chanlocs,'style','contour','plotrad',0.5,'headrad',0,'ccolor','r','numcontour',3);
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
        binnedMeanConn{iGrp} = getMeanConn(ppcConn,elecGroups{elecClusterSide},chanlocs,meanAnglLimits);

        % 3. conn-distance profile
        % "Conn. with distance & color coding"
        subplot(plotTopo{3}(freqBand));
        [h{iGrp},binned_Y{iGrp}] = plotIndivConnData(ppcConn,plotBinWidth,elecGroups{elecClusterSide},chanlocs,colorMrkr{iGrp});
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
    set(plotTopo{3}(freqBand),'Xdir','reverse');
    ylim(cLims);
    xlim(distXlim);
    if(freqBand == 3)
        lgdSpread = legend('show');
        lgdSpread.String = concatStrnum(strList,Nsubs);
        lgdSpread.FontSize = 8;
    end
    
    % plotting significance label above each bin (stat tests)
    binEdges = -1:plotBinWidth:1;
    binned_binX = binEdges(1:end-1)+(plotBinWidth/2);
    for binIndx = 1:size(binned_Y{1},2)
        d = [binned_Y{1}(:,binIndx); binned_Y{2}(:,binIndx)];
        group = [zeros(size(binned_Y{1},1),1); ones(size(binned_Y{2},1),1)];
        valid = ~isnan(d);
        pval_bin(binIndx) = kruskalwallis(d(valid),group(valid),'off');
        if(pval_bin(binIndx) < 0.01) % considering alpha level at "0.01"
            baseLevel = max(nanmedian(binned_Y{1}(:,binIndx)),nanmedian(binned_Y{2}(:,binIndx)));
            scatter(binned_binX(binIndx),baseLevel+0.1,30,'filled','p','k','HandleVisibility','off');
        end
    end
    % 4. Bar plots
    slope_elecAvg = cell(1,2);
    for iGrp = 1:2
        slope_elecAvg{iGrp} = nanmean(cell2mat(binnedMeanConn{iGrp}'),2); % averaging across individual electrodes
        nSubs(iGrp) = length(slope_elecAvg{iGrp});    
        if strcmp(methodOptions.comparisonType,'MidVsOld')
        slopes_elec{freqBand} = [slopes_elec{freqBand}; slope_elecAvg{iGrp}]; % for regression
        end
    end
    barWidth = diff(barLims);
    yticking = [barLims(1):(barWidth/2):barLims(2)];
    pval_y = yticking(end)-(yticking(end)/10);
    subplot(plotTopo{4}(freqBand));
    slope_extracted = slope_elecAvg;
    if(strcmp(methodOptions.comparisonType,'MidVsOld'))
        swarmchart([ones(nSubs(1),1); 2*ones(nSubs(2),1)],[slope_extracted{1}; slope_extracted{2}],20,'k','filled','MarkerFaceAlpha',0.5); hold on;
        if useMedianFlag
            bar([1 2],[nanmedian(slope_extracted{1}) nanmedian(slope_extracted{2})],'FaceAlpha',0.2); hold on;
            SEMdata = @(data,bN)(getSEMedian(rmmissing(data),bN));
            errorbar([1 2],[nanmedian(slope_extracted{1}) nanmedian(slope_extracted{2})],[SEMdata(slope_extracted{1},bN) SEMdata(slope_extracted{2},bN)] ,'o','Color','#CB4779','MarkerFaceColor','w','LineWidth',1.5);
            pvals(1,freqBand) = compare_unequalSets(slope_extracted{1},slope_extracted{2},4);
        else
            bar([1 2],[nanmean(slope_extracted{1}) nanmean(slope_extracted{2})],'FaceAlpha',0.2); hold on;
            SEMdata = @(data)std(data,[],1)/sqrt(size(data,1));
            errorbar([1 2],[nanmean(slope_extracted{1}) nanmean(slope_extracted{2})],[SEMdata(slope_extracted{1},bN) SEMdata(slope_extracted{2},bN)] ,'o','Color','#CB4779','MarkerFaceColor','w','LineWidth',1.5);
            [~,pvals(1,freqBand),~,~] = ttest2(slope_extracted{1},slope_extracted{2},'tail','right');
        end
    else % strcmp(methodOptions.comparisonType,'CaseVsControl')
        if useMedianFlag
            bar([1 2],[nanmedian(slope_extracted{1}) nanmedian(slope_extracted{2})],'m','FaceAlpha',0.2); hold on;
%             pvals(1,freqBand) = compare_unequalSets(slope_extracted{1},slope_extracted{2},4);
            pvals(1,freqBand) = signrank(slope_extracted{1},slope_extracted{2},'tail','right');
        else
            bar([1 2],[nanmean(slope_extracted{1}) nanmean(slope_extracted{2})],'m','FaceAlpha',0.2); hold on;
            [~,pvals(1,freqBand),~,~] = ttest(slope_extracted{1},slope_extracted{2},'tail','right');
        end
        for d=1:length(slope_extracted{1})
            plot([1 2], [slope_extracted{1}(d) slope_extracted{2}(d)],'k','Marker','o','MarkerFaceColor','w'); hold on;
        end
        set(plotTopo{4}(freqBand),'box','off');
    end
    
    text(2.3,pval_y,['p = ' num2str(pvals(1,freqBand),'%.3f')]);
    ylim(barLims);
    if(freqBand==3)
        set(plotTopo{4}(freqBand),'XTick',[1 2],'YTick',yticking,'XTicklabel',strList,'YTicklabel',yticking,'TickLength',[0.02, 0.025],'XTickLabelRotation',20);
        ytickformat('%.2f');
    else
        if(freqBand == 1)
            title(['Slope parameter-' ElecGLabels{elecClusterSide}],'FontSize',displaySettings.fontSizeLarge,'FontWeight','normal');
        end
        set(plotTopo{4}(freqBand),'XTick',[1 2],'YTick',yticking,'XTicklabel',[],'YTicklabel',[],'TickLength',[0.02 0.025]);
    end
    set(plotTopo{4}(freqBand),'FontSize',displaySettings.fontSizeLarge);
end

% generating connectivity vs. inter-electrode distance profile for each
% case and corresponding control subject
if(strcmp(methodOptions.powcontrolType,'unmatched')) % the following is same for 'unmatched' & 'matched'. Hence, running it for only 'unmatched'
if(strcmp(methodOptions.comparisonType,'CaseVsControl'))
for freqBand = 1:3
    figR2 = plotIndivSubjConnProfiles(dataForDisplayAllGroups(:,freqBand),elecClusterSide,methodOptions,elecGroups{elecClusterSide},chanlocs); % dataForDisplayAllGroups(:,freqBand), elecGroups{elecClusterSide}
    print(figR2,'-painters',fullfile(pwd,['Fig5_EachSubj_' freqRangeNames{freqBand} '_' methodOptions.connMethod '_' ElecGLabels{elecClusterSide}]),'-dtiff','-r300');
end
else % for 'MidVsOld' condition
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
        xlabel('Diff. power');
        ylabel('Mean connectivity');
        ylim([0 0.5]);
        set(gca,'FontSize',6);
    end
    print(hf,'-painters',fullfile(pwd,['SuppFig1_SlopeVsPow_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{elecClusterSide}]),'-dtiff','-r300');
    
    hf = figure;
    ages = methodOptions.MidVsOld.ageList;
    ages1 = ages<65;
    ages2 = ages>=65;
    ages_inorder = [ages(ages1) ages(ages2)];
    for freqBand = 1:3
        subplot(3,1,freqBand);
        scatter(ages_inorder,slopes_elec{freqBand},30,'filled');
        mdl = fitlm(ages_inorder,slopes_elec{freqBand});
        bCoeffs = table2array(mdl.Coefficients(:,1));
        text(75,0.45,['slope = ' num2str(bCoeffs(1)) '  + ' num2str(bCoeffs(2)) ' *pow'],'FontSize',6,'BackgroundColor','w');
        text(75,0.4,['pval = ' num2str(table2array(mdl.Coefficients(2,4)),'%.3e')],'FontSize',8);
        xlabel('Age');
        ylabel('Mean connectivity');
        ylim([0 0.5]);
        set(gca,'FontSize',6);
    end
    print(hf,'-painters',fullfile(pwd,['SuppFig2_SlopeVsAge_' pow_label '_' methodOptions.connMethod '_' ElecGLabels{elecClusterSide}]),'-dtiff','-r300');
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
function matchedSubjectNameLists = getPowerMatchedSubjectLists(folderSourceString,subjectNameLists,projectName,refType,protocolType,stRange,freqRange,sideToShow,pow_label,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData)

powerMatchingBinWidth = 0.5; % dB

dataForDisplay{1} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{1},projectName,refType,protocolType,stRange,freqRange,[],removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
dataForDisplay{2} = combineAnalyzedDataConn(folderSourceString,subjectNameLists{2},projectName,refType,protocolType,stRange,freqRange,[],removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);

if(strcmp(pow_label,'diff'))
    tmp = cell2mat(dataForDisplay{1}.diffPowerAllSubjects'); dataToMatch{1} = tmp(sideToShow,:);
    tmp = cell2mat(dataForDisplay{2}.diffPowerAllSubjects'); dataToMatch{2} = tmp(sideToShow,:);
elseif(strcmp(pow_label,'st'))
    tmp = cell2mat(dataForDisplay{1}.stPowerAllSubjects'); dataToMatch{1} = tmp(sideToShow,:);
    tmp = cell2mat(dataForDisplay{2}.stPowerAllSubjects'); dataToMatch{2} = tmp(sideToShow,:);
elseif(strcmp(pow_label,'bl'))
    tmp = cell2mat(dataForDisplay{1}.blPowerAllSubjects'); dataToMatch{1} = tmp(sideToShow,:);
    tmp = cell2mat(dataForDisplay{2}.blPowerAllSubjects'); dataToMatch{2} = tmp(sideToShow,:);
end

minVal = min(min(dataToMatch{1}),min(dataToMatch{2}));
maxVal = max(max(dataToMatch{1}),max(dataToMatch{2}));

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
function out = combinedData(input)
dataCase = averageAllSubjects(input,0);
out = structCatAllSubjects(dataCase);
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
function figR2 = plotIndivSubjConnProfiles(dataForDisplayAllGroups,elecClusterSide,methodOptions,refElecs,chanlocs) % dataForDisplayAllGroups(:,freqBand), elecGroups{elecClusterSide}
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
figR2.PaperSize = [16 36];
figR2.PaperOrientation = 'Portrait';
figR2.PaperPosition = [0 0 figR2.PaperSize];
figR2.Color = [1 1 1]; % White background
plotTopoSup = getPlotHandles(numValidCases,3,[0.1 0.05 0.85 0.9],0.06,0.01,1);
% plotting topoplot (1), conn. profile (2)
for caseIndx = 1:numValidCases
    for plotSide = 1:3
        subplot(plotTopoSup(caseIndx,plotSide));
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
            [~,binned_Y{iGrp}] = plotIndivConnData(dataForDisplayAllGroups{iGrp}.connFreqBandsAllSubjects(caseIndx,elecClusterSide),plotBinWidth,refElecs,chanlocs,colorMrkr{iGrp}); hold on;
            set(plotTopoSup(caseIndx,plotSide),'Xdir','reverse');
            ylim(cLims);
            xlim(distXlim);
        end
        xlimsSet = [-1:0.5:1];
        if(plotSide==1 && caseIndx == numValidCases)
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
        if(plotSide == 2)
            ylabel(plotTopoSup(caseIndx,plotSide),[int2str(ageList(caseIndx)) '/' genderList{caseIndx}]);
            hYLabel = get(plotTopoSup(caseIndx,plotSide),'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','middle','Position',[1.4 0.75 -1.5]);
        end
        end
    end
end
end

function out = structCatAllSubjects(data)
nCases = length(data);
out.freqVals = data{1}.freqVals;
for iCase = 1:nCases
    out.logBLPowerVsFreqAllSubjects(iCase,:) = data{iCase}.logBLPowerVsFreqAllSubjects;
    out.logSTPowerVsFreqAllSubjects(iCase,:) = data{iCase}.logSTPowerVsFreqAllSubjects;
    out.diffPowerAllSubjects(iCase,:) = {data{iCase}.diffPowerAllSubjects};
    out.diffPowerTopoAllSubjects(iCase,:) = {data{iCase}.diffPowerTopoAllSubjects};
    out.stPowerTopoAllSubjects(iCase,:) = {data{iCase}.stPowerTopoAllSubjects};
    out.connAllSubjects(iCase,:) = data{iCase}.connAllSubjects;
    out.connFreqBandsAllSubjects(iCase,:) = data{iCase}.connFreqBandsAllSubjects;
end
end
function out = findValidIndx(connFreqBandsAllSubjects,sideToShow)
rejectIndx = [];
elecThres = 40;
for i = 1:length(connFreqBandsAllSubjects) % across subjects
    finalConn = squeeze(nanmean(connFreqBandsAllSubjects{i,sideToShow},2)); % across freqBands
    if(length(find(isnan(finalConn))) > elecThres)
        rejectIndx = [rejectIndx i];
    end
end
out = setdiff(1:length(connFreqBandsAllSubjects),rejectIndx);
end
function sortedSubjs = sortSubjsOnExpDate(tmp,dataCase)
formatIn = 'ddmmyy';
ControlDates = datenum(tmp.expDate,formatIn);
CaseDate = datenum(dataCase.expDate,formatIn);
diffDates = ControlDates - CaseDate;
[~,sortedSubjs] = sort(diffDates);
end
function out = concatStrnum(str,nums)
out = cell(length(str));
for k =1:length(str)
    out{k} = strcat(str{k},' (',num2str(nums(k)),')');
end
end
function out = selectFreqBandConnAvg(input)
out = squeeze(nanmean(input(1,:,:),2)); % averaging across electrodes in this cluster
end

function [h,mean_binY2] = plotIndivConnData(allConn,binWidth,elecGroup,chanlocs,colorMrkr,plotSwitch)
if ~exist('plotSwitch','var');    plotSwitch = 1;       end
connectDiscrete =  '-o';
DiscreteVisibility = 'on';
binEdges = -1:binWidth:1;
nbins = length(binEdges)-1;
binned_binX = binEdges(1:end-1)+(binWidth/2);
loc = getElecLocAngles(chanlocs);
for i= 1:length(allConn)
    for e = 1:length(elecGroup)
        Gconn = squeeze(allConn{i}(1,e,:));
        dist = sqrt((angl_dist(loc.azi(elecGroup(e)),loc.azi,'a')).^2+(angl_dist(loc.ele(elecGroup(e)),loc.ele,'e')).^2);
        fitx = cos((dist/180)*pi);
        binned_fitx = discretize(fitx,binEdges);
        binned_binY = cell(1,nbins);
        for b = 1:nbins % number of bins
            binned_binY{b} = Gconn(binned_fitx == b);
        end
        mean_binY1(e,:) = cellfun(@nanmedian,binned_binY);
    end
    mean_binY2(i,:) = nanmedian(mean_binY1,1);
end
median_binned_binY = nanmedian(mean_binY2,1);
mSEM = @(data)std(bootstrp(1000,@nanmedian,data));
std_binned_binY = mSEM(mean_binY2);

if(plotSwitch)
    if(length(std_binned_binY)==1)
        h = errorbar(binned_binX,median_binned_binY,nan(1,8),connectDiscrete,'Color',colorMrkr,'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility',DiscreteVisibility);
    else
        h = errorbar(binned_binX,median_binned_binY,std_binned_binY,connectDiscrete,'Color',colorMrkr,'MarkerFaceColor','w','LineWidth',1.5,'HandleVisibility',DiscreteVisibility);
    end
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
        elecsInRange = angl_sep > meanAnglLimits(1) & angl_sep < meanAnglLimits(2);
        out{i}(e) = nanmean(Gconn(elecsInRange));
    end
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

