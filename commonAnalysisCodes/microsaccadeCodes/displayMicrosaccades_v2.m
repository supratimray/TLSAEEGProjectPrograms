function figH = displayMicrosaccades_v2(folderSave,eyeDataDegXOrig,eyeDataDegYOrig,trialNums,protNumForTrial,subjectName,protocolType,FsEyes,eyeRangeMS,saveDataFlag,threshold,timeRange,minCutOff,minMSLength,fixationWindowWidth,optimiseFlag)

% Defaults
if ~exist('threshold','var') || isempty(threshold); threshold = 3; end
if ~exist('timeRange','var') || isempty(timeRange); timeRange = [-0.5 0.75]; end
if ~exist('minCutOff','var') || isempty(minCutOff); minCutOff = 10; end
if ~exist('minMSLength','var') || isempty(minMSLength); minMSLength = 0.015; end;
if ~exist('fixationWindowWidth','var') || isempty(fixationWindowWidth); fixationWindowWidth = 5; end
if ~exist('saveDataFlag','var') || isempty(saveDataFlag); saveDataFlag=0; end;
if ~exist('optimiseFlag','var') || isempty(optimiseFlag); optimiseFlag = 1; end % MD 06-JAN-2019: set to 1 as in dualGamma paper

if strcmpi(protocolType,'CON'); protType = 'Contrast'; else protType = protocolType; end
if isempty(folderSave); folderSave = fullfile(subjectFolder,protType,'eyeData'); end; % MD 19-March-2019
% folderSave = fullfile(subjectFolder,protType,'eyeData'); % MD 19-March-2019
makeDirectory(folderSave); 

FsEye = FsEyes{1}; % FsEye is constant across protocols.
if ~isfloat(FsEye)
    if strcmpi(FsEye,'-'); error('Eye-data not available'); else FsEye=str2num(FsEye); end;
end
timeValsEyePos = (eyeRangeMS(1):1000/FsEye:eyeRangeMS(2)-1000/FsEye)/1000;

timeValsEyeIndicesBL = [];
timeValsEyeIndicesStim = [];
eyeDataDegX = [];
eyeDataDegY = [];

MS=[]; eyeSpeedX=[]; eyeSpeedY=[]; cutoffX=[]; cutoffY = []; eyeSpeedMag=[]; trialVals=[];

% % If saveDataFlag, directly save MS data without displaying the figure.
% if saveDataFlag
%     recalculate_Callback;
%     return;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~,subjectName] = fileparts(subjectFolder);
figH = figure; colormap jet;

plotStartPos = 0.05; plotFullWidth = 0.9; plotHalfWidth = 0.43; plotQuarterWidth = 0.21;
plotFullHeight = plotFullWidth; plotHalfHeight = plotHalfWidth; plotQuarterHeight = plotQuarterWidth; %#ok<NASGU>
dynamicStart = plotStartPos; dynamicHeight = 0.18; dynamicGap=0.01; dynamicTextWidth = 1; backgroundColor = 'w'; %#ok<NASGU>
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; fontSizeTiny = 8; %#ok<NASGU>

plotEyePlots = getPlotHandles(3,1,[plotStartPos 0.05 plotHalfWidth 0.27],0.01,0.01,0);
plotEyePos = plotEyePlots(1,1);
plotEyeVel = plotEyePlots(2,1);
plotEyeSpeed = plotEyePlots(3,1);
plotEyeXY = subplot('position',[plotStartPos 0.38 plotQuarterWidth-2*dynamicGap plotHalfWidth-0.03]); axis square; 
plotEyeVelXY = subplot('position',[plotStartPos+plotQuarterWidth+3*dynamicGap 0.38 plotQuarterWidth-2*dynamicGap plotHalfWidth-0.03]); axis square; 
        
[plotEyeVelAllTrials,~,plotEyeVelAllTrialsPos] = getPlotHandles(4,1,[plotStartPos+plotHalfWidth+6*dynamicGap plotStartPos plotQuarterWidth-4*dynamicGap plotFullHeight],0.075,0.075,0,0);
plotSaccadesAllTrials = getPlotHandles(2,1,[plotEyeVelAllTrialsPos{2}(1) plotEyeVelAllTrialsPos{2}(2) plotEyeVelAllTrialsPos{2}(3) plotEyeVelAllTrialsPos{1}(2)+plotEyeVelAllTrialsPos{1}(4)-plotEyeVelAllTrialsPos{2}(2)],0.05,0.002,0,0);
plotNumMSPos = get(plotEyeVelAllTrials(4,1),'position');
MSFreqPlot = axes('position',[plotNumMSPos(1)+0.09 plotNumMSPos(2)+0.09 0.07 0.07]); 

plotOtherPlots = getPlotHandles(4,1,[plotStartPos+plotHalfWidth+plotQuarterWidth+8*dynamicGap plotStartPos plotQuarterWidth-4*dynamicGap plotFullHeight],0.075,0.075,0,0);

subplot(plotEyePos); set(gca,'xticklabel',[]); ylabel('deg');
subplot(plotEyeVel); set(gca,'xticklabel',[]); ylabel('deg/s');
subplot(plotEyeSpeed); xlabel('Time (s)'); ylabel('deg/s');
subplot(plotEyeXY); xlabel('Deg'); ylabel('Deg');

controlPanel1 = uipanel(figH,'unit','normalized','position',[plotStartPos 0.8 plotQuarterWidth 0.15]);
controlPanel2 = uipanel(figH,'unit','normalized','position',[plotStartPos+plotQuarterWidth+dynamicGap 0.8 plotQuarterWidth 0.15]);

uicontrol('Parent',controlPanel1,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart 1-1*(dynamicHeight+2*dynamicGap) plotFullWidth  dynamicHeight], ...
    'Style','text','String',[subjectName ': ' protocolType],'FontSize',fontSizeMedium); % Changed the text content to include protocolType: MD 06-JAN-2019

uicontrol('Parent',controlPanel1,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart 1-2*(dynamicHeight+2*dynamicGap) plotQuarterWidth  dynamicHeight], ...
    'Style','pushbutton','String','<','FontSize',fontSizeMedium,'Callback',{@previousTrial_Callback});

hTrialEdit = uicontrol('Parent',controlPanel1,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotQuarterWidth+2*dynamicGap 1-2*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','edit','String','1','FontSize',fontSizeMedium,'Callback',{@plotEyeData_Callback});

uicontrol('Parent',controlPanel1,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotQuarterWidth+4*dynamicGap+plotHalfWidth 1-2*(dynamicHeight+2*dynamicGap) plotQuarterWidth  dynamicHeight], ...
    'Style','pushbutton','String','>','FontSize',fontSizeMedium,'Callback',{@nextTrial_Callback});

uicontrol('Parent',controlPanel1,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart 1-3*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','pushbutton','String','Animate','FontSize',fontSizeMedium,'Callback',{@animate_Callback});

hAnimateSpeedDropdown = uicontrol('Parent',controlPanel1,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotHalfWidth+2*dynamicGap 1-3*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','popup','String','1X|2X|3X|5X','FontSize',fontSizeMedium);


%%%%%%%%%%%%%
uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart 1-1*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','text','String','Time range (s)','FontSize',fontSizeMedium);

hTimeRangeMin = uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotHalfWidth+2*dynamicGap 1-1*(dynamicHeight+2*dynamicGap) plotQuarterWidth  dynamicHeight], ...
    'Style','edit','String',num2str(timeRange(1)),'FontSize',fontSizeMedium);

hTimeRangeMax = uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotHalfWidth+plotQuarterWidth+4*dynamicGap 1-1*(dynamicHeight+2*dynamicGap) plotQuarterWidth  dynamicHeight], ...
    'Style','edit','String',num2str(timeRange(2)),'FontSize',fontSizeMedium);

uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart 1-2*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','text','String','Threshold','FontSize',fontSizeMedium);

hThreshold = uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotHalfWidth+2*dynamicGap 1-2*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','edit','String',num2str(threshold),'FontSize',fontSizeMedium);

uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart 1-3*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','text','String','Min cutoff (deg/s)','FontSize',fontSizeMedium);

hMinCutOff = uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotHalfWidth+2*dynamicGap 1-3*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','edit','String',num2str(minCutOff),'FontSize',fontSizeMedium);

uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart 1-4*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','text','String','Min MS length (s)','FontSize',fontSizeMedium);

hMinMSLength = uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotHalfWidth+2*dynamicGap 1-4*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','edit','String',num2str(minMSLength),'FontSize',fontSizeMedium);

uicontrol('Parent',controlPanel2,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicStart+plotQuarterWidth+2*dynamicGap 1-5*(dynamicHeight+2*dynamicGap) plotHalfWidth  dynamicHeight], ...
    'Style','pushbutton','String','Recalculate','FontSize',fontSizeMedium,'Callback',{@recalculate_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

recalculate_Callback;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function recalculate_Callback(~,~)
        
        if ~saveDataFlag
            clear timeRangeMin timeRangeMax timeRange threshold minCutOff minMSLength
            timeRangeMin = str2num(get(hTimeRangeMin,'string')); %#ok<*ST2NM>
            timeRangeMax = str2num(get(hTimeRangeMax,'string'));
            timeRange = [timeRangeMin timeRangeMax];
            threshold = str2num(get(hThreshold,'string'));
            minCutOff = str2num(get(hMinCutOff,'string'));
            minMSLength = str2num(get(hMinMSLength,'string'));
        end
        
        % do baseline correction   
        clear timeValsEyeIndicesBL timeValsEyeIndicesStim
        timeValsEyeIndicesBL = ((timeValsEyePos>=timeRange(1)) & (timeValsEyePos<=0));
        timeValsEyeIndicesStim = ((timeValsEyePos>=0) & (timeValsEyePos<=timeRange(2)));
        
        clear eyeDataDegX eyeDataDegY
        eyeDataDegX = eyeDataDegXOrig - repmat(mean(eyeDataDegXOrig(:,timeValsEyeIndicesBL|timeValsEyeIndicesStim),2),1,size(eyeDataDegXOrig,2));
        eyeDataDegY = eyeDataDegYOrig - repmat(mean(eyeDataDegYOrig(:,timeValsEyeIndicesBL|timeValsEyeIndicesStim),2),1,size(eyeDataDegYOrig,2));

        % Reject bad trials
        clear xTrialsBeyondFixWindow yTrialsBeyondFixWindow xTrialsNoSignals yTrialsNoSignals badEyeTrials
        xTrialsBeyondFixWindow = sum(abs(eyeDataDegX(:,timeValsEyeIndicesBL|timeValsEyeIndicesStim))>(fixationWindowWidth/2),2);
        yTrialsBeyondFixWindow = sum(abs(eyeDataDegY(:,timeValsEyeIndicesBL|timeValsEyeIndicesStim))>(fixationWindowWidth/2),2);
        xTrialsNoSignals = sum(abs(eyeDataDegX(:,timeValsEyeIndicesBL|timeValsEyeIndicesStim)),2);
        yTrialsNoSignals = sum(abs(eyeDataDegY(:,timeValsEyeIndicesBL|timeValsEyeIndicesStim)),2);

        badEyeTrials = find(xTrialsBeyondFixWindow>0 | yTrialsBeyondFixWindow>0 | xTrialsNoSignals==0 | yTrialsNoSignals==0);
        if badEyeTrials==0; badEyeTrials=[]; end;
        eyeDataDegX(badEyeTrials,:) = [];
        eyeDataDegY(badEyeTrials,:) = [];
        
        % Get microsaccades
        if optimiseFlag
            thresholds = 3:0.1:6;
            minMSLengths = [0.010 0.015];
            rAmpVSVelMSToCheck = NaN(length(thresholds),length(minMSLengths));
            for iThr = 1:length(thresholds)
                clear thresholdToCheck
                thresholdToCheck = thresholds(iThr);
                for iDur = 1:length(minMSLengths)                    
                    clear minMSLengthToCheck
                    minMSLengthToCheck = minMSLengths(iDur);
                    disp([num2str(thresholdToCheck) '; ' num2str(minMSLengthToCheck)]);
                    
                    clear numMSInRange peakVelocityAll MSMagAllToCheck rCheck msFreqToCheck
                    [~,~,numMSInRange,~,~,~,~,~,~,peakVelocityAllToCheck,~,MSMagAllToCheck] = getMS(eyeDataDegX,eyeDataDegY,timeValsEyePos,FsEye,timeRange,thresholdToCheck,minCutOff,minMSLengthToCheck);
                    
                    msFreqToCheck = numMSInRange./diff(timeRange);
                    if (mean(msFreqToCheck)<=0.5 || mean(msFreqToCheck)>=3); continue; end;
                    rCheck = corrcoef(log10(MSMagAllToCheck)',log10(peakVelocityAllToCheck)');
                    if numel(rCheck)==1; rCheck = rCheck.*eye(2); end;
                    rAmpVSVelMSToCheck(iThr,iDur) = rCheck(2,1);
                    
                end
            end
            
            clear rAmpVSVelMSOptimal row col
            rAmpVSVelMSOptimal = unique(max(rAmpVSVelMSToCheck(:)));
            [row,col] = find(rAmpVSVelMSToCheck == rAmpVSVelMSOptimal);
            if ~isempty(row) && ~isempty(col)
                clear threshold minMSLength
                threshold = min(thresholds(row)); minMSLength = min(minMSLengths(col));
                set(hThreshold,'string',num2str(threshold));
                set(hMinMSLength,'string',num2str(minMSLength));
            else
                disp(['Microsaccade-rate too low for subject ' subjectName '. Hence, parameters could not be optimised.']);
            end
            optimiseFlag = 0;
            recalculate_Callback;
            return;
        else
            clear MS MSLength numMSInRange eyeSpeedX eyeSpeedY cutoffX cutoffY eyeSpeedMag MSAllData peakVelocityAll MSLengthAll MSMagAll trialVals nMS rAmpVSVelMS
            [MS,MSLength,numMSInRange,eyeSpeedX,eyeSpeedY,cutoffX,cutoffY,eyeSpeedMag,MSAllData,peakVelocityAll,MSLengthAll,...
                MSMagAll,trialVals,nMS] = getMS(eyeDataDegX,eyeDataDegY,timeValsEyePos,FsEye,timeRange,threshold,minCutOff,minMSLength); %#ok<ASGLU>
            rAmpVSVelMS = corrcoef(log10(MSMagAll)',log10(peakVelocityAll)');
        end
%         assignin('base','MSAllDataExtended',MSAllDataExtended);

        clear msFreq meanMSFreq
        msFreq = numMSInRange./diff(timeRange);
        meanMSFreq = mean(msFreq); stdMSFreq = std(msFreq);

        numTrials = size(eyeSpeedX,1); 
        
%         % MD 29-04-2019
%         badTrial = find(MSMagAll>5);
%         s = cumsum(numMSInRange);
%         for iT=1:length(badTrial)
%             disp(['bad MS for trial: ' num2str(find(s>=badTrial(iT),1,'first')-1) ' or ' num2str(find(s>=badTrial(iT),1,'first'))]);
%         end

        if saveDataFlag
            protNumForTrial(badEyeTrials) = [];% MD 29-03-2019
            trialNums(badEyeTrials) = [];            
            numMSInRangePerProtocol = cell(1,length(unique(protNumForTrial))); trialNumsForMSAnalysisPerProtocol = cell(1,length(unique(protNumForTrial)));
            for iP = unique(protNumForTrial)'
                numMSInRangePerProtocol{iP} = numMSInRange(protNumForTrial==iP);
                trialNumsForMSAnalysisPerProtocol{iP} = trialNums(protNumForTrial==iP);
            end            
            
            save(fullfile(folderSave,'microsaccades.mat'),'MS','nMS','MSLength','numMSInRange','trialNumsForMSAnalysisPerProtocol','numMSInRangePerProtocol','cutoffX','cutoffY','timeRange','threshold','minCutOff','minMSLength','badEyeTrials');
            save(fullfile(folderSave,'microsaccadesStats.mat'),'MSAllData','peakVelocityAll','MSLengthAll','MSMagAll','trialVals','meanMSFreq','stdMSFreq','rAmpVSVelMS');
%             return;
        end
        
        subplot(plotSaccadesAllTrials(1,1)); cla; 
        bar(timeValsEyePos,sum(trialVals,1)); 
        xlim(timeRange);
        set(gca,'xticklabels',[]);
        ylabel('No. of MS');

        subplot(plotSaccadesAllTrials(2,1)); cla;  
        imagesc(timeValsEyePos,1:numTrials,trialVals); 
        set(gca,'ydir','normal'); caxis([-1 1]);
        ylabel('Trials'); xlabel('Time (s)');
        xlim(timeRange);
        % assignin('base','trialVals',trialVals);
        box off;

        subplot(plotEyeVelAllTrials(3,1)); cla; 

        for iTrial = 1:numTrials
            clear MSForTrial MSToDiscard
            MSForTrial = MS{iTrial};
            MSToDiscard = MSForTrial<timeRange(1)|MSForTrial>timeRange(2);
            MSForTrial(MSToDiscard)=[];
            for iMS = 1:length(MSForTrial)
                clear MSData
                MSData = MSAllData{iTrial}{iMS};
                plot(MSData(1,:),MSData(2,:)); hold on;
            end
        end
        ylabel('deg'); xlabel('deg');
        axis square; axis tight;

        subplot(plotEyeVelAllTrials(4,1)); cla;  
        hist(numMSInRange,min(numMSInRange):max(numMSInRange));
        xlim([min(numMSInRange)-0.5 max(numMSInRange)+0.5]);
        text(0.8,0.9,['tr = ' num2str(length(numMSInRange))],'unit','normalized','parent',plotEyeVelAllTrials(4,1));
        text(0.8,0.8,['n = ' num2str(sum(numMSInRange))],'unit','normalized','parent',plotEyeVelAllTrials(4,1));
        ylabel('# Trials'); xlabel('No. of MS');
        
        axes(MSFreqPlot); cla; hold on; 
        hist(msFreq,min(msFreq):1/diff(timeRange):max(msFreq));
        yLimsMSFreqPlot = ylim; 
        plot(meanMSFreq+zeros(1,100),linspace(yLimsMSFreqPlot(1),yLimsMSFreqPlot(2),100),'r')
        makeBox(gca,[meanMSFreq-stdMSFreq meanMSFreq+stdMSFreq],yLimsMSFreqPlot,'r',[],':','V');
        text(0.1,0.9,[num2str(round(meanMSFreq,2),2) ', ' num2str(round(stdMSFreq,2),2)],'unit','normalized','parent',gca);
        xlim([min(msFreq)-0.2 max(msFreq)+0.2]);
        ylabel('# Trials'); xlabel('Hz'); set(gca,'color','none');
        

        subplot(plotOtherPlots(1,1)); cla;  
        hist(MSLengthAll,min(MSLengthAll):1/FsEye:max(MSLengthAll));
        xlim([min(MSLengthAll)-1/FsEye max(MSLengthAll)+1/FsEye]);
        ylabel('No. of MS'); xlabel('Time (s)');

        subplot(plotOtherPlots(2,1)); cla;  
        hist(MSMagAll,min(MSMagAll):0.01:max(MSMagAll));
        xlim([min(MSMagAll)-0.1 max(MSMagAll)+0.1]);
        ylabel('No. of MS'); xlabel('Max displacement (deg)');

        subplot(plotOtherPlots(3,1)); cla;  
        hist(peakVelocityAll,min(peakVelocityAll):max(peakVelocityAll));
        xlim([min(peakVelocityAll)-0.5 max(peakVelocityAll)+0.5]);
        ylabel('No. of MS'); xlabel('Peak velocity (deg/s)');

        subplot(plotOtherPlots(4,1)); cla;  
        plot(log10(MSMagAll),log10(peakVelocityAll),'.'); axis square;
        xlim([-2 2]); ylim([0 3]);     
        text(0.1,0.9,['m=' num2str(round(meanMSFreq,2),2) ', sd=' num2str(round(stdMSFreq,2),2) ' /s'],'unit','normalized','parent',plotOtherPlots(4,1));
        text(0.1,0.8,['n = ' num2str(nMS)],'unit','normalized','parent',plotOtherPlots(4,1));
        text(0.1,0.7,['r = ' num2str(round(rAmpVSVelMS(2,1),2),2)],'unit','normalized','parent',plotOtherPlots(4,1));
%         text(0.1,0.2,['th = ' num2str(round(threshold,2),2)],'unit','normalized','parent',plotOtherPlots(4,1));
%         text(0.1,0.1,['min MS Dur = ' num2str(round(minMSLength,3),3)],'unit','normalized','parent',plotOtherPlots(4,1));
        ylabel('log10(Peak velocity)'); xlabel('log10(Max displacement)');

        plotEyeData_Callback;
        
        if saveDataFlag
            
        end
    end

    function previousTrial_Callback(~,~)
        oldTrial = str2num(get(hTrialEdit,'string'));
        if oldTrial==1; return; end;
        newTrial = oldTrial-1;
        set(hTrialEdit,'string',newTrial);
        plotEyeData_Callback;
    end

    function nextTrial_Callback(~,~)
        oldTrial = str2num(get(hTrialEdit,'string'));
        if oldTrial==size(eyeDataDegX,1); return; end;
        newTrial = oldTrial+1;
        set(hTrialEdit,'string',newTrial);
        plotEyeData_Callback;
    end

    function plotEyeData_Callback(~,~)
        trial = str2num(get(hTrialEdit,'string'));
        subplot(plotEyeXY); cla
        subplot(plotEyeVelXY); cla
        subplot(plotEyePos); cla
        subplot(plotEyeVel); cla
        subplot(plotEyeSpeed); cla
        for iPeriod = 1:4
            switch iPeriod
                case 1
                    tVals = true(1,length(timeValsEyePos));
                    tVals(find(timeValsEyeIndicesBL==true,1,'first'):end)=false;
                    lineColour = [0 0 0];
                case 2
                    tVals = timeValsEyeIndicesBL;
                    lineColour = [1 0 0];
                case 3
                    tVals = timeValsEyeIndicesStim;
                    lineColour = [0 0 1];
                case 4
                    tVals = true(1,length(timeValsEyePos));
                    tVals(1:find(timeValsEyeIndicesStim==true,1,'last'))=false;
                    lineColour = [0 0 0];
            end

            subplot(plotEyeXY); hold on;
            plot(eyeDataDegX(trial,tVals),eyeDataDegY(trial,tVals),'color',lineColour); 
            plot(0,0,'+r'); makeBox(gca,[-2.5 2.5],[-2.5 2.5],'m',1,'-','B');
            lims = 1.2*fixationWindowWidth/2;
            xlim([-lims lims]); ylim([-lims lims]);
            axis square; hold off;
            
            subplot(plotEyeVelXY); hold on;
            plot(eyeSpeedX(trial,tVals),eyeSpeedY(trial,tVals),'color',lineColour); 
            plot(0,0,'+r'); 
            xEllipse=cutoffX(trial)*cos(-pi:0.01:pi); yEllipse=cutoffY(trial)*sin(-pi:0.01:pi);
            plot(xEllipse,yEllipse,'m');            
            axis square; hold off;

            subplot(plotEyePos); hold on
            plot(timeValsEyePos(1,tVals),eyeDataDegX(trial,tVals),'-','color',lineColour); hold on;   
            plot(timeValsEyePos(1,tVals),eyeDataDegY(trial,tVals),':','color',lineColour); hold on;
            xlim([timeValsEyePos(1) timeValsEyePos(end)]); ylim([-1 1]); 

            subplot(plotEyeVel); hold on
            plot(timeValsEyePos(1,tVals),eyeSpeedX(trial,tVals),'-','color',lineColour); hold on;   
            plot(timeValsEyePos(1,tVals),eyeSpeedY(trial,tVals),':','color',lineColour); hold on;
            xlim([timeValsEyePos(1) timeValsEyePos(end)]); ylim('auto'); 
            
            subplot(plotEyeSpeed); 
            plot(timeValsEyePos(1,tVals),eyeSpeedMag(trial,tVals),'-','color',lineColour); hold on;
            plot(timeValsEyePos(1,tVals),cutoffX(trial)+zeros(1,sum(tVals)),'-','color','m'); hold on;
            plot(timeValsEyePos(1,tVals),cutoffY(trial)+zeros(1,sum(tVals)),':','color','m'); hold on;
            xlim([timeValsEyePos(1) timeValsEyePos(end)]);
            text(0.7,0.6,['cutoff = ' num2str(round(cutoffX(trial),2),4) ', ' num2str(round(cutoffY(trial),2),4)],'unit','normalized','parent',plotEyeSpeed);
        end
            
        if ~isempty(MS{trial})
            disp([num2str(trial) ': ' num2str(length(MS{trial})) ': ' num2str(MS{trial})]);
            s=find(logical(trialVals(trial,:)));
            subplot(plotEyeXY); hold on; plot(eyeDataDegX(trial,s),eyeDataDegY(trial,s),'.','color','k','linewidth',2); hold on;
            subplot(plotEyeVelXY); hold on; plot(eyeSpeedX(trial,s),eyeSpeedY(trial,s),'.','color','k','linewidth',2); hold on;
            subplot(plotEyePos); hold on; plot(timeValsEyePos(1,s),zeros(1,length(s)),'.','color','k','linewidth',2); hold on; 
            subplot(plotEyeVel); hold on; plot(timeValsEyePos(1,s),zeros(1,length(s)),'.','color','k','linewidth',2); hold on;
            subplot(plotEyeSpeed); hold on; plot(timeValsEyePos(1,s),zeros(1,length(s)),'.','color','k','linewidth',2); hold on;
        else
            disp(num2str(trial));
        end
        
        subplot(plotEyePos); set(gca,'xticklabel',[]); ylabel('deg');
        subplot(plotEyeVel); set(gca,'xticklabel',[]); ylabel('deg/s');
        subplot(plotEyeSpeed); xlabel('Time (s)'); ylabel('deg/s');
        subplot(plotEyeXY); xlabel('Deg'); ylabel('Deg');
        subplot(plotEyeVelXY); xlabel('Deg/s'); ylabel('Deg/s');

    end

    function animate_Callback(~,~)
        trial = str2num(get(hTrialEdit,'string'));
        stepSizeString = get(hAnimateSpeedDropdown,'string');
        stepSize = stepSizeString(get(hAnimateSpeedDropdown,'val'),:);
        stepSize(strfind(stepSize,'X'))=[];
        stepSize = str2num(stepSize);
        
        subplot(plotEyeXY); cla
        subplot(plotEyeVelXY); cla
        subplot(plotEyePos); cla
        subplot(plotEyeVel); cla        
        subplot(plotEyeSpeed); cla        
        
        subplot(plotEyeXY);
        plot(0,0,'*r'); makeBox(gca,[-2.5 2.5],[-2.5 2.5],'m',1,'-','B'); hold on;   
        
        subplot(plotEyeVelXY);
        xEllipse=cutoffX(trial)*cos(-pi:0.01:pi); yEllipse=cutoffY(trial)*sin(-pi:0.01:pi);
        plot(xEllipse,yEllipse,'m'); hold on;        
        
        subplot(plotEyeSpeed);
        plot(timeValsEyePos,cutoffX(trial)+zeros(1,length(timeValsEyePos)),'-','color','m'); hold on;
        plot(timeValsEyePos,cutoffY(trial)+zeros(1,length(timeValsEyePos)),':','color','m'); hold on;
        text(0.7,0.6,['cutoff = ' num2str(round(cutoffX(trial),2),4) ', ' num2str(round(cutoffY(trial),2),4)],'unit','normalized','parent',plotEyeSpeed);
        
        for s = (stepSize+1):stepSize:size(eyeDataDegX,2)
            if timeValsEyeIndicesStim(s)==0 && timeValsEyeIndicesBL(s)==0
                lineColour = [0 0 0];
            elseif timeValsEyeIndicesBL(s)==1
                lineColour = [1 0 0];
            elseif timeValsEyeIndicesStim(s)==1
                lineColour = [0 0 1];
            end

            subplot(plotEyeXY);
            plot(eyeDataDegX(trial,s-stepSize:s),eyeDataDegY(trial,s-stepSize:s),'-','color',lineColour); hold on;
            
            subplot(plotEyeVelXY);
            plot(eyeSpeedX(trial,s-stepSize:s),eyeSpeedY(trial,s-stepSize:s),'-','color',lineColour); hold on;
            
            subplot(plotEyePos); xlim([timeValsEyePos(1) timeValsEyePos(end)]); hold on
            plot(timeValsEyePos(1,s-stepSize:s),eyeDataDegX(trial,s-stepSize:s),'-','color',lineColour); hold on;   
            plot(timeValsEyePos(1,s-stepSize:s),eyeDataDegY(trial,s-stepSize:s),':','color',lineColour); hold on;   

            subplot(plotEyeVel); xlim([timeValsEyePos(1) timeValsEyePos(end)]); hold on
            plot(timeValsEyePos(1,s-stepSize:s),eyeSpeedX(trial,s-stepSize:s),'-','color',lineColour); hold on;   
            plot(timeValsEyePos(1,s-stepSize:s),eyeSpeedY(trial,s-stepSize:s),':','color',lineColour); hold on;   

            subplot(plotEyeSpeed); xlim([timeValsEyePos(1) timeValsEyePos(end)]); hold on
            plot(timeValsEyePos(1,s-stepSize:s),eyeSpeedMag(trial,s-stepSize:s),'-','color',lineColour); hold on;  

            clear MSIndic
            MSIndic = logical(trialVals(trial,s-stepSize:s));
            if any(MSIndic); 
                clear MSPos
                MSPos = s-stepSize:s; MSPos(MSIndic==false)=[];
                subplot(plotEyeXY); plot(eyeDataDegX(trial,MSPos),eyeDataDegY(trial,MSPos),'.','color','k','linewidth',2); hold on; 
                subplot(plotEyeVelXY); plot(eyeSpeedX(trial,MSPos),eyeSpeedY(trial,MSPos),'.','color','k','linewidth',2); hold on; 
                subplot(plotEyePos); plot(timeValsEyePos(1,MSPos),zeros(1,length(MSPos)),'.','color','k','linewidth',2); hold on;
                subplot(plotEyeVel); plot(timeValsEyePos(1,MSPos),zeros(1,length(MSPos)),'.','color','k','linewidth',2); hold on;
                subplot(plotEyeSpeed); plot(timeValsEyePos(1,MSPos),zeros(1,length(MSPos)),'.','color','k','linewidth',2); hold on;                  
            end

            drawnow;
        end
        
        subplot(plotEyePos); set(gca,'xticklabel',[]); ylabel('deg');
        subplot(plotEyeVel); set(gca,'xticklabel',[]); ylabel('deg/s');
        subplot(plotEyeSpeed); xlabel('Time (s)'); ylabel('deg/s');
        subplot(plotEyeXY); xlabel('Deg'); ylabel('Deg');
        subplot(plotEyeVelXY); xlabel('Deg/s'); ylabel('Deg/s');

    end
    
end

function [MS,MSLength,numMSInRange,eyeSpeedX,eyeSpeedY,cutoffX,cutoffY,eyeSpeedMag,MSAllData,peakVelocityAll,MSLengthAll,MSMagAll,trialVals,nMS,MSAllDataExtended]...
    = getMS(eyeDataDegX,eyeDataDegY,timeValsEyePos,FsEye,timeRange,threshold,minCutOff,minMSLength)
    clear MS MSLength numMSInRange eyeSpeedX eyeSpeedY cutoff eyeSpeedMag
    [MS,MSLength,numMSInRange,eyeSpeedX,eyeSpeedY,cutoffX,cutoffY,eyeSpeedMag] = getMicroSaccadesFromEyePositionData(eyeDataDegX,eyeDataDegY,timeValsEyePos,FsEye,timeRange,threshold,minCutOff,minMSLength);

    % Get amplitude and velocity of microsaccades
    clear MSAllData peakVelocityAll MSLengthAll MSMagAll trialVals nMS numTrials
    [MSAllData,~,~,peakVelocityAll,MSLengthAll,MSMagAll,trialVals,nMS,MSAllDataExtended] = getMSAmplitudeAndPeakVelocity(MS,MSLength,eyeSpeedMag,eyeDataDegX,eyeDataDegY,timeRange,timeValsEyePos,FsEye);
end
