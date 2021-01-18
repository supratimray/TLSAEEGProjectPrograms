function [allProtocolsBLData,stPowerVsFreq,blPowerVsFreq,freqVals,tfPower,timeValsTF,freqValsTF,erp,timeVals,numGoodTrials,numAnalysedElecs]=getDataSingleSubject(cleanDataFolder,fileLists,capType,electrodeList,stRange,TFFlag,numMSInRangePerProtocol,condVals,params,discardBadElecFlag,spatialFrequenciesToRemove)

if ~exist('TFFlag','var') || isempty(TFFlag); TFFlag= 1; end
if ~exist('stRange','var') || isempty(stRange); stRange = [0.25 0.75]; end
if ~exist('params','var'); params=[]; end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end

if exist('numMSInRangePerProtocol','var') && ~isempty(numMSInRangePerProtocol)
    if length(numMSInRangePerProtocol)~=length(fileLists); error('Remove bad protocols first...'); end
end

iLoop = 1;
for iProt = 1:length(fileLists)
    sessionData = load(fullfile(cleanDataFolder,fileLists{iProt}),'eegData','badElecs','trialConditionVals','timeVals');
    eegData = sessionData.eegData;
    badElecs = sessionData.badElecs;
    trialConditionVals = sessionData.trialConditionVals;
    timeVals = sessionData.timeVals;
    if length(trialConditionVals)~=size(eegData,2)
        error('Wrong number of trials in EEG data or condition values...');
    end
    
    % Check for bad protocols
    badImpedanceElecs = badElecs.badImpedanceElecs;
    noisyElecs = badElecs.noisyElecs;
    flatPSDElecs = badElecs.flatPSDElecs;
    allBadElecs = unique([badImpedanceElecs;noisyElecs;flatPSDElecs])';
    
    electrodeListVis = getElectrodeList(capType,'bipolar');
    clear goodSideFlag
    for iSide = 1:length(electrodeListVis)
        clear goodElecFlag
        for iBipElec = 1:length(electrodeListVis{iSide})
            goodElecFlag(iBipElec) = ~any(ismember(electrodeListVis{iSide}{iBipElec},allBadElecs)); %#ok<*AGROW>
        end
        if any(goodElecFlag) 
            goodSideFlag(iSide) = true; 
        else
            goodSideFlag(iSide) = false; 
        end
    end
    
    % Deal with bad electrodes: discard data unless specified
    if exist('discardBadElecFlag','var') && ~discardBadElecFlag
        warning('Not discarding data from bad electrodes...')
    else
        for iBE = allBadElecs
            eegData(iBE,:,:) = NaN(size(eegData,2),size(eegData,3));
        end
    end
    
    % Remove microsaccde containing trials from analysis
    if exist('numMSInRangePerProtocol','var') && ~isempty(numMSInRangePerProtocol)
        msTrials = numMSInRangePerProtocol{iProt}>0;
        eegData(:,msTrials,:) = [];
        trialConditionVals(msTrials) = [];
    end
    
    if exist('condVals','var') && ~isempty(condVals)
        eegData(:,~ismember(trialConditionVals,condVals),:) = []; % This only works for TFCP protocols
    else
        if ~isempty(spatialFrequenciesToRemove)
            allSFs = trialConditionVals(:,1);
            badSFPos = [];
            for i=1:length(spatialFrequenciesToRemove)
                badSFPos = cat(1,badSFPos,find(spatialFrequenciesToRemove(i)==allSFs));
            end
            eegData(:,badSFPos,:)=[];
        end
    end

    if all(goodSideFlag)
        goodProtFlag(iLoop)=true; 
    else
        goodProtFlag(iLoop)=false;
    end
    
    [stPowerVsFreq(iLoop,:,:),blPowerVsFreq(iLoop,:,:),freqVals,tfPower(iLoop,:,:,:),timeValsTF,freqValsTF,erp(iLoop,:,:,:),numAnalysedElecs(iLoop,:)]=...
        getDataSingleSession(eegData,timeVals,electrodeList,stRange,TFFlag,params);
    numGoodTrials(iLoop) = size(eegData,2);
    iLoop = iLoop+1;
end

% Baseline data averaged across all eligible visual electrodes. Don't calculate for topoplots
if TFFlag && ~(exist('numMSInRangePerProtocol','var') && ~isempty(numMSInRangePerProtocol))
    allProtocolsBLData.blPowerVsFreq = removeDimIfSingleton(mean(blPowerVsFreq,2),2);
    allProtocolsBLData.goodProtFlag = goodProtFlag;
else
    allProtocolsBLData.goodProtFlag = goodProtFlag;
end

if any(goodProtFlag)
    stPowerVsFreq = removeDimIfSingleton(mean(stPowerVsFreq(goodProtFlag,:,:),1),1);
    blPowerVsFreq = removeDimIfSingleton(mean(blPowerVsFreq(goodProtFlag,:,:),1),1);
    tfPower = removeDimIfSingleton(mean(tfPower(goodProtFlag,:,:,:),1),1);
    erp = removeDimIfSingleton(mean(erp(goodProtFlag,:,:),1),1);
    numGoodTrials = sum(numGoodTrials(goodProtFlag));
else
    stPowerVsFreq = [];
    blPowerVsFreq = [];
    tfPower = [];
    erp = [];
    numGoodTrials = [];
end
end

function [stPowerVsFreq,blPowerVsFreq,freqVals,tfPower,timeValsTF,freqValsTF,erp,numAnalysedElecs]=getDataSingleSession(eegData,timeVals,electrodeList,stRange,TFFlag,params)

blRange = [-diff(stRange) 0];

% Get good positions
Fs = round(1/(timeVals(2)-timeVals(1)));
if round(diff(blRange)*Fs) ~= round(diff(stRange)*Fs)
    disp('baseline and stimulus ranges are not the same');
else
    range = blRange;
    rangePos = round(diff(range)*Fs);
    blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
    stPos = find(timeVals>=stRange(1),1)+ (1:rangePos);
end

% Set MT parameters
if ~exist('params','var') || isempty(params)
    params.tapers   = [1 1];
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 100];
    params.trialave = 1;
end

movingwin = [0.25 0.025];

% Initialize
numElectrodes = length(electrodeList{1});
numSides = length(electrodeList);

numAnalysedElecs = repmat(numElectrodes,numSides,1);

% hW = waitbar(0,'Analysing electrodes...'); iEs = 0;
for iElec=1:numElectrodes % For each electrode or electrode pair
    for iSide=1:numSides    % For each side
        
%         iEs = iEs+1;
%         hW = waitbar(iEs/(numElectrodes*numSides),hW,'Analysing electrodes...');
        
        if length(electrodeList{iSide}{iElec})==1 % Single Electrode
            eeg = squeeze(eegData(electrodeList{iSide}{iElec},:,:));
        else % for bipolar referencing
            chan1 = electrodeList{iSide}{iElec}(1);
            chan2 = electrodeList{iSide}{iElec}(2);
            eeg = squeeze(eegData(chan1,:,:) - eegData(chan2,:,:));
        end
        erp(iSide,iElec,:) = mean((eeg - repmat(mean(eeg(:,blPos),2),1,size(eeg,2))),1); % Correct for DC Shift (baseline correction)
        
        if TFFlag
            [tfPower(iSide,iElec,:,:),timeValsTF0,freqValsTF] = mtspecgramc(eeg',movingwin,params);
            timeValsTF = timeValsTF0 + timeVals(1);
        else
            tfPower = [];
            timeValsTF = [];
            freqValsTF = [];
        end
        stPowerVsFreq(iSide,iElec,:)= mtspectrumc(squeeze(eeg(:,stPos))',params);
        blPowerVsFreq(iSide,iElec,:)= mtspectrumc(squeeze(eeg(:,blPos))',params);
        
        if all(isnan(eeg(:))) % Discard bad electrodes
            erp(iSide,iElec,:) = NaN(1,size(erp,3));
            if TFFlag; tfPower(iSide,iElec,:,:) = NaN(size(tfPower,3),size(tfPower,4)); end
            stPowerVsFreq(iSide,iElec,:) = NaN(1,size(stPowerVsFreq,3));
            blPowerVsFreq(iSide,iElec,:) = NaN(1,size(blPowerVsFreq,3));
            numAnalysedElecs(iSide) = numAnalysedElecs(iSide)-1;
        end
    end
end
% close(hW);

% Find mean across electrodes
erp = removeDimIfSingleton(nanmean(erp,2),2);
tfPower = removeDimIfSingleton(nanmean(tfPower,2),2);
stPowerVsFreq = removeDimIfSingleton(nanmean(stPowerVsFreq,2),2);
blPowerVsFreq = removeDimIfSingleton(nanmean(blPowerVsFreq,2),2);

if ~isempty(eeg)
    [~,freqVals]= mtspectrumc(squeeze(eeg(1,stPos))',params);
else
    [~,freqVals]= mtspectrumc(squeeze(eeg(:,stPos))',params);
end

end