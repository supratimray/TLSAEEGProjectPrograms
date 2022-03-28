% This program gets data in a format that is used in fieldtrip.
% Written by W Santosh. Modified by Supratim Ray.

% Ideally, whether input data is cleanData or decimatedData should be
% recorded somewhere. Also, which SFs are removed should be recorded. But
% that would increase the ftData folder size a lot. To keep things simple,
% we are not recording this information while generating ftData. It is
% meant to be a temporary folder, which should be deleted and regenerated
% when these parameters need to be changed.

function saveFTData(subjectName,expDates,protocolNames,capType,dataFolder,ftDataFolder,spatialFrequenciesToRemove)

elec = getElectrodeDetails_ft(capType);

if ~isempty(elec)
    for i=1:length(expDates)
        fileName = [subjectName '-' expDates{i} '-' protocolNames{i} '.mat'];
        dataFileName = fullfile(dataFolder,fileName);
        ftDataFileName = fullfile(ftDataFolder,fileName);
        [data,goodProtFlag]=getFTDataType(dataFileName,capType,spatialFrequenciesToRemove,elec);
        numGoodTrials = data.numTrials;
        save(ftDataFileName,'data','goodProtFlag','numGoodTrials');
    end
end
end
function [data,goodProtFlag] = getFTDataType(dataFileName,capType,spatialFrequenciesToRemove,elec)
% returns 'data' of the configuration used by fieldtrip functions
% BAD trials are removed so as to keep these file size to minimum
% trials corresponding to spatialFrequenciesToRemove are also removed
% Ideally, whether input data is cleanData or decimatedData should be
% recorded somewhere. Also, which SFs are removed should be recorded. But
% that would increase the ftData folder size a lot. To keep things simple,
% we are not recording this information. 

inputData = load(dataFileName);
% inputData has fields
% -badElecs
% -eegData
% -eyeData
% -eyeRangeMS
% -timeVals
% -timeValsEye
% -trialConditionLabels
% -trialConditionVals

timeVals=inputData.timeVals;
numElecs=size(inputData.eegData,1);

sig_all_in = inputData.eegData; % size: numOfElecs X numTrials X length(timeVals)

% applying powerline correction % DO THIS IN NEXT VERSION
% sig_all = powerline_correct(sig_all_in);
sig_all = sig_all_in;

if ~isempty(spatialFrequenciesToRemove)
    allSFs = inputData.trialConditionVals(:,1);
    badSFPos = [];
    for i=1:length(spatialFrequenciesToRemove)
        badSFPos = cat(1,badSFPos,find(spatialFrequenciesToRemove(i)==allSFs));
    end
    sig_all(:,badSFPos,:)=[];
end
% update numTrials
numTrials = size(sig_all,2);

eeg_all=cell(1,numTrials);
for tr=1:numTrials
    eeg_all{tr}=squeeze(sig_all(:,tr,:));
end

chantype=cell(numElecs,1);
chanunit=cell(numElecs,1);
for i=1:numElecs
    chantype{i}='eeg';
    chanunit{i}='uV';
end

hdr.Fs=round(1/(timeVals(2)-timeVals(1)));
hdr.nChans=numElecs;
hdr.nSamples = length(timeVals); % samples in each trial
hdr.nTrials=numTrials;
hdr.label=elec.label;
hdr.chantype=chantype;
hdr.chanunit=chanunit;

times=cell(1,numTrials);
for i=1:length(times)
    times{i}=timeVals;
end

data.hdr=hdr;
data.label=hdr.label;
data.trial=eeg_all; % cell of each trial as (chans x timpts)
data.time=times; % cell of each trial timepoints
data.fsample=hdr.Fs;
data.cfg=[];
data.elec=elec;

data.timeVals = timeVals;
data.numTrials = numTrials;
data.trialConditionLabels = inputData.trialConditionLabels;
data.trialConditionVals = inputData.trialConditionVals;

data.badElecs = unique([inputData.badElecs.badImpedanceElecs; inputData.badElecs.flatPSDElecs; inputData.badElecs.noisyElecs]);

% goodProtFlag
electrodeListVisFull = getElectrodeList(capType,'unipolar',0,1);
electrodeListVis = electrodeListVisFull(2:4);
clear goodSideFlag goodElecFlag
for iSide = 1:length(electrodeListVis)
    for iUniElec = 1:length(electrodeListVis{iSide})
        goodElecFlag(iUniElec) = ~any(ismember(electrodeListVis{iSide}{iUniElec},data.badElecs)); %#ok<*AGROW>
    end
    if any(goodElecFlag)
        goodSideFlag(iSide) = true;
    else
        goodSideFlag(iSide) = false;
    end
end
goodProtFlag=all(goodSideFlag);
end
function elec = getElectrodeDetails_ft(capType)

if strcmp(capType,'actiCap64')
    labels = load('actiCap64Labels.mat');
    chanlocs = load('actiCap64.mat'); %#ok<NASGU>
    
    dimLabels = 'XYZ';
    for dim = 1:3
        for e = 1:64
            clocs(e,dim) = eval(['chanlocs.chanlocs(' int2str(e) ').' dimLabels(dim)]);
        end
    end
    elec.chanpos = clocs;
    elec.chantype = repmat({'eeg'},64,1);
    elec.elecpos = clocs;
    elec.label = labels.montageLabels(:,2);
    elec.type = 'eeg1010';
    elec.unit = 'cm';
else
    disp([capType{1} ' not specified']);
    elec = [];
end
end