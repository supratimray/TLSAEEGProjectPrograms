
% Program to store decimated data from cleanData file for each protocol for
% each subject. The file format for saving is:
% xDataParentFolder/xData/projectName/protocolType/subjectName-expDate-protocolName.mat
% where x is either 'clean' or 'decimated' as the case may be
%
% Murty V P S Dinavahi: 21-May-2020

% Note: For two files (N001-050716-GAV-0001 and -0003), the program crashes
% because these names are not in the subjectNamesAll. These files
% correspond to '049VK (N0001)'. We don't worry about it because this
% dataset is anyway not used for future analysis since eye data was not
% collected for these two datasets.


clear; clc;
projectName = 'ADGammaProject';%'ADGammaProject';%'DualGamma';%AgeProjectRound1
protocolTypes = {'SF_ORI' 'TFCP'};
if strcmp(projectName,'ADGammaProject')
    cleanDataParentFolder = 'N:\Projects\Murty_ADGammaProject'; % Parent folder for storing the cleanData Folder
elseif strcmp(projectName,'DualGammaProject')
    cleanDataParentFolder = 'N:\Projects\Done\Murty_DualGammaProject';
elseif strcmp(projectName,'AgeProjectRound1')
    cleanDataParentFolder = 'N:\Projects\Murty_AgeProject';
end

% decimatedDataParentFolder = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\TataADProject'; % Parent folder for storing the decimatedData Folder.
decimatedDataParentFolder = 'N:\Projects\TataADProject';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

decimationFactor = 10; % decimate by a factor of 10 (should be a factor for original Fs)

% The current cleanData file doesn't have timeVals for eye data. Adding it in
% this code, to be used if required later.
if strcmpi(projectName,'DualGammaProject')
    [subjectNamesAll,expDatesAll,protocolNamesAll,~,~,~,~,FsEyeAll] = eval('allProtocolsVisualGammaEEG');
else
    [subjectNamesAll,expDatesAll,protocolNamesAll,~,~,~,~,~,FsEyeAll] = loadProjectDetails(projectName);       
end

for iProt = 1:length(protocolTypes)
    protocolType = protocolTypes{iProt};
    cleanDataFolder = fullfile(cleanDataParentFolder,'cleanData',protocolType);
    decimatedDataFolder = fullfile(decimatedDataParentFolder,'decimatedData',projectName,protocolType);
    makeDirectory(decimatedDataFolder);
    listND = dir(cleanDataFolder);

    for iL = 3:length(listND)
        disp(num2str(iL));
        clear dataStruct; dataStruct = load(fullfile(cleanDataFolder,listND(iL).name));
        
        %%%% Decimate EEG data %%%%
        clear eegData eegDataDecimated
        eegData = dataStruct.eegData;        
        for iElec = 1:size(eegData,1)
            eegDataDecimated(iElec,:,:) = resample(squeeze(eegData(iElec,:,:))',1,decimationFactor)'; %#ok<SAGROW>
        end
        dataStruct.eegData = eegDataDecimated;
        dataStruct.timeVals = downsample(dataStruct.timeVals',decimationFactor)';
        
        %%%% Add trialCondition labels %%%%
        if strcmpi(protocolType,'SF_ORI'); trialConditionLabels = {'SF' 'ORI'};
        else trialConditionLabels = protocolType;
        end
        dataStruct.trialConditionLabels = trialConditionLabels;
        
        %%%% Get timeVals for eyeData %%%%
        % Get FsEye
        fileDetails = strsplit(listND(iL).name,'-');
        subjectName = fileDetails{1};
        expDate = fileDetails{2};
        protocolName = fileDetails{3}(1:end-4);
        indeX = strcmpi(subjectNamesAll,subjectName) & strcmpi(expDatesAll,expDate) & strcmpi(protocolNamesAll,protocolName);        
        FsEye = FsEyeAll{indeX}; if ischar(FsEye); if ~strcmp(FsEye,'-'); FsEye = str2num(FsEye); else FsEye = []; end; end;
        
        % get timeValsEye
        if ~isempty(FsEye)
            eyeRangeMS = dataStruct.eyeRangeMS;
            timeValsEye = (eyeRangeMS(1):1000/FsEye:eyeRangeMS(2)-1000/FsEye)/1000;
        else
            timeValsEye = [];
        end        
        dataStruct.timeValsEye = timeValsEye;
        
        % Save
        save(fullfile(decimatedDataFolder,listND(iL).name),'-struct','dataStruct')       
    end
end