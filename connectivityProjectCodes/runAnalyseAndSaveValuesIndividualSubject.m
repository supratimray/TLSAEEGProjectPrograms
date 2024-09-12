% Modified from the program written by Murty Dinavahi (MD)
% Modified by Santosh (16-Sep-2020)
% Data avalable in the decimatedData folder needs to be properly extracted
% for analysis. This program extracts the part that is used for final
% analysis and saves it under analyzedData.

clc; clear; 

% Mandatory fixed options
folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
% folderSourceString = 'E:\Santosh\Project codes\TataADProject'; % of decimated data
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
subProjectName = 'ConnectivityProject';
stRange = [0.25 0.75];

% Choose one of these options
refType = 'unipolar'; % 'unipolar' or 'laplacian'
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.
spatialFrequenciesToRemove = []; % This is not changed, as others are implemented later
useCleanData = 0;

%%%%%%%%%%%%%%%%%%%%%%%%% Get Good Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjects = getGoodSubjectsProjectwise(projectName,1); % program will consider current projectName as ADGammaProject. No diff.
uniqueSubjectNames = getGoodFileNamesForSubjects(goodSubjects{1});

%%%%%%%%%%%%%%%%%%%%%%%%%% Save FieldTrip Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFTDataFlag=0; % if set to 1, data in FT format is saved. 
if saveFTDataFlag
    % Input Data Folder
    if ~useCleanData %#ok<*UNRCH>
        dataFolder = fullfile(folderSourceString,'decimatedData',projectName,protocolType); % Decimated data - decimated version of cleanData
    else
        dataFolder = fullfile(folderSourceString,'cleanData',protocolType);
    end
    
    % Output Fieldtrip Folder
    ftDataFolder = fullfile(folderSourceString,'ftData',projectName,protocolType);
    makeDirectory(ftDataFolder);
    
    for iSub = 1:length(uniqueSubjectNames)
        subjectName = uniqueSubjectNames{iSub};
        [expDates,protocolNames,capType,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);
        
        disp([num2str(iSub) ': ' subjectName]);
        saveFTData(subjectName,expDates,protocolNames,capType,dataFolder,ftDataFolder,spatialFrequenciesToRemove); % option "spatialFrequenciesToRemove" is not being used here
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Save Analyzed Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
SFs = {[], 1}; % this is implemented on the data after loading (no separate files)
connMethods = {'coh','plv','ppc'};
for connCond = 3 %1:length(connMethods)
    for iSF = 1 %1:2
        spatialFrequenciesToRemove = SFs{iSF};
        connMethod = connMethods{connCond}; % coh, plv, ppc
        wb = waitbar(0,['Computing connectivity.. ' connMethod]);
        for iSub = 1:length(uniqueSubjectNames)
            waitbar(iSub/length(uniqueSubjectNames),wb,['Processing subj: ' int2str(iSub)]);
            subjectName = uniqueSubjectNames{iSub};
            disp([num2str(iSub) ': ' subjectName]);
            analyseAndSaveValuesIndividualSubjectConn(folderSourceString,subjectName,projectName,subProjectName,refType,protocolType,connMethod,stRange,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
        end
        close(wb);
    end
end