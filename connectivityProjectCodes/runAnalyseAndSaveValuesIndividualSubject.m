% Modified from the program written by Murty Dinavahi (MD)
% Modified by Santosh (16-Sep-2020)
% Data avalable in the decimatedData folder needs to be properly extracted
% for analysis. This program extracts the part that is used for final
% analysis and saves it under analyzedData.

clc; clear; 

% Mandatory fixed options
folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
subProjectName = 'ConnectivityProject';
stRange = [0.25 0.75];

% Choose one of these options
refType = 'unipolar'; % 'unipolar' or 'laplacian'
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.
spatialFrequenciesToRemove = 1;
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
        saveFTData(subjectName,expDates,protocolNames,capType,dataFolder,ftDataFolder,spatialFrequenciesToRemove);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Analyzed Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
connMethod = 'ppc'; % coh, plv
for iSub = 1:length(uniqueSubjectNames)
    subjectName = uniqueSubjectNames{iSub};
    disp([num2str(iSub) ': ' subjectName]);
    analyseAndSaveValuesIndividualSubjectConn(folderSourceString,subjectName,projectName,subProjectName,refType,protocolType,connMethod,stRange,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData);
end