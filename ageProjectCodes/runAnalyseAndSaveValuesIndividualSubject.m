% Modified from the program written by Murty Dinavahi (MD)

% Data avalable in the decimatedData folder needs to be properly extracted
% for analysis. This program extracts the part that is used for final
% analysis and saves it under analyzedData.

clear; clc;

% Mandatory fixed options
folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
stRange = [0.25 0.75];

% Choose one of these options
refType = 'bipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subProjectName = 'age';
goodSubjects = getGoodSubjectsProjectwise(subProjectName,1);
uniqueSubjectNames = getGoodFileNamesForSubjects(goodSubjects{1});

for iSub = 1:length(uniqueSubjectNames)
    subjectName = uniqueSubjectNames{iSub};
    disp([num2str(iSub) ': ' subjectName]);
    analyseAndSaveValuesIndividualSubject(folderSourceString,subjectName,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag); % Save data in analyzedData
end