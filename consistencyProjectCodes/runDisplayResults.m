clear; clc;

% Mandatory fixed options
folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
subProjectName = 'consistencyProject';
stRange = [0.25 0.75];

% Choose one of these options
refType = 'bipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjectsWithAnalyzableBlocks = getSubjectAndBlocksStatistics(folderSourceString,protocolType,subProjectName);

for i=1:length(subjectsWithAnalyzableBlocks)
    uniqueSubjectNames = getGoodFileNamesForSubjects(subjectsWithAnalyzableBlocks{i});
    [ageList(i,:),genderList{i},~,expDateList{i}] = getDemographicDetails(projectName,uniqueSubjectNames);
end

numSubjects = length(subjectsWithAnalyzableBlocks{1});
for i=1:numSubjects
    x = expDateList{1}{i}; y = expDateList{2}{i};
    diffDateNums(i) = datenum([y(1:2) '/' y(3:4) '/' y(5:6)],'dd/mm/yy') - datenum([x(1:2) '/' x(3:4) '/' x(5:6)],'dd/mm/yy');
end
