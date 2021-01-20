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
spatialFrequenciesToRemove = 1;
useCleanData = 1; % cleanData refers to the data before decimation. You must have cleanData folder for this option to work.
temporalFrequencyToUse = 0; % 0 or 16; Only used for TFCP protocol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjectsList = getGoodSubjectsProjectwise(projectName);
uniqueSubjectNames = getGoodFileNamesForSubjects(goodSubjectsList{1});
[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);

ageLim = 1;
healthyPos = strcmp(cdrList,'HV');
casePos = ~healthyPos;
caseList = setdiff(uniqueSubjectNames(casePos),[{'217SK'} {'225SK'}]);

controlList = [];
for i=1:length(caseList)
    subjectName = caseList{i};
    pos = find(strcmp(subjectName,uniqueSubjectNames));
    age = ageList(pos); gender = genderList(pos);
    
    ageMatchPos = (ageList<= age+ageLim) & (ageList>= age-ageLim);
    genderMatchPos = strcmp(gender,genderList);
    controlPos = healthyPos & ageMatchPos & genderMatchPos;
    controls = uniqueSubjectNames(controlPos);
    controlList = cat(2,controlList,controls);
    
    disp([num2str(i) '. ' subjectName ' (' num2str(age) ',' gender{1} '): ' num2str(length(controls)) ' controls.']);
end

subjectNameList{1} = unique(controlList); strList{1} = 'Controls';
subjectNameList{2} = caseList; strList{2} = 'Cases';

goodSubjectsAll = [subjectNameList{1} subjectNameList{2}];

for iSub = 1:length(goodSubjectsAll)
    subjectName = goodSubjectsAll{iSub};
    disp([num2str(iSub) ': ' subjectName]);
    analyseAndSaveValuesIndividualSubject(folderSourceString,subjectName,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData,temporalFrequencyToUse); % Save data in analyzedData
end
