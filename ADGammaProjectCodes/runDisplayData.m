% Displays data

clear; clc;

% Mandatory fixed options
% folderSourceString = 'Users/supratimray/Supratim/Projects/TLSAEEGProject';
folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work

% Choose one of these options
refType = 'bipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.
spatialFrequenciesToRemove = 1;
useCleanData = 1; % cleanData refers to the data before decimation.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjectsWithAnalyzableBlocks = getSubjectAndBlocksStatistics(protocolType);
uniqueSubjectNames = getGoodFileNamesForSubjects(subjectsWithAnalyzableBlocks{1});

[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);

ageLim = 1;
healthyPos = strcmp(cdrList,'HV');
casePos = strcmp(cdrList,'MCI'); %~healthyPos;
caseList = setdiff(uniqueSubjectNames(casePos),[{'217SK'} {'225SK'}]);

controlList = []; controlListCaseNumber = [];
for i=1:length(caseList)
    subjectName = caseList{i};
    pos = find(strcmp(subjectName,uniqueSubjectNames));
    age = ageList(pos); gender = genderList(pos);
    
    ageMatchPos = (ageList<= age+ageLim) & (ageList>= age-ageLim);
    genderMatchPos = strcmp(gender,genderList);
    controlPos = healthyPos & ageMatchPos & genderMatchPos;
    controls = uniqueSubjectNames(controlPos);
    controlList = cat(2,controlList,controls);
    controlListCaseNumber = cat(2,controlListCaseNumber,i+zeros(1,length(controls)));
    disp([num2str(i) '. ' subjectName ' (' num2str(age) ',' gender{1} '): ' num2str(length(controls)) ' controls.']);
end

% Method 1
subjectNameList{1} = unique(controlList); strList{1} = 'Controls';
subjectNameList{2} = caseList; strList{2} = 'Cases';

% Method 2
casesWithControls = unique(controlListCaseNumber);
numValidCases = length(casesWithControls);

subjectNameListMatched = cell(1,2);
for i=1:numValidCases
    subjectNameListMatched{1}{i} = controlList(casesWithControls(i)==controlListCaseNumber);
    subjectNameListMatched{2}{i} = caseList(casesWithControls(i));
end

stRange = [0.25 0.75]; gamma1Range = [20 34]; gamma2Range = [36 66]; alphaRange = [8 12];
displayAnalyzedData(pwd,subjectNameListMatched,strList,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range,alphaRange,1,spatialFrequenciesToRemove,useCleanData); % Save data in analyzedData
