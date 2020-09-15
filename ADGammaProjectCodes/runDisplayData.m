% Displays data

clear; clc;

% Mandatory fixed options
folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
stRange = [0.25 0.75];

% Choose one of these options
refType = 'bipolar'; % 'unipolar' % Set reference type here.
protocolType = 'TFCP'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjectsWithAnalyzableBlocks = getSubjectAndBlocksStatistics(folderSourceString,protocolType);
uniqueSubjectNames = getGoodFileNamesForSubjects(subjectsWithAnalyzableBlocks);

[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);

%%%%%%%%%%%%%%%%%%%% Choose One List from here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% AgeProject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% healthyPos = strcmp(cdrList,'HV');
% malePos = strcmp(genderList,'M');
% femalePos = strcmp(genderList,'F');
% 
% ageGroup1Pos = (ageList<65) & healthyPos; % ageGroup1Pos = (ageList<65) & healthyPos & femalePos;
% ageGroup2Pos = (ageList>=65) & healthyPos; % ageGroup2Pos = (ageList>=65) & healthyPos & femalePos;
% 
% clear subjectNameList
% subjectNameList{1} = uniqueSubjectNames(ageGroup1Pos); strList{1} = 'Mid';
% subjectNameList{2} = uniqueSubjectNames(ageGroup2Pos); strList{2} = 'Old';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADGammaProject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ageLim = 1;
healthyPos = strcmp(cdrList,'HV');
casePos = ~healthyPos; % casePos = strcmp(cdrList,'MCI');

caseList = uniqueSubjectNames(casePos);
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

%getSubjectAndBlocksStatistics(folderSourceString,protocolType);
displayAnalyzedData(folderSourceString,subjectNameList,strList,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag); % Save data in analyzedData