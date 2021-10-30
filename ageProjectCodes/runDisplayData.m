% Displays data

clear; clc;

% Mandatory fixed options
% folderSourceString = 'Users/supratimray/Supratim/Projects/TLSAEEGProject';
folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
% stRange = [0.25 0.75];

% Choose one of these options
refType = 'bipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% All Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%
subProjectName = 'age';
goodSubjects = getGoodSubjectsProjectwise(subProjectName,1);
uniqueSubjectNames = getGoodFileNamesForSubjects(goodSubjects{1});
[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);

%%%%%%%%%%%%%%%%%%%% Choose One List from here %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% AgeProject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
healthyPos = strcmp(cdrList,'HV');
malePos = strcmp(genderList,'M');
femalePos = strcmp(genderList,'F');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose option by uncommenting it out and commenting others 

% Both males and females
ageGroup1Pos = (ageList<65) & healthyPos;
ageGroup2Pos = (ageList>=65) & healthyPos;
clear subjectNameList
subjectNameList{1} = uniqueSubjectNames(ageGroup1Pos); strList{1} = 'Mid-All';
subjectNameList{2} = uniqueSubjectNames(ageGroup2Pos); strList{2} = 'Old-All';

% % Males
% ageGroup1Pos = (ageList<65) & healthyPos & malePos;
% ageGroup2Pos = (ageList>=65) & healthyPos & malePos;
% clear subjectNameList
% subjectNameList{1} = uniqueSubjectNames(ageGroup1Pos); strList{1} = 'Mid-Males';
% subjectNameList{2} = uniqueSubjectNames(ageGroup2Pos); strList{2} = 'Old-Males';

% Females
% ageGroup1Pos = (ageList<65) & healthyPos & femalePos;
% ageGroup2Pos = (ageList>=65) & healthyPos & femalePos;
% clear subjectNameList
% subjectNameList{1} = uniqueSubjectNames(ageGroup1Pos); strList{1} = 'Mid-Females';
% subjectNameList{2} = uniqueSubjectNames(ageGroup2Pos); strList{2} = 'Old-Females';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stRange = [0.25 0.75];
gamma1Range = [20 34];
gamma2Range = [36 66];
displayAnalyzedData(pwd,subjectNameList,strList,projectName,refType,protocolType,stRange,removeMicroSaccadesFlag,gamma1Range,gamma2Range); % Save data in analyzedData