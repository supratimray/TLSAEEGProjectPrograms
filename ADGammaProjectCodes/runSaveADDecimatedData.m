% This program simply saves the DecimatedData for the ADGammaProject cases
% and control subjects only

clear; clc;

% Mandatory fixed options
folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
folderDestinationString = 'C:\Users\Supratim Ray\Desktop';

projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP

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

dataFolderIn = fullfile(folderSourceString,'decimatedData',projectName,protocolType);
dataFolderOut = fullfile(folderDestinationString,'decimatedData',projectName,protocolType);
makeDirectory(dataFolderOut);

for iSub = 1:length(goodSubjectsAll)
    subjectName = goodSubjectsAll{iSub};
    disp([num2str(iSub) ': ' subjectName]);

    [expDates,protocolNames,capLayout,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);

    if usableDataFlag && ~isempty(expDates)
        for iProt = 1:length(expDates)
            fileNameToCopy = fullfile(dataFolderIn,[subjectName '-' expDates{iProt} '-' protocolNames{iProt} '.mat']);
            fileNameDestination = fullfile(dataFolderOut,[subjectName '-' expDates{iProt} '-' protocolNames{iProt} '.mat']);
            copyfile(fileNameToCopy,fileNameDestination);
        end
    end
end