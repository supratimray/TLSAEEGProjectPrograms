% Modified from runDisplayData. Code simplified to plot connectivity as a
% function of frequency

useBLConnData = 1;
comparisonType = 'MidVsOld'; % 'CaseVsControl' or 'MidVsOld';

% Mandatory fixed options
projectName = 'ADGammaProject';
subProjectName = 'ConnectivityProject';
stRange = [0.25 0.75]; 

useMedianFlagData = true;

freqRanges{1} = [8 12]; freqRangeNames{1} = 'Alpha'; % alpha
freqRanges{2} = [28 34]; freqRangeNames{2} = 'Slow gamma'; % slow gamma
freqRanges{3} = [40 66]; freqRangeNames{3} = 'Fast gamma'; % fast gamma
numFreqRanges = length(freqRanges);

%%%%%%%%%%%%%%%%%%%%%%%% Choose one of these options %%%%%%%%%%%%%%%%%%%%%%
refType = 'unipolar'; % 'unipolar' or 'laplacian'
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Get Good Subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

goodSubjects = getGoodSubjectsProjectwise(projectName,1);
uniqueSubjectNames0 = getGoodFileNamesForSubjects(goodSubjects{1});
%%%%%%%%%%%%%% Find indices for which the correct capType was used %%%%%%%%
capTypeToUse = 'actiCap64';
goodIndices = [];
for i=1:length(uniqueSubjectNames0)
    [expDates,~,capType,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,uniqueSubjectNames0{i},protocolType);
    if usableDataFlag && ~isempty(expDates) && strcmp(capType{1},capTypeToUse)
        goodIndices = cat(2,goodIndices,i);
    end
end
disp([num2str(length(goodIndices)) ' subjects with correct capType chosen for further analysis']);
uniqueSubjectNames = uniqueSubjectNames0(goodIndices);
[ageList,genderList,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);

%%%%%%%%%%%%%%%%%%%%%% Use this for mid vs early %%%%%%%%%%%%%%%%%%%%%%%%%%
healthyPos = strcmp(cdrList,'HV');

if strcmp(comparisonType,'MidVsOld')
    ageGroup1Pos = (ageList<65) & healthyPos;
    ageGroup2Pos = (ageList>=65) & healthyPos;
    subjectNameListFinal{1} = uniqueSubjectNames(ageGroup1Pos); strList{1} = 'MiddleAged';
    subjectNameListFinal{2} = uniqueSubjectNames(ageGroup2Pos); strList{2} = 'Elderly';

else
    ageLim = 1;
    caseType = 'MCI';
    casePos = strcmp(cdrList,caseType);
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
    subjectNameList{1} = unique(controlList); strList{1} = 'Controls';
    subjectNameList{2} = caseList; strList{2} = 'Cases';

    option = 1;
    if option==1 % Unmatched case
        subjectNameListFinal = subjectNameList;
    
    elseif option==2 
        casesWithControls = unique(controlListCaseNumber);
        numValidCases = length(casesWithControls);
        subjectNameListFinal = cell(1,2);

        % This will not work for now since the subjectNameListFinal has to be a cell array of strings. 
        for i=1:numValidCases
            subjectNameListFinal{1}{i} = controlList(casesWithControls(i)==controlListCaseNumber);
            subjectNameListFinal{2}{i} = caseList(casesWithControls(i));
        end

        % Get a single control per case
        subjectNameListMatched = getNewControlSubjectList(projectName,subjectNameListFinal,1);
        subjectNameListFinal = cell(1,2);

        for i=1:2
            x = [];
            for j=1:length(subjectNameListMatched{i})
                x = cat(2,x,subjectNameListMatched{i}{j});
            end
            subjectNameListFinal{i} = x;
        end
    end
end

folderSourceString = pwd;
displayAnalyzedDataConnVsFreq(folderSourceString,subjectNameListFinal,strList,subProjectName,refType,protocolType,stRange,freqRanges,freqRangeNames,removeMicroSaccadesFlag,[],0,useMedianFlagData,useBLConnData);