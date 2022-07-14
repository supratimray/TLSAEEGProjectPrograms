clc; clear;
% Mandatory fixed options
% folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
folderSourceString = 'E:\Santosh\Project codes\TataADProject'; % of decimated data
projectName = 'ADGammaProject';
subProjectName = 'ConnectivityProject';
stRange = [0.25 0.75]; 
pow_label = 'diff'; % 'st' or 'diff' or 'bl' [used for power-matching]
numControls = []; % for 'unmatched' case in CaseVsControl, [] indicates all the available controls
combineOppSide = true;
newfigpath = 'E:\Santosh\Project codes\TataADProject\Plots\AgeConnectivity\NEWnewPlots';
medianFlag = true;

freqRanges{1} = [8 12]; freqRangeNames{1} = 'Alpha'; % alpha
freqRanges{2} = [20 34]; freqRangeNames{2} = 'Slow gamma'; % slow gamma
freqRanges{3} = [36 66]; freqRangeNames{3} = 'Fast gamma'; % fast gamma
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
healthyPos = strcmp(cdrList,'HV');
caseType = 'MCI';
casePos = strcmp(cdrList,caseType); %strcmp(cdrList,'AD'); % 'MCI', 'AD' or ~healthyPos for MCI+AD
methodOptions.caseType = caseType;
methodOptions.CaseVsControl.ageList = ageList(casePos);
methodOptions.MidVsOld.ageList = ageList(healthyPos);
methodOptions.CaseVsControl.unmatched.numControls = numControls;
methodOptions.CaseVsControl.genderList = genderList(casePos);
%%%%%%%%%%%%%%%%%% Generates subjectNameListFinal %%%%%%%%%%%%%%%%%%%%%%%%%
% Note that subjectNameListFinal is a cell array of size 2 which has different structure depending on comparisonType.
% For MidVsOld - the cell arrays contain the subjectNames in middled aged and elderly groups
% For CaseVsControl - each element of the cell array is itself a cell array
% containing all subjectNames that need to be averaged.
% Power matching is done within the display Program
clear subjectNameListFinal
SFs = {[], 1}; % SFs to remove
comparisonTypes = {'MidVsOld', 'CaseVsControl'};
connMethods = {'coh','plv','ppc'};
powTypes = {'unmatched', 'matched'};
methodOptions.pow_label = pow_label;
for iSF = 1
    spatialFrequenciesToRemove = SFs{iSF};
    useCleanData = 0;
    for compCond = 1
        methodOptions.comparisonType = comparisonTypes{compCond}; % MidVsOld or 'CaseVsControl'
        if strcmp(methodOptions.comparisonType,'MidVsOld')
            % For MidVsOld - we simply take healthy subjects who are less than 65, combining both males and females
            ageGroup1Pos = (ageList<65) & healthyPos;
            ageGroup2Pos = (ageList>=65) & healthyPos;
            subjectNameListFinal{1} = uniqueSubjectNames(ageGroup1Pos); strList{1} = 'MiddleAged';
            subjectNameListFinal{2} = uniqueSubjectNames(ageGroup2Pos); strList{2} = 'Elderly';
        elseif strcmp(methodOptions.comparisonType,'CaseVsControl')
            ageLim = 1;
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
            
            casesWithControls = unique(controlListCaseNumber);
            numValidCases = length(casesWithControls);
            subjectNameListFinal = cell(1,2);
            for i=1:numValidCases
                subjectNameListFinal{1}{i} = controlList(casesWithControls(i)==controlListCaseNumber);
                subjectNameListFinal{2}{i} = caseList(casesWithControls(i));
            end
            strList{1} = 'Controls';
            strList{2} = 'Cases';
        end
        
        useMedianFlag = 1;
        for connCond = 1%:length(connMethods)
            methodOptions.connMethod = connMethods{connCond};
            for powCond = 2%1:2
                methodOptions.powcontrolType = powTypes{powCond};
                for sideToShow = 1%:3
                    methodOptions.sideToShow = sideToShow;
                    methodOptions.combineOppSide = combineOppSide; % combining left & right sides
                    figR = figure('numbertitle','off','name',generateFigName(methodOptions,spatialFrequenciesToRemove));
                    figR.PaperType = 'a4';
                    figR.PaperUnits = 'centimeters';
                    figR.PaperSize = [18 29.7];
                    figR.PaperOrientation = 'Landscape';
                    figR.PaperPosition = [0 0 figR.PaperSize];
                    figR.Color = [1 1 1]; % White background
                    displayAnalyzedDataConn(pwd,subjectNameListFinal,methodOptions,strList,subProjectName,refType,protocolType,stRange,freqRanges,freqRangeNames,removeMicroSaccadesFlag,useMedianFlag,spatialFrequenciesToRemove,useCleanData,medianFlag);
                    print(figR,'-painters',fullfile(newfigpath,generateFigName(methodOptions,spatialFrequenciesToRemove)),'-dtiff','-r300');
                end
            end
        end
    end
end

function figName = generateFigName(methodOptions,spatialFrequenciesToRemove)
ElecGLabels = {'Left', 'Right', 'Back'};
figName = [];
figName = cat(2,figName,[methodOptions.comparisonType '-' methodOptions.powcontrolType '-' methodOptions.pow_label]);
figName = cat(2,figName,['-' methodOptions.connMethod]);
figName = cat(2,figName,['_' ElecGLabels{methodOptions.sideToShow}]);
figName = cat(2,figName,'_RemoveSF[');
for i=1:length(spatialFrequenciesToRemove)
    figName = cat(2,figName,num2str(spatialFrequenciesToRemove(i)));
end
figName = cat(2,figName,']');
end