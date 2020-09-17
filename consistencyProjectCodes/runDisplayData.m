clear; clc;

% Mandatory fixed options
folderSourceString = '/Users/supratimray/Supratim/Projects/TLSAEEGProject';
%folderSourceString = 'C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
subProjectName = 'consistencyProject';
stRange = [0.25 0.75];

% Choose one of these options
refType = 'bipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjectsWithAnalyzableBlocks = getSubjectAndBlocksStatistics(folderSourceString,protocolType,subProjectName);
numGroups = length(subjectsWithAnalyzableBlocks);
ageList = cell(1,numGroups);
genderList = cell(1,numGroups);
expDateList = cell(1,numGroups);
dataForDisplay = cell(1,numGroups);
uniqueSubjectNumbers = cell(1,numGroups);
for i=1:numGroups
    uniqueSubjectNames = getGoodFileNamesForSubjects(subjectsWithAnalyzableBlocks{i});
    uSN = zeros(1,length(uniqueSubjectNames));
    for j=1:length(uniqueSubjectNames)
        uSN(j) = str2double(uniqueSubjectNames{j}(1:3));
    end
    uniqueSubjectNumbers{i} = uSN;
    [ageList{i},genderList{i},~,expDateList{i}] = getDemographicDetails(projectName,uniqueSubjectNames);
    dataForDisplay{i} = combineAnalyzedData(folderSourceString,uniqueSubjectNames,projectName,refType,protocolType,stRange);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Matched subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
commonSubjectNumbers = intersect(uniqueSubjectNumbers{1},uniqueSubjectNumbers{2}); % numGroups=2
numSubjects = length(commonSubjectNumbers);

timeValsTF = dataForDisplay{1}.timeValsTF;
freqValsTF = dataForDisplay{1}.freqValsTF;
freqVals = dataForDisplay{1}.freqVals;

ageListMatched = zeros(numGroups,numSubjects);
genderListMatched = cell(numGroups,numSubjects);
dateNumsMatched = zeros(numGroups,numSubjects);
dTFPowerDBAllSubjects = zeros(numGroups,numSubjects,length(timeValsTF),length(freqValsTF));
dPowerVsFreqAllSubjects = zeros(numGroups,numSubjects,length(freqVals));
powerDBAllSubjects = zeros(numGroups,numSubjects,size(dataForDisplay{1}.powerDBAllSubjects,2));
for i=1:numGroups
    uSN = uniqueSubjectNumbers{i};
    
    for j=1:numSubjects
        pos = find(commonSubjectNumbers(j)==uSN);
        
        ageListMatched(i,j) = ageList{i}(pos);
        genderListMatched(i,j) = genderList{i}(pos);
        x = expDateList{i}{pos};
        dateNumsMatched(i,j) = datenum([x(1:2) '/' x(3:4) '/' x(5:6)],'dd/mm/yy');
        
        % dTF Power
        dTFPowerDBAllSubjects(i,j,:,:) = squeeze(dataForDisplay{i}.dTFPowerDBAllSubjects(:,:,pos));
        dPowerVsFreqAllSubjects(i,j,:) = squeeze(dataForDisplay{i}.logSTPowerVsFreqAllSubjects(pos,:) - dataForDisplay{i}.logBLPowerVsFreqAllSubjects(pos,:));
        powerDBAllSubjects(i,j,:) = squeeze(dataForDisplay{i}.powerDBAllSubjects(pos,:));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
powerForSorting = squeeze(mean(sum(powerDBAllSubjects(:,:,2),3),1)); % Avg for SG and FG
allMalePos = find(strcmp(genderListMatched(2,:),{'F'}));
[~,ids] = sort(powerForSorting(allMalePos),'descend');
displayDataIndividualSubjects(allMalePos(ids(1:10)),dTFPowerDBAllSubjects,dPowerVsFreqAllSubjects,timeValsTF,freqValsTF,freqVals);
