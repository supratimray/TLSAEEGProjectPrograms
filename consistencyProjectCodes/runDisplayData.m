% This program uses the common exlcusion criteria that is applicable for
% all projects. However, additional exclusion criteria were added in the
% paper later, which reduced the number of subjects from 48 to 40. These
% can be found here: https://github.com/wupadrasta/TLSAEEGProjectPrograms.
% To preserve compatibility with other projects, these criteria are not
% used in this program, so we end up with 48 subjects.

clear; clc;

% Mandatory fixed options
% folderSourceString = '/Users/supratimray/Supratim/Projects/TLSAEEGProject';
folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\TLSAEEGProject'; % Indicate the parent folder of decimatedData
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
subProjectName = 'consistencyProject';
stRange = [0.25 0.75];

% Choose one of these options
refType = 'bipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP
removeMicroSaccadesFlag = 0; % 0 or 1.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subjectsWithAnalyzableBlocks = getSubjectAndBlocksStatistics(protocolType,subProjectName);
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
    dataForDisplay{i} = combineAnalyzedData(pwd,uniqueSubjectNames,projectName,refType,protocolType,stRange);
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
numSubjectsPerPlot = 20;

% Females
powerForSorting = squeeze(mean(sum(powerDBAllSubjects(:,:,2),3),1)); % Avg for SG and FG
goodPos = find(strcmp(genderListMatched(2,:),{'F'})); % There is an error in the gender of 1 subject. Using list 2 which has the correct gender.
% goodPos = find(strcmp(genderListMatched(2,:),{'M'}));

[~,ids] = sort(powerForSorting(goodPos),'descend');
displayDataIndividualSubjects(goodPos(ids(1:numSubjectsPerPlot)),dTFPowerDBAllSubjects,dPowerVsFreqAllSubjects,timeValsTF,freqValsTF,freqVals);

% %%%%%%%%%%%%%%%%%%%%%%  Self vs Other Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:numSubjects
%     x = squeeze(dPowerVsFreqAllSubjects(1,i,:));
%     for j=1:numSubjects
%         y = squeeze(dPowerVsFreqAllSubjects(2,j,:));
%         c = corrcoef(x,y);
%         cData(i,j) = c(1,2);
%     end
%     selfCorr(i) = cData(i,i);
%     otherCorr(i) = median(setdiff(cData(i,:),cData(i,i)));
% end
% plot(selfCorr,otherCorr,'o');