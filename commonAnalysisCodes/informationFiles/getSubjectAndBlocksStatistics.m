% This program provides details of the number of subjects and usable
% blocks, using data stored in analysedData.

function subjectsWithAnalyzableBlocks = getSubjectAndBlocksStatistics(folderSourceString,protocolType,subProjectName)

if ~exist('protocolType','var');     protocolType = 'SF_ORI';           end

% Mandatory options
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work

if ~exist('subProjectName','var')
    subProjectName = projectName;
end

analyzedDataFolder = fullfile(folderSourceString,'analyzedData',projectName,protocolType);

subjectNamesList = getGoodSubjectsProjectwise(subProjectName,1,protocolType);
numSubjectList = length(subjectNamesList);
subjectsWithAnalyzableBlocks = cell(1,numSubjectList);
for i=1:numSubjectList
    subjectNames = getGoodFileNamesForSubjects(subjectNamesList{i});
    [allProts,allSubjectIDs] = getProtsAndSubjectIDs(subjectNames,analyzedDataFolder);

    disp(['Protocol Type: ' protocolType]);
    disp(['Total number of subjects: ' num2str(length(subjectNames))]);
    disp(['Total number of blocks: ' num2str(length(allProts))]);
    disp(['Number of rejected blocks: ' num2str(length(find(allProts==0)))]);
    
    subjectsWithNoAnalyzableBlocks = setdiff(subjectNames,subjectNames(allSubjectIDs(allProts==1)));
    subjectsWithAnalyzableBlocks{i} = setdiff(subjectNames,subjectsWithNoAnalyzableBlocks);
    disp(['Number of good subjects: ' num2str(length(subjectsWithAnalyzableBlocks{i}))]);
    describeSessionDetails(subjectsWithAnalyzableBlocks{i},projectName);
end
end

function [allProts,allSubjectIDs,allNumTrials] = getProtsAndSubjectIDs(subjectNames,analyzedDataFolder)

stRange = [0.25 0.75];
refType = 'bipolar'; % 'unipolar' % Set reference type here.

allProts = [];
allSubjectIDs = [];
allNumTrials = [];

for iSub = 1:length(subjectNames)
    subjectName = subjectNames{iSub};
    
    analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
        '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '.mat']);
    
    if exist(analysisDetailsFile,'file')
        x=load(analysisDetailsFile);
        p = [x.allProtocolsBLData.goodProtFlag];
        allProts=cat(2,allProts,p);
        allSubjectIDs = cat(2,allSubjectIDs,iSub+zeros(1,length(p)));
        
        if isempty(x.numGoodTrials)
            allNumTrials=cat(2,allNumTrials,zeros(1,length(p)));
        else
            allNumTrials=cat(2,allNumTrials,x.numGoodTrials + zeros(1,length(p)));
        end
    end
end
end