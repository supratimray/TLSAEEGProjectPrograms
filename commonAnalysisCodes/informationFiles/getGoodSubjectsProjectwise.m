% Get good subjects for each subProject
% ageProject - only healthy subjects (N=236; 9 discarded, finally - 227)
% ADGammaProject - healthy/MCI/AD (N=236/15/6; 9/1/0 discarded, finally - 227/14/6)
% consistencyProject - N=48 healthy subjects and their repeats

function goodSubjectsList = getGoodSubjectsProjectwise(subProjectName,discardNoUsefulSessionsFlag,protocolType)

if ~exist('discardNoUsefulSessionsFlag','var'); discardNoUsefulSessionsFlag=1; end
if ~exist('protocolType','var'); protocolType = 'SF_ORI';               end

projectName = 'ADGammaProject';
[goodSubjects,goodSubjectsY0,goodSubjectsY1] = getGoodSubjects(projectName,discardNoUsefulSessionsFlag,protocolType);

if (~exist('subProjectName','var') ||  strncmpi(subProjectName,'ADGamma',7))
    % Do nothing. Return goodSubjects
    goodSubjectsList{1} = goodSubjects;

elseif strncmpi(subProjectName,'age',3)     % AgeProject - only healthy baseline
    
    d = load([projectName 'Details.mat']);
    demographicDetails = d.demographicDetails;
    sessionNames = demographicDetails(2:size(demographicDetails,1),1);
    labels = demographicDetails(2:size(demographicDetails,1),5);
    
    numSubjects = length(goodSubjects);
    goodLabels = cell(1,numSubjects);
    for i=1:numSubjects
        goodLabels{i} = labels{strcmp(goodSubjects(i),sessionNames)};
    end
    
    goodSubjects = goodSubjects(strcmp(goodLabels,'HV'));
    goodSubjectsList{1} = goodSubjects;

elseif strncmpi(subProjectName,'con',3)
    
    d = load([projectName 'Details.mat']);
    demographicDetails = d.demographicDetails;
    sessionNames = demographicDetails(2:size(demographicDetails,1),1);
    labels = demographicDetails(2:size(demographicDetails,1),5);
    
    numGoodRepeats = length(goodSubjectsY0);
    goodLabelsY0 = cell(1,numGoodRepeats);
    goodLabelsY1 = cell(1,numGoodRepeats);
    for i=1:numGoodRepeats
        goodLabelsY0{i} = labels{strcmp(goodSubjectsY0(i),sessionNames)};
        goodLabelsY1{i} = labels{strcmp(goodSubjectsY1(i),sessionNames)};
    end
    
    goodPos = find(strcmp(goodLabelsY0,'HV') & strcmp(goodLabelsY1,'HV'));
    goodSubjectsHVY0 = goodSubjectsY0(goodPos);
    goodSubjectsHVY1 = goodSubjectsY1(goodPos);
    
    clear goodSubjects
    goodSubjectsList{1} = goodSubjectsHVY0;
    goodSubjectsList{2} = goodSubjectsHVY1;
end