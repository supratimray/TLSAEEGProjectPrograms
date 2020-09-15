function [labelList,genderList] = describeSessionDetails(sList,projectName)

if ~exist('projectName','var');     projectName = 'ADGammaProject';     end

d = load([projectName 'Details.mat']);
demographicDetails = d.demographicDetails;
sessionListAll = demographicDetails(2:size(demographicDetails,1),1);
genderListAll = demographicDetails(2:size(demographicDetails,1),3);
labelListAll = demographicDetails(2:size(demographicDetails,1),5);

% Convert both lists to a consistent format
sList = getGoodFileNamesForSubjects(sList);
sessionListAll = getGoodFileNamesForSubjects(sessionListAll);

% Find subjectNumber for each session
numSessions = length(sList);

% characterize sessions as HV/MCI/AD and M/F
genderList = cell(1,numSessions);
labelList = cell(1,numSessions);

for i=1:numSessions
    sID = find(strcmp(sessionListAll,sList{i}));
    genderList{i} = genderListAll{sID};
    labelList{i} = labelListAll{sID};
end

disp(['NumSessions=' num2str(numSessions) ', HV/MCI/AD=' num2str(sum(strcmp(labelList,'HV'))) ...
    '/' num2str(sum(strcmp(labelList,'MCI'))) '/' num2str(sum(strcmp(labelList,'AD')))...
    ', F=' num2str(sum(strcmp(genderList,'F')))]);

end