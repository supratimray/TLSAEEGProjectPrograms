function [expDates,protocolNames,capLayout,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType)

% disp(['subjectName: ' subjectName]); % MD 12-09-2019
usableDataFlag = 1;
if strcmpi(projectName,'VisualGamma')
    [subjectNamesAll,expDatesAll,protocolNamesAll,~,~,capLayoutsAll,protocolTypesAll] = allProtocolsVisualGammaEEG;
else
    [subjectNamesAll,expDatesAll,protocolNamesAll,protocolTypesAll,capLayoutsAll,expDetails] = loadProjectDetails(projectName);
    expDetailsSubjectCol = strcmpi(expDetails(:,1),subjectName) & ~strcmpi(expDetails(:,5),'No'); % Data usable? condition added by MD 10-09-2019
    if sum(expDetailsSubjectCol)==0
        disp(['Data not usable for subject ' subjectName]);
        usableDataFlag = 0;
    end
end

numList = length(subjectNamesAll);
pos = [];
for j=1:numList
    if strcmpi(subjectName,subjectNamesAll{j}) && strncmpi(protocolTypesAll{j},protocolType,6)
        pos = cat(2,pos,j);
    end
end
if isempty(pos); disp(['No protocols listed for the subject ' subjectName]); end
expDates = expDatesAll(pos);
protocolNames = protocolNamesAll(pos);
capLayout = capLayoutsAll(pos);
end