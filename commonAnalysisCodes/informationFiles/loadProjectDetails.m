function [subjectNames,expDates,protocolNames,protocolTypes,capTypes,expDetails,demographicDetails,workbookFolder,FsEyes] = loadProjectDetails(projectName,forcedCapTypeFlag) %#ok<*STOUT>
% dbstop if error
    if ~exist('forcedCapTypeFlag','var') || isempty(forcedCapTypeFlag); forcedCapTypeFlag = 0; end % MD 04-07-2019
    
    clear workbookFile workbookFolder projectStruct 
    workbookFile = which([projectName 'Details.mat']);
    workbookFolder = fileparts(workbookFile);

    load(fullfile(workbookFolder,[projectName 'Details.mat']));
    
    expDetailsCapTypeCol = strcmpi(expDetails(1,:),'capType'); %#ok<NODEF>
    capTypes = cell(length(subjectNames),1);
    
    expDetailsFsEyeCol = strcmpi(expDetails(1,:),'FsEye');
    FsEyes = cell(length(subjectNames),1);
    
    for iSub = 1:length(subjectNames)
        subjectName = subjectNames{iSub};
        expDetailsSubjectCol = strcmpi(expDetails(:,1),subjectName) & ~strcmpi(expDetails(:,5),'No'); % Data usable? condition added by MD 25-03-2019
        if sum(expDetailsSubjectCol)>0             
            capTypes{iSub} = expDetails{expDetailsSubjectCol,expDetailsCapTypeCol};            
            FsEyes{iSub} = expDetails{expDetailsSubjectCol,expDetailsFsEyeCol}; % MD 06-JAN-2019
        end
        if forcedCapTypeFlag % condition added by MD 04-07-2019 for getting info about captype when required even if data not usable
            capTypes{iSub} = expDetails{strcmpi(expDetails(:,1),subjectName),expDetailsCapTypeCol};
        end        
    end
end