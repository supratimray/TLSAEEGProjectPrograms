% projectName - ADGammaProject, AgeProjectRound1 or VisualGamma. Only ADGammaProject option is tested for now.

function [ageList,genderList,cdrList,expDateList] = getDemographicDetails(projectName,subjectNameList)

if ~exist('projectName','var');     projectName = 'ADGammaProject';     end

if strcmpi(projectName,'VisualGamma')
    % This segment not written yet
else
    subjectNameList = getGoodFileNamesForSubjects(subjectNameList);
    
    d = load([projectName 'Details.mat']);
    demographicDetails = d.demographicDetails;
    N = size(demographicDetails,1);
    allSessionNames = getGoodFileNamesForSubjects(demographicDetails(2:N,1));
    allAge = demographicDetails(2:N,2);
    allGender = demographicDetails(2:N,3);
    allCDR = demographicDetails(2:N,5);
    
    expDetails = d.expDetails;
    Nexp = length(expDetails);
    allSessionNamesExp = getGoodFileNamesForSubjects(expDetails(2:Nexp,1));
    allExpDates = expDetails(2:Nexp,2);
    
    numSubjects = length(subjectNameList);
    ageList = zeros(1,numSubjects);
    genderList = cell(1,numSubjects);
    cdrList = cell(1,numSubjects);
    expDateList = cell(1,numSubjects);
    
    for i=1:numSubjects
        pos = strcmp(subjectNameList{i},allSessionNames);
        a = allAge(pos);
        ageList(i) = str2double(a{1}(1:end-1));
        genderList{i} = allGender{pos};
        cdrList{i} = allCDR{pos};
        
        posExp = strcmp(subjectNameList{i},allSessionNamesExp);
        expDateList{i} = allExpDates{posExp};
    end
end
end