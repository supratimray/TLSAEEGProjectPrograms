function subjectNameListMatched2 = getNewControlSubjectList(projectName,subjectNameListMatched,maxNumberOfControls)

clear subjectNameListMatched2
subjectNameListMatched2{2}=subjectNameListMatched{2}; % CaseList
numValidCases = length(subjectNameListMatched{2});

for i=1:numValidCases
    [~,~,~,expDate0] = getDemographicDetails(projectName,subjectNameListMatched{2}{i});
    allControls = subjectNameListMatched{1}{i};
    [~,~,~,expDateList] = getDemographicDetails(projectName,allControls);
    
    x = expDate0{1}; xStr = [x(3:4) '/' x(1:2) '/' x(5:6)];
    dateNumX = datenum(xStr,'mm/dd/yy');
    
    numEntries = length(expDateList);
    numDays = zeros(1,numEntries);
    for j=1:numEntries
        y = expDateList{j}; yStr = [y(3:4) '/' y(1:2) '/' y(5:6)];
        numDays(j) = datenum(yStr,'mm/dd/yy') - dateNumX;
    end
    
    [~,sortOrder] = sort(abs(numDays));
    subjectNameListMatched2{1}{i} = allControls(sortOrder(1:min(numEntries,maxNumberOfControls)));
end
