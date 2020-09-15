% We take subjectNames from demographicDetails. However, for some subjects
% there is descrepancy between names in demographic details and fileNames
% that are used for saving data. This program corrects those names.

function subjectFileNames = getGoodFileNamesForSubjects(subjectNames)

N = length(subjectNames);
subjectFileNames = cell(1,N);

for i=1:N
    
    s = subjectNames{i};
    
    if strcmp(s,'011KH_F1'); subjectFileNames{i} = '011KH';
    elseif strcmp(s,'144BR'); subjectFileNames{i} = '144RP';
    elseif strcmp(s,'227SN'); subjectFileNames{i} = '227MN';
    elseif strcmp(s,'228MN'); subjectFileNames{i} = '228SN';
    elseif strcmp(s,'236SD'); subjectFileNames{i} = '236S';
    elseif strcmp(s,'295R'); subjectFileNames{i} = '295RS';
    else
        subjectFileNames{i} = s;
    end
end
end