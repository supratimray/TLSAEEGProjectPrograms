% projectName - ADGammaProject, AgeProjectRound1 or VisualGamma. Only ADGammaProject option is tested for now.

function [goodSubjects,goodSubjectsY0,goodSubjectsY1] = getGoodSubjects(projectName,discardNoUsefulSessionsFlag,protocolType)

if ~exist('projectName','var');     projectName = 'ADGammaProject';     end
if ~exist('discardNoUsefulSessionsFlag','var'); discardNoUsefulSessionsFlag=1; end
if ~exist('protocolType','var'); protocolType = 'SF_ORI';               end

if strcmpi(projectName,'VisualGamma')
    % This segment not written yet
else
    d = load([projectName 'Details.mat']);
    demographicDetails = d.demographicDetails;
    
    sessionNames = demographicDetails(2:size(demographicDetails,1),1);
    disp(['Total number of sessions: ' num2str(length(sessionNames))]);
    
    rejectedSessions = getRejectedSessionList;
    if discardNoUsefulSessionsFlag
        noUsefulSessions = getNoUsefulSessionList;
        rejectedSessions = cat(1,rejectedSessions,noUsefulSessions);
    end
    sessionNames = setdiff(sessionNames,rejectedSessions);
    numGoodSessions = length(sessionNames);
    disp(['Total sessions discarded: ' num2str(length(rejectedSessions))]);
    disp(['Remaining good sessions: ' num2str(numGoodSessions)]);
        
    % Find subjectNumber for each session
    subjectNumbers = zeros(1,numGoodSessions);
    
    for i=1:numGoodSessions
        subjectNumbers(i) = str2double(sessionNames{i}(1:3));
    end
        
    % Get good subjects
    uniqueSubjectNumbers = unique(subjectNumbers);
    disp(['Number of unique subjects: ' num2str(length(uniqueSubjectNumbers))]);
    
    goodSubjects = []; % This is for Age and ADGammaProjects
    goodSubjectsY0 = []; % This is for consistency project
    goodSubjectsY1 = []; % This is for consistency project

    for i=1:length(uniqueSubjectNumbers)
        pos = find(uniqueSubjectNumbers(i)==subjectNumbers);
        if (length(pos)==1) % Only one subject
            goodSubjects = cat(2,goodSubjects,sessionNames(pos));
        else % could be 2 or more, but in this dataset we have exactly 2 for each repeat case
            if isnan(str2double(sessionNames{pos(1)}(end))) % first one does not end in a number and hence is baseline
                goodSubjects = cat(2,goodSubjects,sessionNames(pos(1)));
                goodSubjectsY0 = cat(2,goodSubjectsY0,sessionNames(pos(1)));
                goodSubjectsY1 = cat(2,goodSubjectsY1,sessionNames(pos(2)));
            elseif isnan(str2double(sessionNames{pos(2)}(end))) % second one does not end in a number and hence is baseline
                goodSubjects = cat(2,goodSubjects,sessionNames(pos(2)));
                goodSubjectsY0 = cat(2,goodSubjectsY0,sessionNames(pos(2)));
                goodSubjectsY1 = cat(2,goodSubjectsY1,sessionNames(pos(1)));
            elseif (str2double(sessionNames{pos(1)}(end)) == 1) % First entry is F1. Second cannot be baseline because that is checked already. Hence must be F2 or F3
                goodSubjects = cat(2,goodSubjects,sessionNames(pos(1)));
                goodSubjectsY0 = cat(2,goodSubjectsY0,sessionNames(pos(1)));
                goodSubjectsY1 = cat(2,goodSubjectsY1,sessionNames(pos(2)));
            elseif (str2double(sessionNames{pos(2)}(end)) == 1) % Second entry is F1. First cannot be baseline because that is checked already. Hence must be F2 or F3
                goodSubjects = cat(2,goodSubjects,sessionNames(pos(2)));
                goodSubjectsY0 = cat(2,goodSubjectsY0,sessionNames(pos(2)));
                goodSubjectsY1 = cat(2,goodSubjectsY1,sessionNames(pos(1)));
            elseif (str2double(sessionNames{pos(1)}(end)) == 2) % First entry is F2. Second cannot be baseline or F1 because that is checked already. Hence must be F3
                goodSubjects = cat(2,goodSubjects,sessionNames(pos(1)));
                goodSubjectsY0 = cat(2,goodSubjectsY0,sessionNames(pos(1)));
                goodSubjectsY1 = cat(2,goodSubjectsY1,sessionNames(pos(2)));
            elseif (str2double(sessionNames{pos(2)}(end)) == 2) % Second entry is F2. First cannot be baseline or F1 because that is checked already. Hence must be F3
                goodSubjects = cat(2,goodSubjects,sessionNames(pos(2)));
                goodSubjectsY0 = cat(2,goodSubjectsY0,sessionNames(pos(2)));
                goodSubjectsY1 = cat(2,goodSubjectsY1,sessionNames(pos(1)));
            end
        end
    end
end

if strcmp(protocolType,'TFCP')
    badSSVEPSubjects = dataNotUsableTFCP;
    disp([num2str(length(badSSVEPSubjects)) ' subjects discarded for SSVEP protocol']);
    goodSubjects = setdiff(goodSubjects,badSSVEPSubjects);
end
end

function rejectedSessions = getRejectedSessionList

rejectedSessions = {
    % Did not complete experiment
    '013JM_F1';
    '019ST';
    '092PS';
    '079AV_F1';
    
    % Removed from analysis - due to various reasons, such as poor vision, errors in recording, subject falling sleep. 
    % Exact details can be found in the excel sheet database maintained by the experimenters.  
    '068TP_F1';
    '069PS_F1';
    '121MM_F1';
    '145SR_F1';
    '165AR_F1';
    '201AS_F1';
    '228MN'; % labeled as 228SN in subjectNamesAll
    '261GB';
    '267TK_F1';
    '268MK';
    '284PP';
    '284PP_F1';
    '286VP';
    '427AN';
    '449KM';
    '466AR';
    '512RS';
    
    % These are also considered bad because of the reasons listed below, but these subjects have done at least one other session which is good, so we get at least one usable dataset from these subjects
    '002_N'; % No eye data
    '049VK'; % No eye data
    '050MG'; % No eye data
    '083MP'; % No eye data
    '133RC'; % No useful protocols
    '154AP'; % label switched between baseline and F1
    
    % Bad Label - either the diagnosis is pending or there is some discrepancy
    '078MB_F1';
    '100JK_F1';
    '104GS_F1';
    '107GR_F1'; % Bad label + Repeat
    '183KP_F1'; % Bad label + Repeat
    '255IS';
    '342RS';
    '368K';
    '390GK';
     
    % Age less than 50
    '509JY';
    
    % These three sessions were repeats, but their original session was later
    % removed from analysis. These should have therefore been used as
    % baseline data for these subjects, but the labels were not changed from "repeat" to "useful" by mistake.
    '261GB_F1';
    '286VP_F1';
    '205MG_F1';
    };
end
function noUsefulSessions = getNoUsefulSessionList

% These subjects were analyzed, but upon removal of good electrodes, no
% useful protocols were left. These can also discarded from further analysis
noUsefulSessions = {
    '053JS_F2';
    '205MG';
    '242VA';
    '323NV';
    '327GH';
    '328MR';
    '343PJ';
    '350SP';
    '402MS';
    '405DM';
    '282AB_F1'; % Repeat. Original can still be used.
    };
end
function badSSVEPSubjects = dataNotUsableTFCP
badSSVEPSubjects{1} = '014SM';
badSSVEPSubjects{2} = '084VB';
end