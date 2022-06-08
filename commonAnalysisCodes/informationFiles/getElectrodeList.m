% This code by default gives a cell array of 10 electrodes when unipolar is
% chosen, and 3 groups each with 3 electrode pairs when bipolar is chosen.
% An option is now added so that for the unipolar case, the 3 separate
% groups are also returned.

function [electrodeList,elecTags,elecNums] = getElectrodeList(capLayout,refType,topoplotFlag,returnSeparateUnipolarGroupsFlag)

if ~exist('topoplotFlag','var');            topoplotFlag=0;             end
if ~exist('returnSeparateUnipolarGroupsFlag','var');  returnSeparateUnipolarGroupsFlag=0;  end

if topoplotFlag
    electrodeList = getElectrodeListForTopoplots(capLayout,refType);
    elecTags = [];
    elecNums = [];
    
else % get selected electrodes
    
    if strcmp(capLayout,'actiCap64')
        if strcmp(refType,'unipolar') || strcmp(refType,'laplacian')
            
            % Left side
            iSide = 1; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = 29;              elecTags{iSide}{iElec} = 'O';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 29+32;           elecTags{iSide}{iElec} = 'PO';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 25+32;           elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 24;              elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            
            % Central electrodes common for both left side and right side
            iElec = iElec+1; electrodeList{iSide}{iElec} = 30;              elecTags{iSide}{iElec} = 'Oz';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 30+32;           elecTags{iSide}{iElec} = 'POz';     elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            
            % Right side
            iElec = iElec+1; electrodeList{iSide}{iElec} = 31;              elecTags{iSide}{iElec} = 'O';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 31+32;           elecTags{iSide}{iElec} = 'PO';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 26+32;           elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 26;              elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            
            if returnSeparateUnipolarGroupsFlag
                % Left side
                iSide = 2; iElec = 0;
                iElec = iElec+1; electrodeList{iSide}{iElec} = 29+32;           elecTags{iSide}{iElec} = 'PO';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                iElec = iElec+1; electrodeList{iSide}{iElec} = 25+32;           elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                iElec = iElec+1; electrodeList{iSide}{iElec} = 24;              elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                
                % Right side
                iSide = 3; iElec = 0;
                iElec = iElec+1; electrodeList{iSide}{iElec} = 31+32;           elecTags{iSide}{iElec} = 'PO';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                iElec = iElec+1; electrodeList{iSide}{iElec} = 26+32;           elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                iElec = iElec+1; electrodeList{iSide}{iElec} = 26;              elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                
                % Central electrodes common for both left side and right side
                iSide = 4; iElec = 0;
                iElec = iElec+1; electrodeList{iSide}{iElec} = 30;              elecTags{iSide}{iElec} = 'Oz';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                iElec = iElec+1; electrodeList{iSide}{iElec} = 30+32;           elecTags{iSide}{iElec} = 'POz';     elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                iElec = iElec+1; electrodeList{iSide}{iElec} = 29;              elecTags{iSide}{iElec} = 'O';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
                iElec = iElec+1; electrodeList{iSide}{iElec} = 31;              elecTags{iSide}{iElec} = 'O';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            end
            
        elseif strcmp(refType,'bipolar')
            
            % Left side
            iSide = 1; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [30+32 29+32];    elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 101;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [29+32 25+32];    elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 94;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [29+32 24];       elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 93;
            
            % Right side
            iSide = 2; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [30+32 31+32];    elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 102;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [31+32 26+32];    elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 96;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [31+32 26];       elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 97;
            
            % Centre posterior
            iSide = 3; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [30 29];          elecTags{iSide}{iElec} = 'O';      elecNums{iSide}{iElec} = 111;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [30 30+32];       elecTags{iSide}{iElec} = 'MO';     elecNums{iSide}{iElec} = 107;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [30 31];          elecTags{iSide}{iElec} = 'O';      elecNums{iSide}{iElec} = 112;
            
        end
        
        
        
    elseif strcmp(capLayout,'actiCap31Posterior')
        if strcmp(refType,'unipolar')
            
            % Left side
            iSide = 1; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = 28;              elecTags{iSide}{iElec} = 'O';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 23;           elecTags{iSide}{iElec} = 'PO';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 16;           elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 15;              elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            
            % Central electrodes common for both left side and right side
            iElec = iElec+1; electrodeList{iSide}{iElec} = 29;              elecTags{iSide}{iElec} = 'Oz';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 24;           elecTags{iSide}{iElec} = 'POz';     elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            %             iElec = iElec+1; electrodeList{iSide}{iElec} = 21+32;    elecTags{iSide}{iElec} = 'CPz';
            
            % Right side
            %             iSide = 2; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = 30;              elecTags{iSide}{iElec} = 'O';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 25;           elecTags{iSide}{iElec} = 'PO';      elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 18;           elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            iElec = iElec+1; electrodeList{iSide}{iElec} = 19;              elecTags{iSide}{iElec} = 'P';       elecNums{iSide}{iElec} = electrodeList{iSide}{iElec};
            
        elseif strcmp(refType,'bipolar')
            
            % Left side
            iSide = 1; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [24 23];    elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 39;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [23 16];    elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 32;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [23 15];       elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 31;
            
            % Right side
            iSide = 2; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [24 25];    elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 40;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [25 18];    elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 34;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [25 19];       elecTags{iSide}{iElec} = 'PO';     elecNums{iSide}{iElec} = 35;
            
            % Centre posterior
            iSide = 3; iElec = 0;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [29 28];          elecTags{iSide}{iElec} = 'O';      elecNums{iSide}{iElec} = 49;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [29 24];       elecTags{iSide}{iElec} = 'MO';     elecNums{iSide}{iElec} = 45;
            iElec = iElec+1; electrodeList{iSide}{iElec} = [29 30];          elecTags{iSide}{iElec} = 'O';      elecNums{iSide}{iElec} = 50;
            
        end
        
    else
        error('Cap invalid or not specified');
    end
end
end

function electrodeList = getElectrodeListForTopoplots(capLayout,refType)

% Get the electrode list
clear cL bL chanlocs iElec electrodeList
switch refType
    case 'unipolar'
        cL = load([capLayout '.mat']);
        chanlocs = cL.chanlocs;
        for iElec = 1:length(chanlocs)
            electrodeList{iElec}{1} = iElec; %#ok<*AGROW>
        end
    case 'bipolar'
        cL = load(['bipolarChanlocs' capLayout '.mat']);
        bL = load(['bipChInfo' capLayout '.mat']);
        chanlocs = cL.eloc;
        for iElec = 1:length(chanlocs)
            electrodeList{iElec}{1} = bL.bipolarLocs(iElec,:);
        end
end
end