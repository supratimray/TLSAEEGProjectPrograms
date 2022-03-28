% Written by Wupadrasta Santosh
% Modified by Supratim Ray

function analyseAndSaveValuesIndividualSubjectConn(folderSourceString,subjectName,projectName,subProjectName,refType,protocolType,connMethod,stRange,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData)

if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('spatialFrequenciesToRemove','var'); spatialFrequenciesToRemove=[];  end
if ~exist('useCleanData','var');    useCleanData=0;                     end

analyzedDataFolder = fullfile(pwd,'analyzedData',subProjectName,protocolType); % analysedFolder now in local project directory
makeDirectory(analyzedDataFolder);

% Analysis file - this is where the final analyzed data is saved
[analysisDetailsFile,~] = getAnalysisDetailsFile(analyzedDataFolder,subjectName,refType,protocolType,stRange,removeMicroSaccadesFlag,spatialFrequenciesToRemove,useCleanData,[]);
analysisDetailsFileConn = [analysisDetailsFile(1:end-4) '_' connMethod '.mat']; % for saving connectivity
analysisDetailsFileFreq = [analysisDetailsFile(1:end-4) '_Freq.mat']; % for saving PSDs etc

% ftData - data in fieldtrip format is saved in ftData
ftDataFolder = fullfile(folderSourceString,'ftData',projectName,protocolType);

[expDates,protocolNames,capType,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);

if usableDataFlag && ~isempty(expDates) && strcmp(capType{1},'actiCap64')
    
    [dataPost,dataPre,goodProtFlag,numGoodTrials] = getDataSingleSubjectFT(subjectName,expDates,protocolNames,refType,ftDataFolder,stRange);
    
    if ~exist(analysisDetailsFileFreq,'file')
        [freqPost,freqPre] = getFreqData(dataPost,dataPre,goodProtFlag);
        save(analysisDetailsFileFreq,'freqPost','freqPre','goodProtFlag','numGoodTrials');
    end
    if ~exist(analysisDetailsFileConn,'file')
        conn = getConnIndividualSubject(dataPost,connMethod,goodProtFlag);
        save(analysisDetailsFileConn,'conn','goodProtFlag','numGoodTrials');
    end
end
end
function [dataPost,dataPre,goodProtFlag,numGoodTrials]=getDataSingleSubjectFT(subjectName,expDates,protocolNames,refType,ftDataFolder,stRange)

numSessions = length(expDates);
dataPre = cell(1,numSessions);
dataPost = cell(1,numSessions);

goodProtFlag = logical([]);
numGoodTrials = zeros(1,numSessions);
for i=1:numSessions
    % get Data
    dataFileName = fullfile(ftDataFolder,[subjectName '-' expDates{i} '-' protocolNames{i} '.mat']);
    x = load(dataFileName);
    data = x.data;
    goodProtFlag = cat(2,goodProtFlag,x.goodProtFlag);
    numGoodTrials(i) = x.numGoodTrials;
    
    % Compute Laplacian if needed
    if(strncmpi(refType,'Lap',3))      
        cfg = [];
        cfg.elec = data.elec;
        data = ft_scalpcurrentdensity(cfg,data);
    end
    
    cfg        = [];
    cfg.toilim = [-diff(stRange) 0]; % just to have data of same length in Pre, Post
    dataPre{i}   = ft_redefinetrial(cfg,data);
    
    cfg        = [];
    cfg.toilim = stRange;
    dataPost{i}   = ft_redefinetrial(cfg,data);
end
end
function [freqPost,freqPre]=getFreqData(dataPost,dataPre,goodProtFlag)

numSessions = length(dataPost);
freqPreList = cell(1,numSessions);
freqPostList = cell(1,numSessions);

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 2;  % +-1Hz smoothing
cfg.foilim    = [0 100];
cfg.keeptrials  = 'no';
    
for i=1:numSessions
    freqPostList{i} = ft_freqanalysis(cfg,dataPost{i});
    freqPreList{i} = ft_freqanalysis(cfg,dataPre{i});
end

freqPreList = freqPreList(goodProtFlag);
freqPostList = freqPostList(goodProtFlag);

cfg = [];
cfg.parameter = 'powspctrm';
freqPost = eval(getCommandStr(length(freqPostList),'Post'));
freqPre = eval(getCommandStr(length(freqPreList),'Pre'));

end

function commandStr = getCommandStr(numSessions,tag)

commandStr = 'ft_freqgrandaverage(cfg,';
for i=1:numSessions
    commandStr = cat(2,commandStr,['freq' tag 'List{' num2str(i) '},']);
end
commandStr = [commandStr(1:end-1) ')'];
end
function connVsElec = getConnIndividualSubject(data,connMethod,goodProtFlag,excludeHighConnElecsFlag)

if ~exist('excludeHighConnElecsFlag','var');    excludeHighConnElecsFlag=0; end

badElecs_accum = [];
for prot = 1:length(data) % across protocols
    temp1=getConnThisCondition(data{prot},connMethod);
    connList(prot,:,:,:) = abs(temp1.([connMethod 'spctrm'])); %#ok<AGROW>
    badElecs_accum = union(badElecs_accum,data{prot}.badElecs);
end
conn = removeDimIfSingleton(nanmean(connList(goodProtFlag,:,:,:),1)); % averaging across protocols

if(excludeHighConnElecsFlag)
    connBadElecs = findBadElecsfromConn(conn);
    totalBadElecs = union(badElecs_accum,connBadElecs);
else
    totalBadElecs = badElecs_accum;
end

conn(totalBadElecs,:,:) = NaN;
conn(:,totalBadElecs,:) = NaN;
connVsElec = conn;
end
function conn_stat = getConnThisCondition(data,method)
cfg              = [];
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq = 2;
if(strcmp(method,'coh'))
    cfg.output = 'powandcsd';
else
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
end
cfg.foilim          = [0 100];
freqPost = ft_freqanalysis(cfg, data);

if(strcmp(method,'coh'))
    TfreqPost = ft_checkdata(freqPost,'cmbrepresentation', 'full');
else
    TfreqPost = freqPost;
end

cfg=[];
cfg.method = method;
conn_stat = ft_connectivityanalysis(cfg, TfreqPost);
end
function connBadElecs = findBadElecsfromConn(conn)
% Very high peaks (as found in wPLI)
freqRanges = {[8 12], [20 34], [36 66]};
f = 0:2:100;
for fr = 1:length(freqRanges)
    freqBins = f>=freqRanges{fr}(1) & f<=freqRanges{fr}(2);
    connBadElecs = [];
    out_band=mean(conn(:,:,freqBins),3);
    for ele = 1:64
        cval = out_band(ele,:);
        cval(ele)=NaN;
        nanindx = isnan(cval);
        connBadElecs = cat(2,connBadElecs,find(cval(~nanindx)>0.95));
    end
end
end