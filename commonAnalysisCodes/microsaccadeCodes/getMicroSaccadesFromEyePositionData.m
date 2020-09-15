%
% getMicroSaccadesFromEyePositionData
% Detect microsaccades from eye data
%
% Inputs:
%       eyeDataDegX & eyeDataDegY: eye data along X and Y axes (Deg)
%       xs: timeVals (s)
%       FsEye: Sampling frequency of eye data (Hz)
%       timeRange: time range of eye data for microsaccade analysis (s) 
%       threshold: no. of times the SD of eye data for determining cutoff
%           of eye velocity
%       minCutOff: minimum cutoff of velocity 
%       minMSLength: minimum duration of eye movement to be designated as
%           microsaccade
%
% Outputs:
%       MS: starting times of microsaccades
%       MSLength: Duration of each microsaccade (s)
%       numMSInRange: number of microsaccades for each trial in analysis
%           period
%       eyeSpeedX & eyeSpeedY: eye data along X and Y axes
%       cutoffX and cutoffY: cutoff of eye-velocities used for microsaccade detection
%       eyeSpeedMag = magnitude of eye-velocities
%
% Murty V P S Dinavahi 15-09-2017
% Based on description of microsaccade detection given in:
% Engbert, R., 2006. Microsaccades: a microcosm for research on oculomotor control, attention, and visual perception. 
%       Progress in Brain Research, Visual Perception 154, 177–192. doi:10.1016/S0079-6123(06)54009-9
%
% This code is meant only for analysing monocular eye data
%

function [MS,MSLength,numMSInRange,eyeSpeedX,eyeSpeedY,cutoffX,cutoffY,eyeSpeedMag] = getMicroSaccadesFromEyePositionData(eyeDataDegX,eyeDataDegY,xs,FsEye,timeRange,threshold,minCutOff,minMSLength)


if ~exist('timeRange','var');     timeRange=[xs(1) xs(end)];               end
if ~exist('minCutOff','var');     minCutOff=10;               end
if ~exist('minMSLength','var');     minMSLength=0.010;               end % in seconds
commonCutOffFlag = 0;

% Calculate eye-speed
xsRange = xs>=timeRange(1) & xs<=timeRange(2);
eyeSpeedX = NaN(size(eyeDataDegX));
eyeSpeedY = NaN(size(eyeDataDegY));
cutoffX = NaN(1,size(eyeDataDegX,1));
cutoffY = NaN(1,size(eyeDataDegY,1));

for j=1:size(eyeDataDegX,1)
    eyeSpeedX(j,:) = getSpeed(eyeDataDegX(j,:),FsEye); 
    eyeSpeedY(j,:) = getSpeed(eyeDataDegY(j,:),FsEye);
    
    if ~commonCutOffFlag
        cutoffX(j) = max([minCutOff threshold*getCutOff(eyeSpeedX(j,xsRange))]);
        cutoffY(j) = max([minCutOff threshold*getCutOff(eyeSpeedY(j,xsRange))]);
    end
end

clear eyeSpeedMag MStmp numMStmp
eyeSpeedMag = sqrt(eyeSpeedX.^2+eyeSpeedY.^2);

% Calculate cutoff
if commonCutOffFlag
    clear cutoffX cutoffY
    eyeDataForCutOFfX = eyeSpeedX(:,xsRange);
    cutoffX = threshold*getCutOff(eyeDataForCutOFfX(:));
    cutoffX = max([minCutOff cutoffX]);
    cutoffX = repmat(cutoffX,1,size(eyeDataDegX,1));

    eyeDataForCutOFfY = eyeSpeedY(:,xsRange);
    cutoffY = threshold*getCutOff(eyeDataForCutOFfY(:));
    cutoffY = max([minCutOff cutoffY]);
    cutoffY = repmat(cutoffY,1,size(eyeDataDegY,1));
end

% Detect microsaccades
N = size(eyeSpeedMag,1);
MS = cell(1,N);
MSLength = cell(1,N);
numMSInRange = zeros(1,N);
    
for i=1:N
    yX = eyeSpeedX(i,:);    
    msX =  sort([find(yX>cutoffX(i)) find(yX<-cutoffX(i))]);  
    
    yY = eyeSpeedY(i,:);    
    msY =  sort([find(yY>cutoffY(i)) find(yY<-cutoffY(i))]);  
    
    % Check
    msAll = union(msX,msY);
    msAllSpeedX = yX(msAll);
    msAllSpeedY = yY(msAll);
    
    msS = (((msAllSpeedX./cutoffX(i)).^2)+((msAllSpeedY./cutoffY(i)).^2))>1;
    ms = msAll(msS);
    
    if ~isempty(ms)
                
        [msI,msLen] = findMSPositionsAndLengthsPerTrial(ms,minMSLength,FsEye);
        
        MS{i} = xs(ms(msI));
        MSLength{i} = msLen;
    else
        MS{i}=[];
        MSLength{i} = [];
    end
    
    numMSInRange(i) = length(intersect(find(MS{i}>=timeRange(1)),find(MS{i}<=timeRange(2))));
end
end

function v=getSpeed(x,FsEye)
    deltaT = 1/FsEye;
    
    v=NaN(1,length(x));
    for n=3:length(x)-2
        v(n) = (x(n+2)+x(n+1)-x(n-1)-x(n-2))./(6*deltaT);
    end
end

function c=getCutOff(x)
    medX = nanmedian(x);
    medX2 = nanmedian(x.^2);    
    c=sqrt(medX2-medX.^2);    
end

function [msIndex,msLengthInSeconds] = findMSPositionsAndLengthsPerTrial(ms,minMSLength,FsEye)

    clear dms dms2 dmsIn msIn
    dms = [0 diff(ms)];
    minMSSamples = ceil(minMSLength*FsEye)-1;
    if length(ms)<=minMSSamples; msIndex=[]; msLengthInSeconds=[]; return; end;
    dms2 = zeros(1,length(ms)-minMSSamples);
    for iD = 1:length(ms)-minMSSamples
        dms2(iD) = ms(iD)-ms(iD+minMSSamples);
    end
    
    msIndex=find(dms2==-minMSSamples);

    msIn = ms(msIndex);
    dmsIn = [0 diff(msIn)];
    msIndex(dmsIn==1) = [];

    il = 1;
    msLen = ones(1,length(msIndex));
    for iL = msIndex
        for id = iL+1:length(dms)
            if dms(id)>1; 
                break; 
            end
            msLen(il) = msLen(il)+1;
        end
        il = il + 1;
    end
    msLengthInSeconds = msLen./FsEye;
end

