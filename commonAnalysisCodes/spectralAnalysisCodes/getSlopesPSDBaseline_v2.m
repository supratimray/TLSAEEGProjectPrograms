% DisplayAllDataSRC is the main program used for plotting and comparing the
% slopes of the power spectra for the SRC protocol. It is modified from displayDataV0 

% This is modified from displayAllDataSRC. We take *all* highRMS and good
% electrodes each day. This works only for baseline data because stimuli
% are centered differently each day.

% v2 is used for Reference Project. Slopes are computed as a function of
% centre frequency. Default baseline duration is 500 ms. Program is written
% such that once intermediate files are available, original data is not
% required.

% modified from displayAllDataBaselineAllElectrodesv2
% MD/07-Feb-2020: added outputs: amps and slopes
% MD/02-July-2020: Comment: Output 'out' added by Santosh; works only when one freq
% band is given as inout (computeSlopesAtFreqs) for slope calculation

function  [out,amps,slopes] = getSlopesPSDBaseline_v2(log10PSD,fAxis,computeSlopesAtFreqs,fitDataLengthHz,freqsToAvoid)

% if size(log10PSD,1)>size(log10PSD,2); log10PSD = log10PSD'; end; % NON_SENSE

if ~exist('computeSlopesAtFreqs','var') || isempty(computeSlopesAtFreqs); computeSlopesAtFreqs = 20:10:100; end;
if ~exist('fitDataLengthHz','var') || isempty(fitDataLengthHz); fitDataLengthHz=30+zeros(1,length(computeSlopesAtFreqs)); end;

for trialNum=1:size(log10PSD,1)
    logPO{trialNum} = log10PSD(trialNum,:); % MD 24-12-2017
%     logPO{trialNum} = conv2Log(rawPSDBL(trialNum,:)); 
%     colorNames(trialNum,:) = getColorRGB(trialNum); % Commented out by MD: 06-05-2019
end
[amps,slopes] = getSlopes(logPO,fAxis,fitDataLengthHz,computeSlopesAtFreqs,freqsToAvoid);
out = [amps,slopes]; % MD: added by Santosh; works only if 
end

function [amps,slopes] = getSlopes(dataToPlot,ys,fitDataLengthHz,computeSlopesAtFreqs,freqsToAvoid)

if iscell(computeSlopesAtFreqs) % MD
    fitRanges = computeSlopesAtFreqs; % MD
else % MD
    for i=1:length(computeSlopesAtFreqs)
    %     fitRanges{i} = computeSlopesAtFreqs(i) + [-fitDataLengthHz/2 fitDataLengthHz/2]; %#ok<*AGROW>
        fitRanges{i} = computeSlopesAtFreqs(i) + [-fitDataLengthHz(i)/2 fitDataLengthHz(i)/2]; %#ok<*AGROW> % MD
    end
end % MD

numDataEntries = length(dataToPlot);
numRanges     = length(fitRanges);

% disp('Fitting curve to data...');
clear signal N
for i=1:numDataEntries
    
    for j=1:numRanges
%         disp(['Entry ' num2str(i) ' of ' num2str(numDataEntries) ', range ' num2str(j) ' of ' num2str(numRanges)]);
        [xData,exitFlag{i,j}]=fitPowerFunction2D(dataToPlot{i},ys,fitRanges{j},freqsToAvoid); %#ok<*NASGU>
        amps{i,j} = xData(:,1); % MD 09-08-2018
        slopes{i,j} = xData(:,2);
%         noiseFloor{i,j} = xData(:,3);
    end
end
end

function plotSlopes(hPlot,slopes,computeSlopesAtFreqs,colorNames,freqLim)

numDataEntries = size(slopes,1);
numRanges = size(slopes,2);
yLim = [0 3];

for i=1:numDataEntries
    for j=1:numRanges
%         mSlopeVals(i,j) = median(xData{i,j}(:,2));
%         seSlopeVals(i,j) = getSEMedian(xData{i,j}(:,2),1000);
        mSlopeVals(i,j) = mean(slopes{i,j});
        seSlopeVals(i,j) = std(slopes{i,j})/sqrt(length(slopes{i,j}));
    end
    plot(hPlot,computeSlopesAtFreqs,mSlopeVals(i,:),'color',colorNames(i,:),'LineWidth',2);
    hold(hPlot,'on');
    plot(hPlot,computeSlopesAtFreqs,mSlopeVals(i,:)+seSlopeVals(i,:),'color',colorNames(i,:),'linestyle','--');
    plot(hPlot,computeSlopesAtFreqs,mSlopeVals(i,:)-seSlopeVals(i,:),'color',colorNames(i,:),'linestyle','--');
end

% axis(hPlot,[freqLim yLim]);
% set(hPlot,'XScale','log');
set(hPlot,'XScale','linear');

for j=1:numRanges
    clear slopeVals slopeValsForAnova groupValsForAnova
    
    slopeValsForAnova = [];
    groupValsForAnova = [];
    
    for i=1:numDataEntries
        slopeVals{i} = slopes{i,j};
        N = length(slopeVals{i});
        slopeValsForAnova = [slopeValsForAnova;slopeVals{i}];
        groupValsForAnova = [groupValsForAnova;i+zeros(N,1)];
    end
    
    % Anova test to compare means
    [pValues(j),anovatab{j},stats{j}]=anova1(slopeValsForAnova,groupValsForAnova,'off');
end

% sigPoints = find(pValues<0.05);
% if ~isempty(sigPoints)
%     showSignificance(hPlot,computeSlopesAtFreqs,sigPoints,'g',yLim(1),0.05);
% end
sigPointsBonferroni = find(pValues<0.05/length(computeSlopesAtFreqs));
if ~isempty(sigPointsBonferroni)
    showSignificance(hPlot,computeSlopesAtFreqs,sigPointsBonferroni,'k',yLim(1),0.05);
end
end

function [x,exitFlag]=fitPowerFunction2D(logP,ys,fitRange,freqsToAvoid)
N=size(logP,1);
for i=1:N
    [x(i,:),exitFlag(i)]=fitPowerFunction(logP(i,:),ys,fitRange,freqsToAvoid);
end
end
function [x,exitFlag]=fitPowerFunction(logP,ys,fitRange,freqListToRemove)

% logP = logP(:)'; % MD
% ys = ys(:)'; % MD

if ~exist('freqListToRemove','var') || isempty(freqListToRemove) % MD
    freqListToRemove{1}=[98 102]; % [97 103]; %MD
    for i=1:8
        freqListToRemove{i+1}=[50*i-2 50*i+2];
    end
end % MD

indicesToRemove = [];

for i=1:length(freqListToRemove)
    indicesToRemove = cat(2,indicesToRemove,intersect(find(ys>=freqListToRemove{i}(1)),find(ys<=freqListToRemove{i}(2))));
end 

fPos = setdiff(intersect(find(ys>=fitRange(1)),find(ys<=fitRange(2))),indicesToRemove);
ys   = ys(fPos);
logP = logP(fPos);

opts = optimset('TolX',1e-6,'TolFun',1e-6,'MaxIter',5000,...
    'Display','off','LargeScale','off','MaxFunEvals',5000);

% Optimize Power
P = 10.^logP;
X0 = 1.5; 
% C0 = P(length(P));
A0 = P(1);%-C0;
startPt = [A0 X0];% C0];

% Optimize Power
% x = fminsearch(@(x) mseFunction1(x,P,ys),startPt,opts);

% Optimize Log Power
[x,~,exitFlag] = fminsearch(@(x) mseFunction2(x,logP,ys),startPt,opts);

%     function mse = mseFunction1(x,yData,xData)
%         yHat = powerFunction(x,xData);
%         mse = sum((yData-yHat).^2);
%     end
    function mse = mseFunction2(x,yData,xData)
        yHat = conv2Log(powerFunction(x,xData));
        mse = sum((yData-yHat).^2);
    end
end
function y=powerFunction(x,xData)
A=x(1);X=x(2);%C=x(3);

y = A*(xData).^(-X);% + max(C,0); % Commented out to remove noise floor WSK: 26-08-2019
% y=A*(xData).^(-X) + C;
end
function showSignificance(hPlot,xVals,sigPoints,colorName,yVal,dY)

if ~exist('dY','var');                   dY = 1;                        end

if ~isempty(sigPoints)
    dX = xVals(2)-xVals(1);
    hold(hPlot,'on');
    for i=1:length(sigPoints)
        patchX = xVals(sigPoints(i))-dX/2;
        patchY = yVal;
        patchLocX = [patchX patchX patchX+dX patchX+dX];
        patchLocY = [patchY patchY+dY patchY+dY patchY];
        patch(patchLocX,patchLocY,colorName,'Parent',hPlot,'EdgeColor',colorName);
    end
end
end