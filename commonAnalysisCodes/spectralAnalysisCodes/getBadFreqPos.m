% Remove line noise and it's harmonics.
% For India, line noise is at 50 Hz.

function badFreqPos = getBadFreqPos(freqVals,deltaF)
badFreqs = 50:50:max(freqVals);
if nargin<2; deltaF = 1; end;

badFreqPos = [];
for i=1:length(badFreqs)
    badFreqPos = cat(2,badFreqPos,intersect(find(freqVals>=badFreqs(i)-deltaF),find(freqVals<=badFreqs(i)+deltaF)));
end
end