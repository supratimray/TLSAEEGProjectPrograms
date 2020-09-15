% 19-05-2017 MD: Added alphaPos, stPower and blPower for alpha
% 11-Sep-2020 SR: Modified to include position list of any length

function diffPower = getDiffPowerSingleSubject(stPowerVsFreq,blPowerVsFreq,posList)
            
    for i=1:length(posList)
        stPower(:,i) = sum(stPowerVsFreq(:,posList{i}),2);  %#ok<*AGROW>
        blPower(:,i) = sum(blPowerVsFreq(:,posList{i}),2);
    end
    
    diffPower = 10*log10(stPower./blPower);
end