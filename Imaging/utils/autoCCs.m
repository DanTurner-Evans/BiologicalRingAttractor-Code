function [tLags, CCVals] = autoCCs(tPts,velVals,actVals,numPts)
% [tLags, CCVals] = autoCCs(tPts,velVals,actVals,numPts)
%   Find the autocorrelation between velocity and calcium activity
%
%   Input:
%     tPts          the time points for the data
%     velVals       the velocity values at those times
%     actVals       the activity values at those times
%     numPts        the number of points to shift for the autocorrelation
%
%   Output:
%     tLags         the time lags for the output correlations
%     CCVals        the correlation values

    CCVals = zeros(numPts,1);
    
    for lag=1:numPts
        lagCCs = corrcoef(...
            velVals(floor(numPts/2):end-ceil(numPts/2)),...
            actVals(numPts-lag+1:end-lag+1));
        CCVals(lag) = lagCCs(2,1);
    end
    
    tLags = mean(diff(tPts))*([1:numPts] - 1 - numPts + floor(numPts/2));

end