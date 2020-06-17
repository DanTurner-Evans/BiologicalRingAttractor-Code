function [tLags, CCVals] = autoCCs(tPts,velVals,actVals,numPts)

    CCVals = zeros(numPts,1);
    
    for lag=1:numPts
        lagCCs = corrcoef(...
            velVals(floor(numPts/2):end-ceil(numPts/2)),...
            actVals(numPts-lag+1:end-lag+1));
        CCVals(lag) = lagCCs(2,1);
    end
    
    tLags = mean(diff(tPts))*([1:numPts] - 1 - numPts + floor(numPts/2));

end