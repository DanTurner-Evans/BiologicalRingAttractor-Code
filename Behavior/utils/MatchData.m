function matchedData = MatchData(rawData,positionDat)

    matchedData = zeros(positionDat.maxFG-positionDat.minFG+1,1);
    for i = positionDat.minFG:positionDat.maxFG
        tMatch = find(positionDat.t >=...
            (positionDat.t(1) +...
            (positionDat.tFrameGrab((i-1)*positionDat.num_planes+1)-...
            positionDat.tVR(1))/10000)); 
        matchedData(i-positionDat.minFG+1) = rawData(tMatch(1));
    end
end
