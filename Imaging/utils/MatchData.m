function matchedData = MatchData(rawData,positionDat)
% matchedData = MatchData(rawData,positionDat)
%   Matches behavioral position data to the framegrab times from calcium
%   imaging
%
%   Input:
%     rawData       the data to be matched
%     positionDat   a structure holding the experimental data
%
%   Output:
%     matchedData	the matched data

    matchedData = zeros(positionDat.maxFG-positionDat.minFG+1,1);
    % maxFG is the last frame and minFG is the first frame that overlap
    % with the behaviroal recordings
    for i = positionDat.minFG:positionDat.maxFG
        tMatch = find(positionDat.t >=...
            (positionDat.t(1) +...
            (positionDat.tFrameGrab((i-1)*positionDat.num_planes+1)-...
            positionDat.tVR(1))/10000)); 
        % the TTL pulses are recorded via a DAQ at 10000 Hz
        matchedData(i-positionDat.minFG+1) = rawData(tMatch(1));
    end
end
