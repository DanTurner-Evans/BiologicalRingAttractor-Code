function [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(data)

    % Find the dark period
    darkPer = find(data.Direction == 0);

    % Find the open loop period
    OLPer = find(data.OLGain > 0);

    % Find the clockwise rotations
    CWPer = find(data.Direction > 0);
    CWPer = intersect(OLPer,CWPer);
    % Find the counter-clockwise rotations
    CCWPer = find(data.Direction < 0);
    CCWPer = intersect(OLPer,CCWPer);

    % Find the closed loop period
    CLPer = find(data.Direction > 0);
    CLPer = setdiff(CLPer,CWPer);

end