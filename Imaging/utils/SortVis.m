function [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(data)
% [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(data)
%   Find the different time periods in a trial where a given visual
%   condition is shown
%
%   Input:
%     data       a structure holding the experimental data
%
%   Output:
%     darkPer    the periods where the fly is in the dark
%     OLPer      the periods where the fly is in open-loop
%     CLPer      the periods where the fly is in closed-loop
%     CWPer      the periods where the fly is in open-loop with CW movement
%     CCWPer     the periods where the fly is in open-loop with CCW movement

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