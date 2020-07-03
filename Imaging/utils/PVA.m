function [angs, PVA, PVAStren] = PVA(act)
% [angs, PVA, PVAStren] = PVA(act)
%   Calculate the angles, population vector average (PVA), and PVA strength
%   from an activity matrix
%
%   Input:
%     act       the calcium activity matrix (# of ROIS x # of time points)
%
%   Output:
%     ans       the ROI angles
%     PVA       the PVA at each time point
%     PVAStren  the PVA strength at each time point

    num_ROIs = size(act,1);
    
    angs = (1:num_ROIs)*2*pi/num_ROIs-pi;
    angs = angs';
    
    PVA = zeros(size(act,2),1);
    PVAStren = zeros(size(act,2),1);
    for ts = 1:size(act,2)
        PVA(ts) = circ_mean(angs,...
            squeeze(act(:,ts)));
        PVAStren(ts) = circ_r(angs,...
            squeeze(act(:,ts)));
    end

end