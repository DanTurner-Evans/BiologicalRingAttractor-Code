function [angs, PVA, PVAStren] = PVA(act)

    % Plot the heading vs. the stripe position
    % Find the PVA
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