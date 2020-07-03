function WidAmpDat = BumpAmpAndWidth(cond,trialName,PVAThresh) 
% WidAmpDat = BumpAmpAndWidth(cond,trialName,PVAThresh) 
%   Calculate the bump width and amplitude at each point
%
%   Input:
%     cond          a structure with experimental data
%     trialName     the type of trial to consider
%     PVAThresh     the PVA threshold to use
%
%   Output:
%     WidAmpDat 

    % Specify filtering parameters
    sgolayOrder = 3;
    sgolayFrames = 11;

    % Initialize the structure
    WidAmpDat = {};
    
    % Step through the different conditions
    for condID = 1:length(cond)
        WidAmpDat{condID}.Dark = {};
        WidAmpDat{condID}.CL = {};
        WidAmpDat{condID}.OL = {};

        % Step throught the different flies
        for flyID = 1:cond{condID}.numFlies

            % Initialize the arrays
            WidAmpDat{condID}.Dark{flyID}.allvR = [];
            WidAmpDat{condID}.Dark{flyID}.allWid = [];
            WidAmpDat{condID}.Dark{flyID}.allAmp = [];

            WidAmpDat{condID}.CL{flyID}.allvR = [];
            WidAmpDat{condID}.CL{flyID}.allWid = [];
            WidAmpDat{condID}.CL{flyID}.allAmp = [];

            WidAmpDat{condID}.OL{flyID}.allvR = [];
            WidAmpDat{condID}.OL{flyID}.allWid = [];
            WidAmpDat{condID}.OL{flyID}.allAmp = [];

            % Step through the trials
            for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
                % Pull out the data for this given trial
                datNow = cond{condID}.allFlyData{flyID}.(trialName){trialID};

                if isempty(datNow)
                    continue;
                end
                
                % Filter the rotational velocity
                vR = datNow.positionDatMatch.vRot;
                vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);

                % Calculate the PVA, max, and min
                DF = datNow.ROIaveMax;
                [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));
                DFMax = max(DF);
                DFMin = min(DF);

                % Sort the visual conditions
                [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                darkJumps = find(diff(darkPer)>1);
                if ~isempty(darkJumps)
                    darkPer(darkJumps(1)+1:end) = [];
                end

                % Step through the different visual conditions
                for visID = 1:3
                    if visID == 1
                        tRange = darkPer;
                        visCond = 'Dark';
                    elseif visID == 2
                        tRange = CLPer;
                        visCond = 'CL';
                    elseif visID == 3
                        tRange = OLPer;
                        visCond = 'OL';
                    end
                    if isempty(tRange)
                        continue;
                    end
                    if tRange(end) > length(vR)
                        tRange(end) = [];
                    end
                    
                    % Populate the rotational velocity and max value
                    WidAmpDat{condID}.(visCond){flyID}.allvR = ...
                        vertcat(WidAmpDat{condID}.(visCond){flyID}.allvR, vR(tRange));
                    WidAmpDat{condID}.(visCond){flyID}.allAmp = ...
                        vertcat(WidAmpDat{condID}.(visCond){flyID}.allAmp, DFMax(tRange)');
                    
                    % Calculate the FWHM at each time point
                    allWids = [];
                    for tStep = tRange(1):tRange(end)
                        if PVAStren(tStep+1) < PVAThresh
                            allWids = vertcat(allWids,NaN);
                            continue;
                        end
                        DFNow = DF(:,1+tStep);
                        DFNow = circshift(DFNow,8-find(DFNow == DFMax(tStep+1)));

                        widMin = find(DFNow(1:8)<(0.5*(DFMax(tStep+1)+DFMin(tStep+1))));
                        widMax = 7+find(DFNow(8:16)<(0.5*(DFMax(tStep+1)+DFMin(tStep+1))));

                        if isempty(widMin) || isempty(widMax)
                            allWids = vertcat(allWids,NaN);
                            continue;
                        end
                        widNow = pi/8*(widMax(1)-widMin(end)-1);
                        allWids = vertcat(allWids,widNow);
                    end
                    
                    % Populate the bump width
                    WidAmpDat{condID}.(visCond){flyID}.allWid = ...
                        vertcat(WidAmpDat{condID}.(visCond){flyID}.allWid, allWids);
                end
            end
        end
    end
end