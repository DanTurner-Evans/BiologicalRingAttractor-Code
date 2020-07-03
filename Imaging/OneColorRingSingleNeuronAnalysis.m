%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\OneColor\R4d-single\';

%% Load the data
cond{1}.name = 'RN-FLP';
cond{1}.dirs{1} = strcat(dataDir,'20171103');
cond{1}.dirs{2} = strcat(dataDir,'20171106');
cond{1}.dirs{3} = strcat(dataDir,'20171107');
cond{1}.dirs{4} = strcat(dataDir,'20171108');
cond{1}.dirs{5} = strcat(dataDir,'20171129');
cond{1}.dirs{6} = strcat(dataDir,'20171205');

cond = FlyDatLoad(2,cond);

%% Find the total number of ring neurons
total_RNs = 0;
for flyID = 1:cond{1}.numFlies
    for trialID = 1
        if ndims(cond{1}.allFlyData{flyID}.All{trialID}.GROIaveMax) == 2
            total_RNs = total_RNs + 1;
        else
            total_RNs = total_RNs +...
                size(cond{1}.allFlyData{flyID}.All{trialID}.GROIaveMax,1);
        end
    end
end

%% For the open loop part of the trials, average the mean activity vs. position - Figure S8 A

act = figure('units','normalized','outerposition',[0 0 1 1]);

% Specify the current ring neuron
RNNow = 1;

% Discretize the stripe position
inpRng = linspace(-pi,pi,43);

% Create a matrix to hold the 1D receptive fields

% Step through the flies
RFs = zeros(total_RNs,length(inpRng));
for flyID = 1:cond{1}.numFlies
    if ndims(cond{1}.allFlyData{flyID}.All{1}.GROIaveMax) == 2
        num_RNs = 1;
    else
        num_RNs = size(cond{1}.allFlyData{flyID}.All{1}.GROIaveMax,1);
    end
    
    % Step through the ring neurons
    for RN = 1:num_RNs
        allCWRot = {};
        allCWAct = {};
        allCCWRot = {};
        allCCWAct = {};
        
        % Step through the trials
        for trialID = 1:length(cond{1}.allFlyData{flyID}.All)
        
            % Get the stripe position
            RotNow = cond{1}.allFlyData{flyID}.All{trialID}.positionDatMatch.OffsetRotMatch(:,2);
            
            % Find the open loop period
            OLPer = find(cond{1}.allFlyData{flyID}.All{trialID}.positionDatMatch.OLGain > 0);
            % Find the clockwise rotations
            CWPer = find(cond{1}.allFlyData{flyID}.All{trialID}.positionDatMatch.Direction > 0);
            CWPer = intersect(OLPer,CWPer);
            posCWPer = find(RotNow(CWPer)>=0);
            stepCWPer = diff(posCWPer);
            startPtsCW = posCWPer(find(stepCWPer>1))+1;
            endPtsCW = vertcat(startPtsCW(2:end)-1,posCWPer(end));
            
            % Find the counter-clockwise rotations
            CCWPer = find(cond{1}.allFlyData{flyID}.All{trialID}.positionDatMatch.Direction < 0);
            CCWPer = intersect(OLPer,CCWPer);
            negCCWPer = find(RotNow(CCWPer)<=0);
            stepCCWPer = diff(negCCWPer);
            startPtsCCW = negCCWPer(find(stepCCWPer>1))+1;
            endPtsCCW = vertcat(startPtsCCW(2:end)-1,negCCWPer(end));
            
            % Find the mean ring activity
            if num_RNs == 1
                meanRingAct = mean(squeeze(cond{1}.allFlyData{flyID}.All{trialID}.GROIaveMax),1);
            else
                meanRingAct = mean(squeeze(cond{1}.allFlyData{flyID}.All{trialID}.GROIaveMax(RN,:,:)),1);
            end

            % Find the CW activity
            CWRot = zeros(length(startPtsCW),max(startPtsCW-endPtsCW));
            CWAct = zeros(length(startPtsCW),max(startPtsCW-endPtsCW));
            for iter = 1:length(startPtsCW)
                CWRot(iter,1:(endPtsCW(iter)-startPtsCW(iter))+1) = RotNow(CWPer(startPtsCW(iter):endPtsCW(iter)));
                CWAct(iter,1:(endPtsCW(iter)-startPtsCW(iter))+1) = meanRingAct(CWPer(startPtsCW(iter):endPtsCW(iter)));
            end
            allCWRot{trialID} = CWRot;
            allCWAct{trialID} = CWAct;
            
            % Find the CCW activity
            CCWRot = zeros(length(startPtsCCW),max(startPtsCCW-endPtsCCW));
            CCWAct = zeros(length(startPtsCCW),max(startPtsCCW-endPtsCCW));
            for iter = 1:length(startPtsCCW)
                CCWRot(iter,1:(endPtsCCW(iter)-startPtsCCW(iter))+1) = RotNow(CCWPer(startPtsCCW(iter):endPtsCCW(iter)));
                CCWAct(iter,1:(endPtsCCW(iter)-startPtsCCW(iter))+1) = meanRingAct(CCWPer(startPtsCCW(iter):endPtsCCW(iter)));
            end
            allCCWRot{trialID} = CCWRot;
            allCCWAct{trialID} = CCWAct;
            
        end
        
        % Find the mean RF
        inpValsCW = zeros(length(allCWRot)*5,43);
        inpValsCCW = zeros(length(allCCWRot)*5,43);
        
        fillNow = 1;
        
        for trialID = 1:length(allCWRot)
            if ~isempty(allCWRot{trialID})
                for iter = 1:5
                    posNowCW = allCWRot{trialID}(iter,:);
                    posNowCCW = allCCWRot{trialID}(iter,:);
                    actNowCW = allCWAct{trialID}(iter,:);
                    actNowCCW = allCCWAct{trialID}(iter,:);
                    if posNowCW(end) == 0
                        posNowCW(end) = [];
                        actNowCW(end) = [];
                    end
                    if posNowCCW(end) == 0
                        posNowCCW(end) = [];
                        actNowCCW(end) = [];
                    end
                    inpValsCW(iter+5*(fillNow-1),:) = interp1(posNowCW,actNowCW,inpRng);
                    inpValsCCW(iter+5*(fillNow-1),:) = interp1(posNowCCW,actNowCCW,inpRng);
                end
                fillNow = fillNow+1;
            else
                inpValsCW(([1:5]+5*(fillNow-1)),:) = [];
                inpValsCCW(([1:5]+5*(fillNow-1)),:) = [];
            end
        end
        subplot(ceil(sqrt(total_RNs)),ceil(sqrt(total_RNs)),RNNow);
        hold on;
        plot(inpRng,mean(inpValsCW,1)-1,'b');
        plot(inpRng,mean(inpValsCCW,1)-1,'m');
        plot(inpRng,0.5*(mean(inpValsCW,1)+mean(inpValsCCW,1))-1,'k','LineWidth',2);
        xlim([-pi pi]);
        ylim([0 1.5]);
        xlabel('pos (rad)');
        ylabel('mean act (DF/F)');
        patch([-pi -pi*2/3 -pi*2/3 -pi],[0 0 1.5 1.5],'k','FaceAlpha',0.1);
        patch([pi pi*2/3 pi*2/3 pi],[0 0 1.5 1.5],'k','FaceAlpha',0.1);
        if RNNow == 1
            legend({'CW','CCW','mean'});
            legend('boxoff');
        end
        
        RFs(RNNow,:) = 0.5*(mean(inpValsCW,1)+mean(inpValsCCW,1))-1;
        RNNow = RNNow + 1;
    end
end

% set(act,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
% print(act,'C:\Users\turnerevansd\Documents\Ring\FLP\Plots\MeanRFs','-dpdf');