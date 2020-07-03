%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\Shi\EPGs\';

%% Clear out old data and load the new data
cond{1}.name = 'Shi and GCaMP in E-PGs';
cond{1}.dirs{1} = strcat(dataDir,'20180122');
cond{1}.dirs{2} = strcat(dataDir,'20180123');
cond{1}.dirs{3} = strcat(dataDir,'20180124');
cond{1}.dirs{4} = strcat(dataDir,'20180223');
cond{1}.dirs{5} = strcat(dataDir,'20180225');
cond{1}.dirs{6} = strcat(dataDir,'20180226');

cond = FlyDatLoad(1,cond);

%% Specify parameters to use throughout the analyses
condNow = 1;

% Thresholds to determine when the fly is moving
vRThresh = pi/10;
vFThresh = 0.1;

% Set S-G filter parameters for the velocity
sgolayOrder = 3;
sgolayFrames = 11;

%% Activity plots, including those from Figure 3F

% Choose a length of time to plot
t2Plt = 60;

% Specify the flies and trials to show
fly2plt{1} = [1:5,7:9]
fly2plt{2} = [3,5,11]
coldTrials{1} = [2,1,2,2,1,0,1,1,1];
coldTrials{2} = [0,0,2,0,1,0,0,0,0,0,2];
hotTrials{1} = [2,2,2,2,3,0,2,2,2];
hotTrials{2} = [0,0,2,0,1,0,0,0,0,0,2];

% Find the angles of the ROIs
angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs + mean(diff(angs));

% Initialize the figures
coldFig = figure('units','normalized','outerposition',[0 0 1 1]);
hotFig = figure('units','normalized','outerposition',[0 0 1 1]);
coldPltID = 0;
hotPltID = 0;

% Step through the flies
for flyID = 1:cond{condID}.numFlies
    if ~ismember(flyID,fly2plt{condID})
        continue;
    end

    % Step through the temperatures
    for hcID = 1:2
        if hcID == 1
            if condID == 1 & (flyID > 4)
                trialName = 'Stripe';
            else
                trialName = 'All';
            end
        elseif hcID == 2
            if condID == 1 & (flyID > 4)
                trialName = 'Stripe_30C';
            else
                trialName = 'All_30C';
            end
        end

        % Step through the trials
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
            if hcID == 1
                if trialID ~= coldTrials{condID}(flyID)
                    continue;
                end
                coldPltID = coldPltID + 1;
                figure(coldFig);
                pltID = coldPltID;
            else
                if trialID ~= hotTrials{condID}(flyID)
                    continue;
                end
                hotPltID = hotPltID + 1;
                figure(hotFig);
                pltID = hotPltID;
            end

            % Pull the data
            datNow = cond{condID}.allFlyData{flyID}.(trialName){trialID};

            % Find the different visual conditions
            [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);

            % Get the calcium data
            DF = datNow.ROIaveMax-1;

            % Get the time points
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            tPts = tPts - tPts(1);

            % Get the heading
            heading = datNow.positionDatMatch.OffsetRotMatch(:,2);

            % Calculate the PVA
            [angs, PVAPlt, PVAStren] = PVA(DF);
            for tStep = 2:length(PVAPlt)
                if abs(PVAPlt(tStep)-PVAPlt(tStep-1))>pi
                    PVAPlt(tStep-1) = NaN;
                end
            end
            PVAWeak = find(PVAStren<0.1);
            PVAPlt(PVAWeak) = NaN;


            % Plot the activity during the closed loop period
            tPlt = tPts(CLPer);
            tPlt = tPlt - tPlt(1);
            if condID == 1
                subplot(4,2,pltID)
            else
                subplot(4,2,2-mod(pltID,2)+4*floor((pltID-1)/2));
            end
            hold on;
            imagesc(tPlt,angs,DF(:,CLPer))
            xlim([tPlt(1) t2Plt]);
            ylim([-pi pi]);
            ylabel('EB act. (rad)');
            set(gca,'YTick',[-pi, 0, pi],'YTickLabels',{'-p','0','p'});
            if condID == 1
                colorbar;
                xlabel('time (s)');
                caxis([0 2]);
            else
                plot(tPlt,PVAPlt(CLPer),'Color',[0.25 0.25 0.5]);
                caxis([0 1.25]);
            end

            if condID == 2
                subplot(4,2,4-mod(pltID,2)+4*floor((pltID-1)/2));
                hold on;
                plot(tPlt,heading(CLPer));
                plot(tPlt,PVAPlt(CLPer),'Color',[0.25 0.25 0.5]);
                xlim([tPlt(1) t2Plt]);
                xlabel('time (s)');
                ylim([-pi pi]);
                ylabel('heading (rad)')
                set(gca,'YTick',[-pi, 0, pi],'YTickLabels',{'-p','0','p'});
            end

        end
        colormap(brewermap(64, 'Blues'));

        if condID == 2 & hotPltID== 3
            subplot(4,2,6);
            colorbar;
            caxis([0 1.25]);
        end
    end
end

%% PVA Strength plot from Figure 3G

% Initialize the figure
PVAStrenSummary = figure('units','normalized','outerposition',[0 0 0.25 0.5]);
subplot(1,3,1);
hold on;

% Step through the flies
for flyID = 1:cond{condNow}.numFlies
    
    % Get the trial names
    trialNames = fields(cond{condNow}.allFlyData{flyID});
    
    % Step through the trials
    for trialType = 2:3
        for trialID = 1:min(3,length(cond{condNow}.allFlyData{flyID}.(trialNames{trialType})))
            
            % Get the data for this trial
            datNow = cond{condNow}.allFlyData{flyID}.(trialNames{trialType}){trialID};
            
            % Calculate the PVA
            [angs, PVAPlt, PVAStren] = PVA(datNow.ROIaveMax-min(min(datNow.ROIaveMax)));
            
            % Get the behavioral data
            vR = datNow.positionDatMatch.vRot;
            vR = sgolayfilt(vR,3,11);
            vF = datNow.positionDatMatch.vF;
            vF = sgolayfilt(vF,3,11);
            flyAct = union(find(abs(vR)>vRThresh),find(vF>vFThresh));
            
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            tPts = tPts-tPts(1);
            
            if length(strsplit(trialNames{trialType},'_'))==2
                if trialID >=2
                    scatter(flyID-1,mean(PVAStren(flyAct)),150,'r','filled');
                end
            else
                scatter(flyID-1,mean(PVAStren(flyAct)),150,'b','filled')
            end
            xlim([0 11]);
            ylim([0 0.5]);
            ax = gca;
            xticks([1:10]);
            xlabel('fly #');
            ylabel('mean PVA strength');
        end
    end
end

%% Plots for Figure S5;

% Specify the fly of interest
flyID2Plt = 3;
trialIDs = [2,2];

% Create a vector to hold the mean EB activity
meanEBact = zeros(2,2,10);

% Bin the velocities
vREdges = linspace(0,pi,17);
vFEdges = linspace(0,2,17);

% Plot example data
actProf = figure('units','normalized','outerposition',[0 0 1 1]);   
for hc = 1:2
    for flyID = flyID2Plt
    
        % Get the trial names
        trialNames = fields(cond{1}.allFlyData{flyID});
        trialNow = trialIDs(hc);
        if hc == 1
            trialName = trialNames{3};
        else
            trialName = trialNames{2};
        end

        % Step through the second and third trials
        for trialID = trialNow

            % Pull out the data for this trial
            datNow = cond{condNow}.allFlyData{flyID}.(trialName){trialID};

            % Get the different visual conditions used in the trial
            [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
            if ~isempty(find(diff(darkPer)>1))
                darkPer(find(diff(darkPer)>1)+1:end) = [];
            end

            % Get the behavioral data
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            vR = datNow.positionDatMatch.vRot;
            vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);

            % Find the mean activity
            meanAct = mean(datNow.ROIaveMax,1)-1;

            % Plot it!
            subplot(5,4,8*(hc-1) +[1:2]);
            hold on;
            plot(tPts(darkPer),meanAct(darkPer),'k');
            plot(tPts(CLPer),meanAct(CLPer),'b');
            line([tPts(darkPer(1)) tPts(darkPer(end))],...
                [mean(meanAct(darkPer)) mean(meanAct(darkPer))],...
                'Color','k');
            line([tPts(CLPer(1)) tPts(CLPer(end))],...
                [mean(meanAct(CLPer)) mean(meanAct(CLPer))],...
                'Color','b');
            xlabel('time (s)');
            ylabel('mean DF/F ');
            ylim([0 2.5]);
            xlim([tPts(1) tPts(end)]);
            if hc == 1
                if flyID == 1
                    legend({'dark','CL stripe'});
                    legend('boxoff');
                end
            end
            
            subplot(5,4,8*(hc-1) +[5:6]);
            hold on;
            plot(tPts(2:end),abs(vR),'Color',[0.5 0 0]);
            xlabel('time (s)');
            ylabel('|vR| (rad/s) ');
            ylim([0 pi]);
            xlim([tPts(1) tPts(end)]);
        end
    end
end

% Get the mean activity over all trials for all flies
for hc = 1:2
    for flyID = 1:cond{condNow}.numFlies
    
        % Get the trial names
        trialNames = fields(cond{1}.allFlyData{flyID});
        if hc == 1
            trialName = trialNames{3};
            trialStart = 1;
            trialEnd = 2;
        else
            trialName = trialNames{2};
            trialStart = 2;
            trialEnd = 3;
        end

        
        actDark = [];
        actCL = [];
        % Step through the second and third trials
        for trialID = trialStart:trialEnd

            % Pull out the data for this trial
            datNow = cond{condNow}.allFlyData{flyID}.(trialName){trialID};

            % Find the visual conditions used in the trial
            [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
            if ~isempty(find(diff(darkPer)>1))
                darkPer(find(diff(darkPer)>1)+1:end) = [];
            end

            % Get the behavioral data
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            vR = datNow.positionDatMatch.vRot;
            vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);
            vF = datNow.positionDatMatch.vF;
            vF = sgolayfilt(vF,sgolayOrder,sgolayFrames);
            
            % Get the periods where the fly is active
            flyAct = union(find(vR>vRThresh),find(vF>vFThresh));
            darkPer = intersect(darkPer,flyAct);
            CLPer = intersect(CLPer,flyAct);
           
            % Get the mean activity
            meanAct = mean(datNow.ROIaveMax,1)-1;
            
            actDark = [actDark meanAct(darkPer)];
            actCL = [actCL meanAct(CLPer)];        

        end
        meanEBact(hc,1,flyID) = mean(actDark);
        meanEBact(hc,2,flyID) = mean(actCL);
    end
end

% Plot the mean activity
for hc = 1:2
    for flyID = 1:cond{condNow}.numFlies
        subplot(5,4,8*(hc-1) +[3 7]);
        hold on;
        plot([1 2],[meanEBact(hc,1,flyID) meanEBact(hc,2,flyID)],'k');
        scatter(1,meanEBact(hc,1,flyID),50,'k','filled');
        scatter(2,meanEBact(hc,2,flyID),50,'b','filled');
        alpha(0.4);
        xlim([0.5 2.5]);
        xticks([1 2]);
        ylim([0 1]);
        ylabel('mean EB activity (DF/F)');
        
        if hc == 1
            title('permissive temp');
        else
            title('restrictive temp');
        end
    end
end

% Find the p values
[h,p] = ttest(meanEBact(1,1,:),meanEBact(1,2,:))
[h,p] = ttest(meanEBact(2,1,:),meanEBact(2,2,:))
[h,p] = ttest(meanEBact(1,1,:)-meanEBact(1,2,:),...
    meanEBact(2,1,:)-meanEBact(2,2,:))

% Look at the activity crosscorrelation with velocity
% Step through the flies
for flyID = 1:cond{condNow}.numFlies
    
    trialNames = fields(cond{1}.allFlyData{flyID});
    for hc = 1:2
        if hc == 1
            trialName = trialNames{3};
            trialStart = 1;
            trialEnd = 2;
        else
            trialName = trialNames{2};
            trialStart = 2;
            trialEnd = 3;
        end

        % Step through the second and third trials
        for trialID = trialStart:trialEnd

            datNow = cond{condNow}.allFlyData{flyID}.(trialName){trialID};

            % Pull out necessary parameters
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            vR = datNow.positionDatMatch.vRot;
            vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);
            vF = datNow.positionDatMatch.vF;
            vF = sgolayfilt(vF,sgolayOrder,sgolayFrames);
            meanAct = mean(cond{condNow}.allFlyData{flyID}.(trialName){trialID}.ROIaveMax,1)-1;

            % Look at the crosscorrelation
            [tLags, CCVals] = autoCCs(tPts,abs(vR),meanAct(2:end),71);
            
            if flyID == flyID2Plt
                subplot(5,4,[4 8]);
                hold on;
                if hc == 1
                    plot(tLags,CCVals,'b');
                else
                    plot(tLags,CCVals,'r');
                end
                ylim([-0.25 1]);
                xlim([tLags(1) tLags(end)]);
                xlabel('time (s)');
                ylabel('cross correlation');
                plot([tLags(1) tLags(end)],[0 0],'k');
                plot([0 0],[-0.25 1],'k');
            end

            subplot(5,4,12);
            hold on;
            if hc == 1
                scatter(flyID,max(CCVals),20,'b','filled');
            else
                scatter(flyID,max(CCVals),20,'r','filled');
            end
            ylim([0 1]);
            xlim([0 11]);
            alpha(0.4);
            xlabel('fly ID');
            ylabel('CC');
            
            subplot(5,4,16);
            hold on;
            maxPos = find(CCVals == max(CCVals));
            lowCCs = find(CCVals<=0.5*(max(CCVals)+min(CCVals)));
            FWHM = max(diff(lowCCs))*mean(diff(tLags));
            
            if hc == 1
                scatter(flyID,FWHM,20,'b','filled');
            else
                scatter(flyID,FWHM,20,'r','filled');
            end
%             ylim([0 1]);
            xlim([0 11]);
            alpha(0.4);
            xlabel('fly ID');
            ylabel('FWHM');
            
        end
    end
end