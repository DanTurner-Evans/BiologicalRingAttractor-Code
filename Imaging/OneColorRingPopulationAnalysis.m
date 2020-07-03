%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\OneColor\R4d-population\';

%% Clear out old data
cond{1}.name = 'R4d';
cond{1}.dirs{1} = strcat(dataDir,'20190704');
cond{1}.dirs{2} = strcat(dataDir,'20190708');
cond{1}.dirs{3} = strcat(dataDir,'20190710');
cond{1}.dirs{4} = strcat(dataDir,'20190712');

cond = FlyDatLoad(1,cond);

%% Specify the parameters to use throughout the analyses
condID = 1;

% Filtering values
sgolayOrder = 3;
sgolayFrames = 11;

%% Plot the activity vs. stripe position or dark - Figure S8 C,D 

% Initialize the figures
ActFig = figure('units','normalized','outerposition',[0 0 1 1]);
RingZoomFig = figure('units','normalized','outerposition',[0 0 1 1]);
SummaryFig = figure('units','normalized','outerposition',[0 0 1 0.5]);
peaks = figure('units','normalized','outerposition',[0 0 0.25 0.5]);

% Initialize an array to hold the stats
allRingStats = zeros(4,cond{condID}.numFlies,10);

% Step through the flies
for flyID = 1:cond{condID}.numFlies
    figure(ActFig);
    
    % Structure to hold the activity across trials
    allRing = {};
    
    % Step through the trials
    for trialID = 1:length(cond{condID}.allFlyData{flyID}.OpenLoop)
        
        % Pull the data for this trial
        datNow = cond{condID}.allFlyData{flyID}.OpenLoop{trialID};
        
        % Get the visual conditions
        [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
        darkPerBreak = find(diff(darkPer)>1);
        darkPer(darkPerBreak(1)+1:end) = [];
        
        % Get the behavioral data
        tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
        stripePos = datNow.positionDatMatch.OffsetRotMatch(:,2);
        
        % Find when the stripe disappears behind the fly
        stripeInvis = find(abs(stripePos(OLPer))>2*pi/3);
        invisEnd = stripeInvis(find(diff(stripeInvis)>1));
        invisEnd = vertcat(invisEnd,stripeInvis(end));
        invisStart = stripeInvis(find(diff(stripeInvis)>1)+1);
        invisStart = vertcat(stripeInvis(1),invisStart);
        
        tPts = tPts-tPts(OLPer(1));
        
        % Plot the stripe position
        subplot(7,10,flyID);
        hold on;
        plot(tPts(darkPer),zeros(length(darkPer),1),'k');
%         for invis = 1:length(invisStart)
%             patch('XData',[tPts(OLPer(invisStart(invis))) tPts(OLPer(invisEnd(invis))) tPts(OLPer(invisEnd(invis))) tPts(OLPer(invisStart(invis)))],...
%                 'YData',[-pi -pi pi pi],'FaceColor','k','FaceAlpha',0.1);
%         end
        plot(tPts(CWPer),stripePos(CWPer),'b');
        plot(tPts(CCWPer),stripePos(CCWPer),'c');
        xlim([-60 60]);
        title(strcat('fly #',num2str(flyID)));
        if flyID == 1
            ylabel('stripe position (rad)');
            legend({'dark','open loop (CW)','open loop (CCW)'},'Location','northwest');
            legend('boxoff');
        end
        ylim([-pi pi]);
        yticks([-pi, 0, pi]);
        
        % Plot the activity
        subplot(7,10,flyID+10*trialID);
        hold on;
        
        imagesc(tPts,1:size(datNow.ROIaveMax,1)-1,sgolayfilt(datNow.ROIaveMax(2:end,:)-1,sgolayOrder,sgolayFrames,[],2));
%         imagesc(tPts,1:size(datNow.ROIaveMax,1)-1,datNow.ROIaveMax(2:end,:)-1);
        xlim([-60 60]);
        ylim([0.5 size(datNow.ROIaveMax,1)-0.5]);
        caxis([-0.1 1]);
        if trialID == 1
            title('glomeruli activity');
        end
        if flyID == 1
            ylabel(strcat('trial#',num2str(trialID)));
        end
        
        subplot(7,10,flyID+10*(length(cond{condID}.allFlyData{flyID}.OpenLoop)+1));
        hold on;
        plot(tPts,sgolayfilt(datNow.ROIaveMax(1,:)-1,sgolayOrder,sgolayFrames));
%         plot(tPts,datNow.ROIaveMax(1,:)-1);

        allRing{trialID,1} = tPts;
        allRing{trialID,2} = sgolayfilt(datNow.ROIaveMax(1,:)-1,sgolayOrder,sgolayFrames);
    end
    
    
    % Find the mean for the times when the fly is in the dark (neg) or when
    % there is an open-loop stripe (pos)
    ring1 = allRing{1,2};
    posMean = ring1(find(allRing{1,1}>=0));
    negMean = ring1(find(allRing{1,1}<=0));
    for trialID  = 2:length(cond{condID}.allFlyData{flyID}.OpenLoop)
        
        ringNow = allRing{trialID,2};
        
        posMeanNow = ringNow(find(allRing{trialID,1}>=0));
        if length(posMeanNow)> length(posMean)
            posMeanNow(length(posMean)+1:end) = [];
        elseif length(posMeanNow)< length(posMean)
            posMean(length(posMeanNow)+1:end) = [];
        end
        posMean = posMean + posMeanNow;
        
        negMeanNow = ringNow(find(allRing{trialID,1}<=0));
        if length(negMeanNow)> length(negMean)
            negMeanNow(1:length(negMeanNow)-length(negMean)) = [];
        elseif length(negMeanNow)< length(negMean)
            negMean(1:length(negMean)-length(negMeanNow)) = [];
        end
        negMean = negMean + negMeanNow;
    end
    posMean = posMean./trialID;
    negMean = negMean./trialID;
    
    % Plot
    for figID = 1:2
        
        % Pull the time points
        tPts = allRing{trialID,1};
        tPos = tPts(find(allRing{trialID,1}>=0));
        tNeg = tPts(find(allRing{trialID,1}<=0));
        
        
        if figID == 1
            subplot(7,10,flyID+10*(length(cond{condID}.allFlyData{flyID}.OpenLoop)+1));
        else
            figure(RingZoomFig);
            subplot(5,4,2*flyID-1);
            hold on;
            for invis = 1:length(invisStart)
                patch('XData',[tPts(OLPer(invisStart(invis))) tPts(OLPer(invisEnd(invis))) tPts(OLPer(invisEnd(invis))) tPts(OLPer(invisStart(invis)))],...
                    'YData',[-pi -pi pi pi],'FaceColor','k','FaceAlpha',0.1);
            end
        end
        
        OLMean = plot(tPos(1:length(posMean)),posMean,'b');
        darkMean = plot(tNeg(length(tNeg)-length(negMean)+1:end),negMean,'k');

        darkAv = plot([tNeg(length(tNeg)-length(negMean)+1) tPos(length(posMean))],[mean(negMean) mean(negMean)],'k','LineStyle','--');
        if flyID == 1
            legend([darkMean,OLMean,darkAv],{'mean across trials - dark','mean across trials - open loop','mean dark activity'});
            legend('boxoff');
        end
        xlim([-60 60]);
        ylim([0 0.25]);
        
        if figID == 1
            xlabel('time (s)');
            if flyID == 1
                ylabel('ring activity (DF/F)');
            end
            title('ring activity - all trials');
        else
            ylabel('ring activity (DF/F)');
            if flyID == 1
                title('ring activity - all trials');
            end
            if flyID == 9 || flyID == 10
                xlabel('time (s)');
            end
        end
        
        if figID == 2
            subplot(5,4,2*flyID);
            hold on;
            for vis = 1:5
                pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
                plot(tPos(pltRng1)-tPos(pltRng1(1)),posMean(pltRng1),'b');
                pltRng2 = [OLPer(invisStart(vis)):OLPer(invisEnd(vis))]-length(tNeg);
                plot(tPos(pltRng2)-tPos(pltRng1(1)),posMean(pltRng2),'Color','k');
            end
            for vis = 7:11
                pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
                plot(tPos(pltRng1)-tPos(pltRng1(1)),fliplr(posMean(pltRng1)),'c');
                pltRng2 = [OLPer(invisStart(vis+1)):OLPer(invisEnd(vis+1))]-length(tNeg);
                plot(tPos(pltRng2)-tPos(pltRng1(end))-(tPos(pltRng2(end))-tPos(pltRng2(1))),fliplr(posMean(pltRng2)),'Color','k');
            end
            xlim([-1.5 3.55]);
            ylim([0 0.25]);
            plot([-1.5 4],[mean(negMean) mean(negMean)],'k','LineStyle','--');
            
            figure(SummaryFig);
            subplot(2,1,1);
            hold on;
            for vis = 1:5
                pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
                plot(5.1*(flyID-1)+tPos(pltRng1)-tPos(pltRng1(1)),posMean(pltRng1)-mean(negMean),'b');
                pltRng2 = [OLPer(invisStart(vis)):OLPer(invisEnd(vis))]-length(tNeg);
                plot(5.1*(flyID-1)+tPos(pltRng2)-tPos(pltRng1(1)),posMean(pltRng2)-mean(negMean),'Color','k');
            end
            for vis = 7:11
                pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
                plot(5.1*(flyID-1)+tPos(pltRng1)-tPos(pltRng1(1)),fliplr(posMean(pltRng1))-mean(negMean),'c');
                pltRng2 = [OLPer(invisStart(vis+1)):OLPer(invisEnd(vis+1))]-length(tNeg);
                plot(5.1*(flyID-1)+tPos(pltRng2)-tPos(pltRng1(end))-(tPos(pltRng2(end))-tPos(pltRng2(1))),fliplr(posMean(pltRng2))-mean(negMean),'Color','k');
            end
            xlim([-1.5 3.55+5.1*9]);
            ylim([-0.25 0.25]);
            plot([5.1*(flyID-1)-1.5 4+5.1*(flyID-1)],[0 0],'k','LineStyle','--');
            
            % Calculate and plot the max and min activity in the different periods
            figure(peaks)        
            for vis = 1:5
                subplot(1,2,1);
                hold on;
                pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
                scatter(flyID,max(posMean(pltRng1))-mean(negMean),20,[0 1 0.5],'filled');
                scatter(flyID,min(posMean(pltRng1))-mean(negMean),20,[1 0 0.5],'filled');
                allRingStats(1,flyID,vis) = max(posMean(pltRng1))-mean(negMean);
                allRingStats(2,flyID,vis) = min(posMean(pltRng1))-mean(negMean);
                
                subplot(1,2,2);
                hold on;
                pltRng2 = [OLPer(invisStart(vis)):OLPer(invisEnd(vis))]-length(tNeg);
                scatter(flyID,max(posMean(pltRng2))-mean(negMean),20,[0 1 0.5],'filled');
                scatter(flyID,min(posMean(pltRng2))-mean(negMean),20,[1 0 0.5],'filled');
                allRingStats(3,flyID,vis) = max(posMean(pltRng2))-mean(negMean);
                allRingStats(4,flyID,vis) = min(posMean(pltRng2))-mean(negMean);
            end
            for vis = 7:11
                subplot(1,2,1);
                pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
                scatter(flyID,max(posMean(pltRng1))-mean(negMean),20,[0 1 0.5],'filled');
                scatter(flyID,min(posMean(pltRng1))-mean(negMean),20,[1 0 0.5],'filled');
                allRingStats(1,flyID,vis-1) = max(posMean(pltRng1))-mean(negMean);
                allRingStats(2,flyID,vis-1) = min(posMean(pltRng1))-mean(negMean);

                subplot(1,2,2);
                pltRng2 = [OLPer(invisStart(vis+1)):OLPer(invisEnd(vis+1))]-length(tNeg);
                scatter(flyID,max(posMean(pltRng2))-mean(negMean),20,[0 1 0.5],'filled');
                scatter(flyID,min(posMean(pltRng2))-mean(negMean),20,[1 0 0.5],'filled');
                allRingStats(3,flyID,vis-1) = max(posMean(pltRng2))-mean(negMean);
                allRingStats(4,flyID,vis-1) = min(posMean(pltRng2))-mean(negMean);
            end
        end
    end
end
figure(peaks);
subplot(1,2,1);
title('stripe visible');
xlim([0 11]);
ylim([-0.25 0.25]);
xlabel('fly #');
ylabel('\Delta DF/F');
plot([0 11],[0 0],'Color','k');
alpha(0.4);

subplot(1,2,2);
title('stripe behind');
xlim([0 11]);
ylim([-0.25 0.25]);
xlabel('fly #');
plot([0 11],[0 0],'Color','k');
alpha(0.4);
    
figure(ActFig);
colormap(brewermap(64, 'Blues'));

%% Run the above analyses, now for the glomeruli

% Specify random colors for the glomeruli
glomCols = rand(8,3);

GlomAct = figure('units','normalized','outerposition',[0 0 1 1]);
peaks = figure('units','normalized','outerposition',[0 0 0.25 0.5]);

allGlomStats = cell(cond{condID}.numFlies,1);
for flyID = 1:cond{condID}.numFlies
    
    allGlom = {};
    for trialID = 1:length(cond{condID}.allFlyData{flyID}.OpenLoop)
        
        datNow = cond{condID}.allFlyData{flyID}.OpenLoop{trialID};
        
        [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
        darkPerBreak = find(diff(darkPer)>1);
        darkPer(darkPerBreak(1)+1:end) = [];
        
        tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
        stripePos = datNow.positionDatMatch.OffsetRotMatch(:,2);
        
        stripeInvis = find(abs(stripePos(OLPer))>2*pi/3);
        invisEnd = stripeInvis(find(diff(stripeInvis)>1));
        invisEnd = vertcat(invisEnd,stripeInvis(end));
        invisStart = stripeInvis(find(diff(stripeInvis)>1)+1);
        invisStart = vertcat(stripeInvis(1),invisStart);
        
        tPts = tPts-tPts(OLPer(1));
        allGlom{trialID,1} = tPts;
        
        figure(GlomAct);
        subplot(5,4,2*flyID-1);
        hold on;
        for glom = 2:size(datNow.ROIaveMax,1)
%             plot(tPts,sgolayfilt(datNow.ROIaveMax(glom,:)-1,sgolayOrder,sgolayFrames),'Color',glomCols(glom-1,:));
            allGlom{trialID,glom} = sgolayfilt(datNow.ROIaveMax(glom,:)-1,sgolayOrder,sgolayFrames);
        end
        xlim([tPts(1) tPts(end)]);
        ylim([0 1.5]);
    end
    
    glomPosMean = {};
    glomNegMean = {};
    
    flyGlomStats = zeros(size(datNow.ROIaveMax,1)-1,4,10);
    for glom = 2:size(datNow.ROIaveMax,1)
        glomNow = allGlom{1,glom};
        posMean = glomNow(find(allGlom{1,1}>=0));
        negMean = glomNow(find(allGlom{1,1}<=0));

        for trialID  = 2:length(cond{condID}.allFlyData{flyID}.OpenLoop)

            glomNow = allGlom{trialID,glom};

            posMeanNow = glomNow(find(allGlom{trialID,1}>=0));
            if length(posMeanNow)> length(posMean)
                posMeanNow(length(posMean)+1:end) = [];
            elseif length(posMeanNow)< length(posMean)
                posMean(length(posMeanNow)+1:end) = [];
            end
            posMean = posMean + posMeanNow;

            negMeanNow = glomNow(find(allGlom{trialID,1}<=0));
            if length(negMeanNow)> length(negMean)
                negMeanNow(1:length(negMeanNow)-length(negMean)) = [];
            elseif length(negMeanNow)< length(negMean)
                negMean(1:length(negMean)-length(negMeanNow)) = [];
            end
            negMean = negMean + negMeanNow;
        end
        posMean = posMean./trialID;
        negMean = negMean./trialID;
        
        figure(GlomAct);
        subplot(5,4,2*flyID-1);
        tPts = allGlom{trialID,1};
        tPos = tPts(find(allGlom{trialID,1}>=0));
        tNeg = tPts(find(allGlom{trialID,1}<=0));
        OLMean = plot(tPos(1:length(posMean)),posMean,'b');
        darkMean = plot(tNeg(length(tNeg)-length(negMean)+1:end),negMean,'k');
        
        subplot(5,4,2*flyID);
        hold on;
        for vis = 1:5
            pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
            plot(tPos(pltRng1)-tPos(pltRng1(1)),posMean(pltRng1),'b');
            pltRng2 = [OLPer(invisStart(vis)):OLPer(invisEnd(vis))]-length(tNeg);
            plot(tPos(pltRng2)-tPos(pltRng1(1)),posMean(pltRng2),'Color','k');
        end
        for vis = 7:11
            pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
            plot(tPos(pltRng1)-tPos(pltRng1(1)),fliplr(posMean(pltRng1)),'c');
            pltRng2 = [OLPer(invisStart(vis+1)):OLPer(invisEnd(vis+1))]-length(tNeg);
            plot(tPos(pltRng2)-tPos(pltRng1(end))-(tPos(pltRng2(end))-tPos(pltRng2(1))),fliplr(posMean(pltRng2)),'Color','k');
        end
        xlim([-1.5 3.55]);
        ylim([0 1.5]);
        plot([-1.5 4],[mean(negMean) mean(negMean)],'k','LineStyle','--');
        
        glomPosMean{glom-1} = posMean;
        glomNegMean{glom-1} = negMean;
        
        figure(peaks)        
        for vis = 1:5
            subplot(1,2,1);
            hold on;
            pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
            scatter(flyID,max(posMean(pltRng1))-mean(negMean),20,[0 1 0.5],'filled');
            scatter(flyID,min(posMean(pltRng1))-mean(negMean),20,[1 0 0.5],'filled');
            flyGlomStats(glom-1,1,vis) = max(posMean(pltRng1))-mean(negMean);
            flyGlomStats(glom-1,2,vis) = min(posMean(pltRng1))-mean(negMean);
            
            subplot(1,2,2);
            hold on;
            pltRng2 = [OLPer(invisStart(vis)):OLPer(invisEnd(vis))]-length(tNeg);
            scatter(flyID,max(posMean(pltRng2))-mean(negMean),20,[0 1 0.5],'filled');
            scatter(flyID,min(posMean(pltRng2))-mean(negMean),20,[1 0 0.5],'filled');
            flyGlomStats(glom-1,3,vis) = max(posMean(pltRng2))-mean(negMean);
            flyGlomStats(glom-1,4,vis) = min(posMean(pltRng2))-mean(negMean);
        end
        for vis = 7:11
            subplot(1,2,1);
            pltRng1 = [OLPer(invisEnd(vis)):OLPer(invisStart(vis+1))]-length(tNeg);
            scatter(flyID,max(posMean(pltRng1))-mean(negMean),20,[0 1 0.5],'filled');
            scatter(flyID,min(posMean(pltRng1))-mean(negMean),20,[1 0 0.5],'filled');
            flyGlomStats(glom-1,1,vis-1) = max(posMean(pltRng1))-mean(negMean);
            flyGlomStats(glom-1,2,vis-1) = min(posMean(pltRng1))-mean(negMean);
            
            subplot(1,2,2);
            pltRng2 = [OLPer(invisStart(vis+1)):OLPer(invisEnd(vis+1))]-length(tNeg);
            scatter(flyID,max(posMean(pltRng2))-mean(negMean),20,[0 1 0.5],'filled');
            scatter(flyID,min(posMean(pltRng2))-mean(negMean),20,[1 0 0.5],'filled');
            flyGlomStats(glom-1,3,vis-1) = max(posMean(pltRng2))-mean(negMean);
            flyGlomStats(glom-1,4,vis-1) = min(posMean(pltRng2))-mean(negMean);
        end
    end
    allGlomStats{flyID} = flyGlomStats;
end

figure(peaks);
subplot(1,2,1);
title('stripe visible');
xlim([0 11]);
ylim([-0.5 2]);
xlabel('fly #');
ylabel('\Delta DF/F');
plot([0 11],[0 0],'Color','k');
alpha(0.4);

subplot(1,2,2);
title('stripe behind');
xlim([0 11]);
ylim([-0.5 2]);
xlabel('fly #');
plot([0 11],[0 0],'Color','k');
alpha(0.4);

%% Plot the mean of the max and min across flies - Figure S8 E 

% Make the figures
figure;
glomMaxMeansDark = [];
glomMaxMeansStripeVis = [];
glomMinMeansDark = [];
glomMinMeansStripeVis = [];

% Step through the flies
for flyID = 1:cond{condID}.numFlies
    
    % Plot the stats
    subplot(1,2,1);
    hold on;
    scatter(1,mean(allRingStats(1,flyID,:)),40,[1 0 0.5],'filled');
    scatter(2,mean(allRingStats(3,flyID,:)),40,[1 0 0.5],'filled');
    plot([1 2],[mean(allRingStats(1,flyID,:)) mean(allRingStats(3,flyID,:))],'Color',[1 0 0.5]);
    
    scatter(1,mean(allRingStats(2,flyID,:)),40,[0 0.5 1],'filled');
    scatter(2,mean(allRingStats(4,flyID,:)),40,[0 0.5 1],'filled');
    plot([1 2],[mean(allRingStats(2,flyID,:)) mean(allRingStats(4,flyID,:))],'Color',[0 0.5 1]);
    
    glomStatsNow = allGlomStats{flyID};
    for glom = 1:size(glomStatsNow,1);        
        subplot(1,2,2);
        hold on;
        scatter(1,mean(glomStatsNow(glom,1,:)),40,[1 0 0.5],'filled');
        scatter(2,mean(glomStatsNow(glom,3,:)),40,[1 0 0.5],'filled');
        plot([1 2],[mean(glomStatsNow(glom,1,:)) mean(glomStatsNow(glom,3,:))],'Color',[1 0 0.5]);
        glomMaxMeansDark = [glomMaxMeansDark mean(glomStatsNow(glom,3,:))];
        glomMaxMeansStripeVis = [glomMaxMeansStripeVis mean(glomStatsNow(glom,1,:))];
        
        scatter(1,mean(glomStatsNow(glom,2,:)),40,[0 0.5 1],'filled');
        scatter(2,mean(glomStatsNow(glom,4,:)),40,[0 0.5 1],'filled');
        plot([1 2],[mean(glomStatsNow(glom,2,:)) mean(glomStatsNow(glom,4,:))],'Color',[0 0.5 1]);
        glomMinMeansDark = [glomMinMeansDark mean(glomStatsNow(glom,2,:))];
        glomMinMeansStripeVis = [glomMinMeansStripeVis mean(glomStatsNow(glom,4,:))];
    end
end

subplot(1,2,1);
xlim([0 3]);
xticks([1 2])
xticklabels({'stripe visible','stripe behind'});
ylim([-0.15 0.15]);
ylabel('\Delta DF/F');
plot([0 3],[0 0],'Color','k');
alpha(0.4);

subplot(1,2,2);
xlim([0 3]);
xticks([1 2])
xticklabels({'stripe visible','stripe behind'});
ylim([-0.5 2]);
ylabel('\Delta DF/F');
plot([0 3],[0 0],'Color','k');
alpha(0.4);

% Calculate p values
[h,p] = ttest(mean(allRingStats(1,:,:),3),mean(allRingStats(3,:,:),3)) % max
[h,p] = ttest(mean(allRingStats(2,:,:),3),mean(allRingStats(4,:,:),3)) % min

[h,p] = ttest(glomMaxMeansDark, glomMaxMeansStripeVis) % min
[h,p] = ttest(glomMinMeansDark, glomMinMeansStripeVis) % min