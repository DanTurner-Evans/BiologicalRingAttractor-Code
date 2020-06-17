%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\Kir\PEN2\';

%% Load the data
cond{1}.name = 'Raised at 21C';
cond{1}.dirs{1} = strcat(dataDir,'Line1-Ctrl\20190822');
cond{1}.dirs{2} = strcat(dataDir,'Line1-Ctrl\20190826');

cond{2}.name = 'VT20739, Raised at 30C';
cond{2}.dirs{1} = strcat(dataDir,'Line1\20190819');
cond{2}.dirs{2} = strcat(dataDir,'Line1\20190823');

% % cond{3}.name = 'GR12D09, Raised at 21C';
cond{1}.dirs{3} = strcat(dataDir,'Line2-Ctrl\20191127');
cond{1}.dirs{4} = strcat(dataDir,'Line2-Ctrl\20191204');

% cond{4}.name = 'GR12D09, Raised at 30C';
cond{2}.dirs{3} = strcat(dataDir,'Line2\20191115');
cond{2}.dirs{4} = strcat(dataDir,'Line2\20191118');
cond{2}.dirs{5} = strcat(dataDir,'Line2\20191122');
cond{2}.dirs{6} = strcat(dataDir,'Line2\20191204');

cond = FlyDatLoad(2,cond);

%% Look at the mean DF/F, max DF/F, and PVA strength vs. vR - Figure 7 L,M, Figure S11 E,F,H,H,K

% Set analysis parameters
trialName = 'All';
PVAThresh = 0.075;
sgolayOrder = 3;
sgolayFrames = 11;

% Set the bins for the velocity histograms
vREdges = linspace(0,pi,11);
vRCents = vREdges+0.5*mean(diff(vREdges));
vRCents(end) = [];

% Calculate the mean DF/F, max DF/F, and PVA strength
actPlt = figure('units','normalized','outerposition',[0 0 1 1]);
 for condID = 1:length(cond)
     vRAll_dark = [];
     vRAll_CL = [];
     vFAll_dark = [];
     vFAll_CL = [];
     meanDFAll_dark = [];
     meanDFAll_CL = [];
     maxDFAll_dark = [];
     maxDFAll_CL = [];
     PVAStrenAll_dark = [];
     PVAStrenAll_CL = [];
     meanDFOffPk_dark = [];
     meanDFOffPk_CL = [];
     
    for flyID = 1:cond{condID}.numFlies
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
            % Get the data for this trial
            datNow = cond{condID}.allFlyData{flyID}.(trialName){trialID};

            % Get the behavior
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            tPts = tPts - tPts(1);
            heading = pi/180*datNow.positionDatMatch.PosRotMatch;
            stripePos = datNow.positionDatMatch.OffsetRotMatch(:,2);
            vR = datNow.positionDatMatch.vRot;
            vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);
            vF = datNow.positionDatMatch.vF;
            vF = sgolayfilt(vF,sgolayOrder,sgolayFrames);

            % Get the calcium activity and calculate the relevant metrics
            DF = datNow.RROIaveMax;
            DFMax = max(DF);
            DFMean = mean(DF,1);
            DFMeanAwayFromPeak = [];
            for tPt = 1:length(tPts)
                pkCent = find(DF(:,tPt)==DFMax(tPt));
                pkMin = pkCent-1;
                if pkMin == 0
                    pkMin = 16;
                end
                pkMax = pkCent + 1;
                if pkMax == 16
                    pkMax = 1;
                end
                offPkRng = setdiff([1:16],[pkMin,pkCent,pkMax]);
                DFMeanAwayFromPeak(tPt) = mean(DF(offPkRng,tPt));
            end
            [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));
            angs = angs-0.5*mean(diff(angs));
            PVAPlt(find(PVAStren<PVAThresh)) = NaN;
            for tStep = 2:length(tPts)
                if abs(PVAPlt(tStep)-PVAPlt(tStep-1))> pi
                    PVAPlt(tStep-1) = NaN;
                end
            end

            % Find the different visual conditions over the trial
            [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
            darkJumps = find(diff(darkPer)>1);
            if ~isempty(darkJumps);
                darkPer(darkJumps(1)+1:end) = [];
            end

            % Find when the stripe jumps
            stripeJumps = find(diff(datNow.positionDatMatch.Trans)~=0)+1;
            stripeJumps(end) = [];
            
            % Add the metrics to the population arrays
            vRAll_dark = vertcat(vRAll_dark,vR(darkPer));
            vRAll_CL = vertcat(vRAll_CL,vR(CLPer(1):stripeJumps(1)));
            vFAll_dark = vertcat(vFAll_dark,vF(darkPer));
            vFAll_CL = vertcat(vFAll_CL,vF(CLPer(1):stripeJumps(1)));
            meanDFAll_dark = horzcat(meanDFAll_dark,DFMean(darkPer));
            meanDFAll_CL = horzcat(meanDFAll_CL,DFMean(CLPer(1):stripeJumps(1)));
            maxDFAll_dark = horzcat(maxDFAll_dark,DFMax(darkPer));
            maxDFAll_CL = horzcat(maxDFAll_CL,DFMax(CLPer(1):stripeJumps(1)));
            PVAStrenAll_dark = vertcat(PVAStrenAll_dark,PVAStren(darkPer));
            PVAStrenAll_CL = vertcat(PVAStrenAll_CL,PVAStren(CLPer(1):stripeJumps(1)));
            meanDFOffPk_dark = horzcat(meanDFOffPk_dark,DFMeanAwayFromPeak(darkPer));
            meanDFOffPk_CL = horzcat(meanDFOffPk_CL,DFMeanAwayFromPeak(CLPer(1):stripeJumps(1)));
        end
        
        vRAll_dark(find(isnan(PVAStrenAll_dark))) = [];
        meanDFAll_dark(find(isnan(PVAStrenAll_dark))) = [];
        maxDFAll_dark(find(isnan(PVAStrenAll_dark))) = [];
        PVAStrenAll_dark(find(isnan(PVAStrenAll_dark))) = [];
%         meanDFOffPk_dark(find(isnan(PVAStrenAll_dark))) = [];
        
        % Plot it!
        for it = 1:4
            subplot(6,20,20*(it-1)+flyID+5*(condID-1));
            hold on;
            vRDisc = discretize(abs(vRAll_CL),vREdges);
            Means = [];
            P25 = [];
            P75 = [];
            if it == 0
                dataNow = meanDFAll_CL;
            elseif it == 2
                dataNow = maxDFAll_CL;
            elseif it == 3
                dataNow = PVAStrenAll_CL;
            else
                dataNow = meanDFOffPk_CL;
            end
            for discBin = 1:length(vRCents)
                rngNow = find(vRDisc == discBin);
                Means = [Means mean(dataNow(rngNow))];
                P25 = [P25 prctile(dataNow(rngNow),25)];
                P75 = [P75 prctile(dataNow(rngNow),75)];
            end
            pltRng = find(Means>0);
            p1 = plot(vRCents(pltRng),Means(pltRng),'b');
            patch('XData',[vRCents(pltRng) fliplr(vRCents(pltRng))],...
                'YData',[P25(pltRng) fliplr(P75(pltRng))],...
                'FaceColor','b','FaceAlpha',0.25,'EdgeColor','none');
            xlim([0 pi]);
            xticks([0 pi/2 pi]);

            
            subplot(6,20,76+4*it+[1:2]);
            hold on;
            if condID == 1
                plot(vRCents(pltRng),Means(pltRng),'b');
            elseif condID == 2
                plot(vRCents(pltRng),Means(pltRng),'Color',[0.5 0 1]);
            elseif condID == 3
                plot(vRCents(pltRng),Means(pltRng),'Color',[0 0 0.8]);
            else
                plot(vRCents(pltRng),Means(pltRng),'Color',[0.5 0 0.8]);
            end
            if it == 3
                ylim([0 0.3]);
            elseif it == 4
                ylim([0.1 0.3]);
            end
            xlim([0 pi]);
            xticks([0 pi/2 pi]);
            
            subplot(6,20,20*(it-1)+flyID+5*(condID-1));
            
            vRDisc = discretize(abs(vRAll_dark),vREdges);
            Means = [];
            P25 = [];
            P75 = [];
            if it == 1
                dataNow = meanDFAll_dark;
                ylim([0.1 0.3]);
                title(strcat('fly #',num2str(flyID)));
            elseif it == 2
                dataNow = maxDFAll_dark;
                ylim([0.2 0.7]);
            elseif it == 3
                dataNow = PVAStrenAll_dark;
                ylim([0 0.3]);
            else
                dataNow = meanDFOffPk_dark;
                ylim([0.1 0.25]);
                xlabel('vR (rad/s)');
            end
            for discBin = 1:length(vRCents)
                rngNow = find(vRDisc == discBin);
                datMoment = dataNow(rngNow);
                datMoment(find(isnan(datMoment))) = [];
                Means = [Means mean(datMoment)];
                P25 = [P25 prctile(datMoment,25)];
                P75 = [P75 prctile(datMoment,75)];
            end
            p2 = plot(vRCents,Means,'k');
            patch('XData',[vRCents fliplr(vRCents)],...
                'YData',[P25 fliplr(P75)],...
                'FaceColor','k','FaceAlpha',0.25,'EdgeColor','none');
            if flyID == 1 & condID == 1
                if it == 1
                    ylabel('mean DF');
                    legend([p1, p2],{'dark','CL'});
                    legend('boxoff');
                elseif it == 2
                    ylabel('max DF');
                elseif it == 3
                    ylabel('mean resultant vector length');
                else
                    ylabel('mean DF (off peak)');
                end
            end
            if flyID == 2 & it == 1
                text(pi,0.35,cond{condID}.name,'FontSize',12);
            end
            xlim([0 pi]);
            xticks([0 pi/2 pi]);
            
            subplot(6,20,78+4*it+[1:2]);
            hold on;
            if condID == 1
                plot(vRCents(pltRng),Means(pltRng),'k');
            elseif condID == 2
                plot(vRCents(pltRng),Means(pltRng),'Color',[0.5 0 0]);
            elseif condID == 3
                plot(vRCents(pltRng),Means(pltRng),'Color',[0 0.25 0.25]);
            else
                plot(vRCents(pltRng),Means(pltRng),'Color',[0.5 0.25 0.25]);
            end
            if it == 3
                ylim([0 0.3]);
            elseif it == 4
                ylim([0.1 0.3]);
            end
            xlim([0 pi]);
            xticks([0 pi/2 pi]);
        end
    end
 end
 
%% Create a summary figure for the paper - Figure 7 K, Figure S11 G, I

% Make the figure
paperFig = figure('units','normalized','outerposition',[0 0 1 1]);

% Specify which flies and trials to use as the example
flies = [4,4];
trials = [1,4];

% Set the PVA threshold
PVAThresh = 0.075;

% Set the savitzky-golay filtering parameters
sgolayOrder = 3;
sgolayFrames = 11;

% Plot the DF/F, the PVA, and the heading
for condID = 1:length(cond)
    for flyID = flies(condID)
        for trialID = trials(condID)
            % Get the data for this trial
            datNow = cond{condID}.allFlyData{flyID}.All{trialID};

            % Pull out the behavior
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            tPts = tPts - tPts(1);
            heading = pi/180*datNow.positionDatMatch.PosRotMatch;
            stripePos = 360/240*datNow.positionDatMatch.OffsetRotMatch(:,2);
            vR = datNow.positionDatMatch.vRot;
            vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);
            vF = datNow.positionDatMatch.vF;
            vF = sgolayfilt(vR,sgolayOrder,sgolayFrames);

            % Remove jumps for plotting
            headingPlt = heading;
            stripePosPlt = stripePos;
            for tStep = 2:length(tPts)
                if abs(headingPlt(tStep)-headingPlt(tStep-1))> pi
                    headingPlt(tStep-1) = NaN;
                end
                if abs(stripePosPlt(tStep)-stripePosPlt(tStep-1))> pi
                    stripePosPlt(tStep-1) = NaN;
                end
            end

            % Pull out the Ca activity
            DF = datNow.RROIaveMax;
            
            % Calculate activity properties and the PVA
            meanDF = mean(DF,1);
            maxDF = max(DF);
            [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));
            angs = angs-0.5*mean(diff(angs));
            PVAUW = UnWrap(PVAPlt,1.5,0);
            PVAPlt(find(PVAStren<PVAThresh)) = NaN;
            for tStep = 2:length(tPts)
                if abs(PVAPlt(tStep)-PVAPlt(tStep-1))> pi
                    PVAPlt(tStep-1) = NaN;
                end
            end

            % Find the different trial periods
            [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
            darkJumps = find(diff(darkPer)>1);
            if ~isempty(darkJumps);
                darkPer(darkJumps(1)+1:end) = [];
            end
            % Find when the stripe jumps
            stripeJumps = find(diff(datNow.positionDatMatch.Trans)~=0)+1;
            stripeJumps(end) = [];

            subplot(7,5,15*(condID-1)+[1:2]);
            hold on;
            imagesc(tPts,angs,DF);
            caxis([0 0.6]);
            ylim([-pi pi]);
            yticks([-pi 0 pi]);
            ylabel('EB pos. (rad)');
            xlim([tPts(1) tPts(stripeJumps(1))]);
%             xlim([tPts(1) tPts(darkPer(end))]);
            xlims = get(gca,'XLim');
            text(xlims(1),1.5*pi,cond{condID}.name,'FontSize',14);
            plot([tPts(darkPer(end)) tPts(darkPer(end))],[-pi pi],'--k');
            
            
            PVAStrenCol = zeros(length(PVAStren),3);
            for tPt = 1:length(PVAStren)
                colNow = min(1,2.5*PVAStren(tPt));
                PVAStrenCol(tPt,:) =  [colNow colNow colNow];
            end
            
            subplot(7,5,15*(condID-1)+[6:7]);
            hold on;
            plot(tPts(darkPer),headingPlt(darkPer),'k');
            plot(tPts(CLPer),stripePosPlt(CLPer),'b');
            plot(tPts,PVAPlt,'Color',[0.25 0 0.5]);
            ylim([-pi pi]);
            ylabel('pos. (rad)');
                        xlim([tPts(1) tPts(stripeJumps(1))]);
%             xlim([tPts(1) tPts(darkPer(end))]);
            plot([tPts(darkPer(end)) tPts(darkPer(end))],[-pi pi],'--k');
            if condID == 1
                legend({'heading (dark)','stripe position (CL)'});
                legend('boxoff');
            end
            xlabel('time (s)');
            
%             subplot(7,5,15*(condID-1)+[11:13]);
%             hold on;
% %             tRng = 0.5*mean(diff(tPts));
% %             for tPt = 1:length(PVAStren)
% %                 rectangle('Position',[tPts(tPt)-tRng -0.1 tRng 0.2],...
% %                     'FaceColor',PVAStrenCol(tPt,:),'EdgeColor','none');
% %             end
% %             axis off;
% %             text(tPts(1),0.2,'PVA strength');
% %             plot(tPts,meanDF);
% %             plot(tPts,maxDF);
%             plot(tPts,meanDF);
% %             plot(tPts(2:end),abs(vR)/(2*pi));
%             plot([tPts(1) tPts(stripeJumps(1))],[0.2 0.2],'--k');
%             xlim([tPts(1) tPts(stripeJumps(1))]);
%             xlabel('time (s)');
%             ylim([0 0.4]);

        end
    end
    colormap(brewermap(64, 'Blues'));
end

%% Plot the general bump properties - Figure 1 J, Figure S1

paperFig = figure('units','normalized','outerposition',[0 0 1 1]);

% Rescale the heading to cover a full 360 deg.
headingScale = 360/240;

% Look at the control flies
condID = 1;
trialType = 'All';

% Specify the filtering parameters
filtParam1 = 3;
filtParam2 = 11;
PVAThresh = 0.075;

pltCol = 'k';

% Now plot example of the bump over time
flyID_ex = 5;
trialID_ex =2;

% Pull the data for this trial
datNow = cond{condID}.allFlyData{flyID_ex}.(trialType){trialID_ex};

% Pull the behavior
tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
tPts = tPts - tPts(1);
heading = pi/180*datNow.positionDatMatch.PosRotMatch;
stripePos = headingScale*datNow.positionDatMatch.OffsetRotMatch(:,2);

headingPlt = heading;
stripePosPlt = UnWrap(stripePos,2,0);

% Pull the Ca2+ activity
DF = sgolayfilt(datNow.RROIaveMax,3,11,[],2);
[angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));
angs = angs-0.5*mean(diff(angs));
PVAUW = UnWrap(PVAPlt,1.5,1);
PVAPlt(find(PVAStren<PVAThresh)) = NaN;
for tStep = 2:length(tPts)
    if abs(PVAPlt(tStep)-PVAPlt(tStep-1))> pi
        PVAPlt(tStep-1) = NaN;
    end
end

% Find the different trial conditions
[darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
darkJumps = find(diff(darkPer)>1);
if ~isempty(darkJumps);
    darkPer(darkJumps(1)+1:end) = [];
end
% Find when the stripe jumps
stripeJumps = find(diff(datNow.positionDatMatch.Trans)~=0)+1;
stripeJumps(end) = [];

% Plot the trial data
subplot(5,4,[1:3]);
hold on;
imagesc(tPts,angs,DF);
plot(tPts,PVAPlt,'Color',[0.25 0 0.5]);
caxis([0 0.5]);
ylim([-pi pi]);
yticks([-pi 0 pi]);
ylabel('EB pos. (rad)');
xlim([tPts(1) tPts(stripeJumps(1))]);
xlims = get(gca,'XLim');
text(xlims(1),1.5*pi,cond{condID}.name,'FontSize',14);
plot([tPts(darkPer(end)) tPts(darkPer(end))],[-pi pi],'--k');
colormap(brewermap(64, 'Blues'));

zeroPt = 4*pi;%PVAUW(darkPer(end))-stripePosPlt(darkPer(end));
PVAUW(find(PVAStren<PVAThresh)) = NaN;
subplot(5,4,[5:7]);
hold on;
plot(tPts(darkPer),stripePosPlt(darkPer),'k');
plot(tPts(CLPer),stripePosPlt(CLPer),'b');
plot(tPts,PVAUW-zeroPt,'Color',[0.25 0 0.5]);
xlabel('time (s)');
ylabel('pos. (rad)');
xlim([tPts(1) tPts(stripeJumps(1))]);
yticks(linspace(-20*pi,20*pi,21));
yl = get(gca,'ylim');
plot([tPts(darkPer(end)) tPts(darkPer(end))],[yl(1) yl(2)],'--k');
if condID == 1
    legend({'heading (dark)','stripe position (CL)'});
    legend('boxoff');
end

% Plot examples of the offset distribution
stripeMult = 360/240;
offsetBins = linspace(-pi,pi,17);
offsetCents = offsetBins;
offsetCents(end) = [];
offsetCents = offsetCents + 0.5*mean(diff(offsetCents));

datNow = cond{condID}.allFlyData{flyID_ex}.(trialType){trialID_ex};

tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
tPts = tPts - tPts(1);
heading = pi/180*datNow.positionDatMatch.PosRotMatch;
stripePos = datNow.positionDatMatch.OffsetRotMatch(:,2);

DF = datNow.RROIaveMax;
[angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));

[darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
darkJumps = find(diff(darkPer)>1);
if ~isempty(darkJumps)
    darkPer(darkJumps(1)+1:end) = [];
end

% Find when the stripe jumps
stripeJumps = find(diff(datNow.positionDatMatch.Trans)~=0)+1;
stripeJumps(end) = [];

% Get the closed loop period before the stripe jumps
CLCompPer = [CLPer(1):stripeJumps(1)-1];

% Find the offset distribution around these two periods
darkRng = darkPer(PVAStren(darkPer)>PVAThresh);
darkOffset = stripeMult*stripePos(darkPer)-PVAPlt(darkPer);
darkOffset(darkOffset>pi) = darkOffset(darkOffset>pi)-2*pi;
darkOffset(darkOffset<-pi) = darkOffset(darkOffset<-pi)+2*pi;

CLRng = CLCompPer(PVAStren(CLCompPer)>PVAThresh);
CLOffset = stripeMult*stripePos(CLPer)-PVAPlt(CLPer);
CLOffset(CLOffset>pi) = CLOffset(CLOffset>pi)-2*pi;
CLOffset(CLOffset<-pi) = CLOffset(CLOffset<-pi)+2*pi;

% Plot the distributions
subplot(5,4,9);
hold on;
h1 = histogram(darkOffset,offsetBins);
h1.FaceColor = 'k';
h1.FaceAlpha = 0.1;
xlim([-pi pi]);
xticks([-pi 0 pi]);
ylim([0 600]);
box off;
if condID == 1
    ylabel('counts');
end
xlabel('offset (rad)');

subplot(5,4,10);
h2 = histogram(CLOffset,offsetBins);
h2.FaceColor = 'b';
h2.FaceAlpha = 0.1;
xlim([-pi pi]);
xticks([-pi 0 pi]);
ylim([0 600]);
box off;
xlabel('offset (rad)');

% Now plot the quantification

% Look at the offset distribution in the dark and in the initial closed loop period
trialName = 'All';
PVAThresh = 0.075;
stripeMult = 1;
offsetBins = linspace(-pi,pi,17);
offsetCents = offsetBins;
offsetCents(end) = [];
offsetCents = offsetCents + 0.5*mean(diff(offsetCents));

% Find the circular variance
for flyID = 1:cond{condID}.numFlies
    for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
        datNow = cond{condID}.allFlyData{flyID}.(trialName){trialID};

        tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
        tPts = tPts - tPts(1);
        heading = pi/180*datNow.positionDatMatch.PosRotMatch;
        stripePos = datNow.positionDatMatch.OffsetRotMatch(:,2);

        DF = datNow.RROIaveMax-1;
        [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));

        [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
        darkJumps = find(diff(darkPer)>1);
        if ~isempty(darkJumps)
            darkPer(darkJumps(1)+1:end) = [];
        end

        % Get the closed loop period before the stripe jumps
        CLCompPer = CLPer;

        % Find the offset distribution around these two periods
        darkRng = darkPer(PVAStren(darkPer)>PVAThresh);
        darkOffset = stripeMult*stripePos(darkPer)-PVAPlt(darkPer);
        darkOffset(darkOffset>pi) = darkOffset(darkOffset>pi)-2*pi;
        darkOffset(darkOffset<-pi) = darkOffset(darkOffset<-pi)+2*pi;

        CLRng = CLCompPer(PVAStren(CLCompPer)>PVAThresh);
        CLOffset = stripeMult*stripePos(CLPer)-PVAPlt(CLPer);
        CLOffset(CLOffset>pi) = CLOffset(CLOffset>pi)-2*pi;
        CLOffset(CLOffset<-pi) = CLOffset(CLOffset<-pi)+2*pi;

        [Ndark,edges] = histcounts(darkOffset,offsetBins);           
        [NCL,edges] = histcounts(CLOffset,offsetBins);           

        [S s] = circ_var(offsetCents',Ndark');
        cond{condID}.allFlyData{flyID}.(trialName){trialID}.darkPVA = S;
        [S s] = circ_var(offsetCents',NCL');
        cond{condID}.allFlyData{flyID}.(trialName){trialID}.CLPVA = S;
    end
end     

% Plot the circular variance
subplot(5,4,[4 8]);
hold on;
for flyID = 1:cond{condID}.numFlies
    meanDark = 0;
    meanCL = 0;
    for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)
        meanDark = meanDark + cond{condID}.allFlyData{flyID}.(trialName){trialID}.darkPVA;
        meanCL = meanCL + cond{condID}.allFlyData{flyID}.(trialName){trialID}.CLPVA;
    end
    meanDark = meanDark./length(cond{condID}.allFlyData{flyID}.All);
    meanCL = meanCL./length(cond{condID}.allFlyData{flyID}.All);

    scatter(0,meanDark,40,'k','filled');
    scatter(0.5,meanCL,40,'b','filled');
    alpha(0.5);
    plot([0 0.5],[meanDark meanCL],'k');
end
xlim([-0.5 1]);
ylim([0 1]);
ylabel('circular variance');
xticks([0 0.5 1 1.5]);
xticklabels({'dark','stripe','dark','stripe'});



% Plot the bump statistics - bump width, bump amplitude, bump stability - ctrl at RT only
minT = 5;

vFThresh = 0.1;
vRThresh = pi/10;

RBins = linspace(-1,1,11);
slopeBins = linspace(-4,4,21);

% Calculate the median slopes and slope percentile range
for flyID = 1:cond{condID}.numFlies
    allPerc = zeros(2);
    allMedSlopes = zeros(2);

    for visID = 1:2
        allRs = [];
        allSlopes = [];
        if visID == 1
            visCond = 'Dark';
        else
            visCond = 'CL';
        end

        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialType))

            % Get the data for this trial
            datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};

            % Find the different conditions across the trial
            [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
            darkEnd = find(diff(darkPer)>1);
            if ~isempty(darkEnd)
                darkPer(darkEnd+1:end) = [];
            end

            % Get the behavior
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            heading = headingScale*datNow.positionDatMatch.OffsetRotMatch(:,2);

            % Get the Ca activity
            DF = sgolayfilt(datNow.RROIaveMax, filtParam1, filtParam2, [], 2);
            [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));

            per = {};
            if strcmp(visCond,'Dark')
                per{1} = darkPer;
            elseif strcmp(visCond,'CL')
                CLPer(find(tPts(CLPer)<tPts(CLPer(1))+30)) = [];
                per{1} = CLPer;
            end

            for perNow = 1:length(per)
                PVANow = UnWrap(PVAPlt(per{perNow}),1.5,0);
                PVAStrenNow = PVAStren(per{perNow});
                headingNow = UnWrap(heading(per{perNow}),2,0);
                tPtsNow = tPts(per{perNow});

                wind = round(minT/mean(diff(tPtsNow)));

                for tStep = 1:length(per{perNow})-wind+1
                    range = [tStep:tStep+wind-1];
                    if ~isempty(find(PVAStrenNow(range)<PVAThresh))
                        continue;
                    else
                        CorrGain= robustfit( headingNow(range), PVANow(range));
                        allSlopes = [allSlopes CorrGain(2)];

                        % Get correlation coefficient
                        [R,P] = corrcoef( headingNow(range), PVANow(range));
                        allRs = [allRs R(1,2)];
                        
                        if flyID == flyID_ex & visID == 2 & (mod(tStep,30) == 0)
                            subplot(5,4,[16 20]);
                            hold on;
                            scatter(headingNow(range)-headingNow(range(1)),...
                                PVANow(range)-PVANow(range(1)),20,pltCol,'filled');
                            alpha(0.2);
                        end
                    end
                end
            end
            % Plot the data for the eample fly
            if flyID == flyID_ex & trialID == trialID_ex & visID == 2
                subplot(5,4,[16 20]);
                hold on;
                plot([-1.5*pi 1.5*pi],[-1.5*pi 1.5*pi],'Color','k');
                plot([-1.5*pi 1.5*pi],[0 0],'Color','k');
                plot([0 0],[-1.5*pi 1.5*pi],'Color','k');
                yticks([-2*pi -pi 0 pi 2*pi]);
                xlim([-2*pi 2*pi]);
                xticks([-2*pi -pi 0 pi 2*pi]);
                xlabel('heading (rad)');
                ylabel('PVA (rad)');
                axis equal;

                subplot(5,4,12);
                hold on;
                histogram(allSlopes,slopeBins,'FaceColor',pltCol);
                yl = ylim;
                line([prctile(allSlopes,25) prctile(allSlopes,25)],[0 yl(2)],'Color','k','LineStyle','--');
                line([prctile(allSlopes,75) prctile(allSlopes,75)],[0 yl(2)],'Color','k','LineStyle','--');
                line([1 1],[0 yl(2)],'Color','k');
                xlim([-4 4]);
                xlabel('slope');
            end
        end

        % Plot the median and percentile data
        allSlopes = allSlopes(find(~isnan(allSlopes)));
        subplot(5,4,17);
        hold on;
        scatter(visID,prctile(allSlopes,75)-prctile(allSlopes,25),50,pltCol,'filled');
        alpha(0.5);
        ylim([0 1]);
        xlim([0.5 2.5]);
        ylabel('percentile diff');

        subplot(5,4,18);
        hold on;
        scatter(visID,median(allSlopes),50,pltCol,'filled');
        alpha(0.5);
        ylim([0 1.5]);
        xlim([0.5 2.5]);
        ylabel('median slope');

        allPerc(visID) = prctile(allSlopes,75)-prctile(allSlopes,25);
        allMedSlopes(visID) = median(allSlopes);
    end
    subplot(5,4,17);
    plot([1 2],allPerc,'k');
    xticks([1 2]);
    xticklabels({'dark','CL'});
    subplot(5,4,18);
    plot([1 2],allMedSlopes,'k');
    plot([0.5 2.5],[1 1],'k');
    xticks([1 2]);
    xticklabels({'dark','CL'});
end


% Get (and plot) the bmp width
allWids = zeros(2,2,10);

WidAmpDat = BumpAmpAndWidth2Cols(cond,'All',PVAThresh);

visCond = 'Dark';
for flyID = 1:cond{condID}.numFlies
    allWidsNow = WidAmpDat{condID}.(visCond){flyID}.allWid_R;
    allWids(condID,1,flyID) = mean(allWidsNow(~isnan(allWidsNow)));
end

for flyID = 1:cond{condID}.numFlies
    subplot(5,4,19)
    hold on;
    scatter(1,allWids(condID,1,flyID),50,'b','filled');
    yticks([0 pi/4 pi/2 3*pi/4 pi]);
    ylim([0 pi]);
    text(2*(condID-1)+1,-0.5,cond{condID}.name);
    if condID == 1
        ylabel('mean(FWHM)');
    end
    set(gca,'FontSize',10);
    alpha(0.4);
end

visCond = 'CL';
for flyID = 1:cond{condID}.numFlies
    allWidsNow = WidAmpDat{condID}.(visCond){flyID}.allWid_R;
    allWids(condID,2,flyID) = mean(allWidsNow(~isnan(allWidsNow)));
end

for flyID = 1:cond{condID}.numFlies
    subplot(5,4,19)
    hold on;
    scatter(2,allWids(condID,2,flyID),50,'m','filled');
    plot([1 2],allWids(condID,:,flyID),'Color','k');
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'dark','CL'});
    yticks([0 pi/4 pi/2 3*pi/4 pi]);
    ylim([0 pi]);
    if condID == 1
        ylabel('mean(FWHM)');
    end
    set(gca,'FontSize',10);
    alpha(0.4);
end