%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\Shi\';

%% Load the data
cond{1}.name = 'Empty Gal4';
cond{1}.dirs{1} = strcat(dataDir,'Delta7-Ctrl\20190701');
cond{1}.dirs{2} = strcat(dataDir,'Delta7-Ctrl\20190702');
cond{1}.dirs{3} = strcat(dataDir,'Delta7-Ctrl\20190715');
cond{1}.dirs{4} = strcat(dataDir,'Delta7-Ctrl\20191003');
cond{1}.dirs{5} = strcat(dataDir,'Delta7-Ctrl\20191007');
cond{1}.dirs{6} = strcat(dataDir,'Delta7-Ctrl\20191120');
cond{1}.dirs{7} = strcat(dataDir,'Delta7-Ctrl\20191203');

cond{2}.name = 'Delta 7s';
cond{2}.dirs{1} = strcat(dataDir,'Delta7\20190715');
cond{2}.dirs{2} = strcat(dataDir,'Delta7\20190716');
cond{2}.dirs{3} = strcat(dataDir,'Delta7\20191003');
cond{2}.dirs{4} = strcat(dataDir,'Delta7\20191119');
cond{2}.dirs{5} = strcat(dataDir,'Delta7\20191203');

cond = FlyDatLoad(1,cond);

%% Plot the behavior - Figure S12 E,F

% Specify the histogram bins
pHistEdges = linspace(-pi,pi,32);

% Set a threshold on movement
vRThresh = 0.1*pi;
vFThresh = 0.1;

% Specify colors for each fly
colCode = hsv(2*5);

% Set up a vector to hold the resultant vector lengths and angles
angles = cell(2,2);

for condID = 1:length(cond)
    for flyID = 1:cond{condID}.numFlies
        for hc = 1:2
            if hc == 1
                trialType = 'All';
                pltCol = 'b';
                trial1 = 1;
            else
                trialType = 'All_30C';
                pltCol = 'r';
                trial1 = 2;
            end
            for trialID = trial1:length(cond{condID}.allFlyData{flyID}.(trialType));
                
                % Get the data 
                datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};
                
                % Pull the behavioral info
                tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                tPts = tPts-tPts(1);
                heading = datNow.positionDatMatch.OffsetRotMatch(:,2);
                vR = datNow.positionDatMatch.vRot;
                vF = datNow.positionDatMatch.vF;
                
                % Sort the trial periods
                [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                darkSkip = darkPer(find(diff(darkPer)>1));
                darkPer(darkSkip(1)+1:end) = [];
                
                CLPer(find(tPts(CLPer)<tPts(CLPer(1))+30)) = [];
                
                % Find the fly's orientation while moving
                flyMov = intersect(find(vR(CLPer)>vRThresh),find(vF(CLPer)>vFThresh));
                [N,angEdge] = histcounts(heading(CLPer(flyMov)),pHistEdges);
                arrowAng = circ_mean((angEdge(1:end-1)+mean(diff(angEdge))/2)',N');
                arrowL = circ_r((angEdge(1:end-1)+mean(diff(angEdge))/2)',N');
                angles{condID,hc} = [angles{condID, hc},arrowAng];            
            end
        end
    end
end

% Plot the headings for the control and D7 shi flies at both temperatures
angleFig = figure('units','normalized','outerposition',[0 0 0.25 0.5]);

subplot(2,2,1);
hold on;
[counts,bins] = hist(angles{1,1},[-pi:pi/4:pi]); %# get counts and bin locations
bar(bins,counts,'b')
ylim([0 15]);
xticks([-pi:pi/2:pi]);
xlabel('mean angle');
ylabel('counts');

subplot(2,2,3);
hold on;
[counts,bins] = hist(angles{1,2},[-pi:pi/4:pi]); %# get counts and bin locations
bar(bins,counts,'r')
ylim([0 15]);
xticks([-pi:pi/2:pi]);
xlabel('mean angle');
ylabel('counts');

subplot(2,2,2);
hold on;
[counts,bins] = hist(angles{2,1},[-pi:pi/4:pi]); %# get counts and bin locations
bar(bins,counts,'b')
ylim([0 15]);
xticks([-pi:pi/2:pi]);
xlabel('mean angle');
ylabel('counts');

subplot(2,2,4);
hold on;
[counts,bins] = hist(angles{2,2},[-pi:pi/4:pi]); %# get counts and bin locations
bar(bins,counts,'r')
ylim([0 15]);
xticks([-pi:pi/2:pi]);
xlabel('mean angle');
ylabel('counts');

%% Make a summary figure - Figure 3 K, P-S, Figure S6 D

% Set a PVA threshold
PVAThresh = 0.075;

% Set a time window
minT = 5;

% Set savitzky-golay filter parameters 
filtParam1 = 3;
filtParam2 = 11;

% Set thresholds to determine if the fly is moving
vFThresh = 0.1;
vRThresh = pi/10;

% Set up the figure and choose example flies
SummaryFig = figure('units','normalized','outerposition',[0 0 1 1]);
condIDs = [1,2];
flyIDs = [5,1];
trialIDs_Cold = [2,2];
trialIDs_Hot = [1,3];
visCond = 'CL';

for condID = condIDs
    flyID = flyIDs(condID);    
    for hc = 1:2
        if hc == 1
            trialID = trialIDs_Cold(condID);
            trialType = 'All';
        else
            trialID = trialIDs_Hot(condID);
            trialType = 'All_30C';
        end
        
        % Get the data for this fly and this trial
        datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};
                
        % Pull the behavioral info
        tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
        tPts = tPts-tPts(1);
        heading = datNow.positionDatMatch.OffsetRotMatch(:,2);

        % Determine if the fly is moving
        vR = datNow.positionDatMatch.vRot;
        vF = datNow.positionDatMatch.vF;
        flyStop = intersect(find(vR<vRThresh),find(vF<vFThresh));
        
        % Determine the different trail periods
        [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
        darkSkip = darkPer(find(diff(darkPer)>1));
        darkPer(darkSkip(1)+1:end) = [];

        CLPer(find(tPts(CLPer)<tPts(CLPer(1))+30)) = [];
        
        if strcmp(visCond,'dark')
            plotRng = darkPer;
            pltCol = 'k';
        elseif strcmp(visCond,'CL');
            plotRng = CLPer;
            pltCol = 'b';
            tPts = tPts - tPts(CLPer(1));
        end
        
        % Filter the imaging data for plotting
        DF = sgolayfilt(datNow.ROIaveMax-1, filtParam1, filtParam2, [], 2);
  
        % Get the PVA
        [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));
        
        PVAPaused = PVAPlt;
        for tPt = 1:length(flyStop)
            if PVAStren(flyStop(tPt))<PVAThresh & flyStop(tPt)>1
                PVAPaused(flyStop(tPt)) = PVAPaused(flyStop(tPt)-1);
            end
        end
        
        % Unwrap the data
        PVAUnwrap = UnWrap(PVAPlt,2,0);
        headingUnwrap = UnWrap(heading,2,0);

        for tPt = 2:length(tPts) 
            if abs(PVAPlt(tPt) - PVAPlt(tPt-1)) > pi
                PVAPlt(tPt-1) = NaN;
            end
            if abs(heading(tPt) - heading(tPt-1)) > pi
                heading(tPt-1) = NaN;
            end
        end
        PVAPlt(find(PVAStren<PVAThresh)) = NaN;

        % Plot the trial examples
        subplot(6,6,3*condID-2+6*(hc-1)+[0:2]);
        hold on;
        imagesc(tPts(plotRng),angs,DF(:,plotRng));
        cbh = colorbar;
        set(cbh,'YTick',[0:0.25:0.75])
        plot(tPts(plotRng),heading(plotRng),'Color',pltCol);
        caxis([0 0.75]);
%         xlim([tPts(plotRng(1)) tPts(plotRng(end))]);
        xlim([0 60]);
        ylim([-pi pi])
        xlabel('time (s)');
        ylabel('position (rad)');
        yticks([-pi 0 pi]);
    end
end
colormap(brewermap(64, 'Blues'));     

% Look at the data within a sliding window

% First, chunk up the data by bouts where the bump can be tracked
RBins = linspace(-1,1,11);
slopeBins = linspace(-4,4,21);

for condID = 1:length(cond)  
    for flyID = 1:cond{condID}.numFlies
        allPerc = zeros(2,2);
        allMedSlopes = zeros(2,2);
        for hc = 1:2
            if hc == 1
                trialType = 'All';
                trialStart = 1;
                trialStop = 2;
                pltCol = 'b';
            else
                trialType = 'All_30C';
                trialStart = 2;
                trialStop = 3;
                pltCol = 'r';
            end
            
            for visID = 1:2
                allRs = [];
                allSlopes = [];
                if visID == 1
                    visCond = 'Dark';
                else
                    visCond = 'CL';
                end
                
                chunkedDat.Fly = flyID;
                chunkedDat.FlyPositionUnwrapNorm_All = {};
                chunkedDat.BumpPositionUnwrapNorm_All = {};
                boutCount = 0;

                for trialID = trialStart:trialStop

                    % Skip trials where the fly tracking failed
                    if trialID > length(cond{condID}.allFlyData{flyID}.(trialType))
                        continue;
                    end
                    
                    if visID == 1
                        if condID == 1 & hc == 1 & flyID == 6 & trialID == 2
                            continue;
                        end
                        if condID == 2 & hc == 2 & flyID == 7 & trialID == 1
                            continue;
                        end
                    end
                    datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};

                    % Pull out the different trial periods
                    [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                    darkEnd = find(diff(darkPer)>1);
                    if ~isempty(darkEnd)
                        darkPer(darkEnd+1:end) = [];
                    end

                    % Pull the behavior
                    tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                    heading = datNow.positionDatMatch.OffsetRotMatch(:,2);

                    % Pull the imaging data
                    DF = sgolayfilt(datNow.ROIaveMax, filtParam1, filtParam2, [], 2);
                    [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));

                    per = {};
                    if strcmp(visCond,'Dark')
                        per{1} = darkPer;
                    elseif strcmp(visCond,'CL')
                        CLPer(find(tPts(CLPer)<tPts(CLPer(1))+30)) = [];
                        per{1} = CLPer;
                    end

                    % Look at each window
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
                                
                                % Plot the example fly data
                                if flyID == flyIDs(condID) & visID == 2 & (mod(tStep,30) == 0)
                                    subplot(6,6,12+3*condID-2+[0 6]);
                                    hold on;
                                    scatter(headingNow(range)-headingNow(range(1)),...
                                        PVANow(range)-PVANow(range(1)),20,pltCol,'filled');
%                                     plot([min(headingNow(range))-headingNow(range(1)) max(headingNow(range))-headingNow(range(1))],...
%                                         [CorrGain(2)*min(headingNow(range))+CorrGain(1)-PVANow(range(1)) ...
%                                         CorrGain(2)*max(headingNow(range))+CorrGain(1)-PVANow(range(1))],...
%                                         'Color',pltCol);
                                    alpha(0.2);
                                end
                            end
                        end
                    end
                end
                % Set plotting params for the example fly data
                if flyID == flyIDs(condID) & visID == 2
                    subplot(6,6,12+3*condID-2+[0 6]);
                    plot([-1.5*pi 1.5*pi],[-1.5*pi 1.5*pi],'Color','k');
                    plot([-1.5*pi 1.5*pi],[0 0],'Color','k');
                    plot([0 0],[-1.5*pi 1.5*pi],'Color','k');
                    yticks([-2*pi -pi 0 pi 2*pi]);
                    xlim([-1.5*pi 1.5*pi]);
                    xticks([-1.5*pi -pi -0.5*pi 0 0.5*pi pi 1.5*pi]);
                    xlabel('heading (rad)');
                    ylabel('PVA (rad)');
                    axis equal;
                end
                    
                % Plot the histograms of the slopes for the example flies
                if flyID == flyIDs(condID) & visID == 2
%                     subplot(6,6,12+3*condID-2+6*(hc-1));
%                     hold on;
%                     histogram(allRs,RBins,'FaceColor',pltCol);
%                     yl = ylim;
%                     line([median(allRs) median(allRs)],[0 yl(2)],'Color','k');
%                     title({['Fly ' num2str(flyID) '-' visCond]});
%                     xlim([-1 1]);
%                     if visID == 2 & hc == 2
%                         xlabel('r');
%                     end
%                     if flyID == 1 & visID == 1 & hc == 1
%                         text(-1.25, 1.25*yl(2),cond{condID}.name);
%                     end

                    subplot(6,6,14+3*(condID-1)+6*(hc-1));
                    hold on;
                    histogram(allSlopes,slopeBins,'FaceColor',pltCol);
                    yl = ylim;
%                     line([median(allSlopes) median(allSlopes)],[0 yl(2)],'Color','k');
                    line([prctile(allSlopes,25) prctile(allSlopes,25)],[0 yl(2)],'Color','k','LineStyle','--');
                    line([prctile(allSlopes,75) prctile(allSlopes,75)],[0 yl(2)],'Color','k','LineStyle','--');
                    line([1 1],[0 yl(2)],'Color','k');
                    title({['Fly ' num2str(flyID) '-' visCond]});
                    xlim([-4 4]);
                    if visID == 2 & hc == 2
                        xlabel('slope');
                    end
                end

%                 subplot(6,6,25+3*(condID-1)+6*(visID-1));
%                 hold on;
%                 scatter(flyID+10*(hc-1),median(allRs),20,pltCol,'filled');
%                 title(['R for ' visCond '-' cond{condID}.name]);
%                 ylim([-1 1]);
% 
%                 subplot(6,6,26+3*(condID-1)+6*(visID-1));
%                 hold on;
%                 scatter(flyID+10*(hc-1),median(allSlopes),20,pltCol,'filled');
%                 title(['slope for ' visCond '-' cond{condID}.name]);
%                 ylim([-0.5 1.5]);
%                 plot([1 18],[1 1],'k');
                
    
                allSlopes = allSlopes(find(~isnan(allSlopes)));
                subplot(6,6,25+6*(visID-1));
                hold on;
                scatter(hc+2*(condID-1),prctile(allSlopes,75)-prctile(allSlopes,25),50,pltCol,'filled');
%                 if flyID == flyIDs(condID)
%                     scatter(hc+2*(condID-1),prctile(allSlopes,75)-prctile(allSlopes,25),50,'green','filled');
%                 end
                alpha(0.5);
                ylim([-1 12]);
                xlim([0.5 4.5]);
                ylabel('percentile diff');
                
                subplot(6,6,28+6*(visID-1));
                hold on;
                scatter(hc+2*(condID-1),median(allSlopes),50,pltCol,'filled');
%                 if flyID == flyIDs(condID)
%                     scatter(hc+2*(condID-1),median(allSlopes),50,'green','filled');
%                 end
                alpha(0.5);
                ylim([-1 5]);
                xlim([0.5 4.5]);
                ylabel('median slope');
                plot([0.5 4.5],[1 1],'k');
                
                allPerc(hc,visID) = prctile(allSlopes,75)-prctile(allSlopes,25);
                allMedSlopes(hc,visID) = median(allSlopes);
            end
        end
        subplot(6,6,25);
        plot([1 2]+2*(condID-1),allPerc(:,1),'k');
        subplot(6,6,31);
        plot([1 2]+2*(condID-1),allPerc(:,2),'k');
        
        subplot(6,6,28);
        plot([1 2]+2*(condID-1),allMedSlopes(:,1),'k');
        subplot(6,6,34);
        plot([1 2]+2*(condID-1),allMedSlopes(:,2),'k');
        
        if flyID == 1
            for pltNow = [1:4]
                if pltNow < 3
                    visCond = 'dark';
                else
                    visCond = 'CL';
                end
                subplot(6,6,25+3*(pltNow-1));
                text(1+2*(condID-1),-3,[visCond '-' cond{condID}.name]);
                xticks([1 2 3 4]);
                xticklabels({'RT','30C','RT','30C'});
            end
        end
    end
end


% Plot the behavior
pHistEdges = linspace(-pi,pi,32);

% Set a threshold on movement
vRThresh = 0.1*pi;
vFThresh = 0.1;

for condID = 1:length(cond)
    for flyID = 1:cond{condID}.numFlies
        for hc = 1:2
            if hc == 1
                trialType = 'All';
                trialStart = 1;
                trialStop = 2;
                pltCol = 'b';
            else
                trialType = 'All_30C';
                trialStart = 2;
                trialStop = 3;
                pltCol = 'r';
            end
            for trialID = trialStart:trialStop
                
                if trialID > length(cond{condID}.allFlyData{flyID}.(trialType))
                    continue;
                end
                if visID == 1
                    if condID == 1 & hc == 1 & flyID == 6 & trialID == 2
                        continue;
                    end
                    if condID == 2 & hc == 2 & flyID == 7 & trialID == 1
                        continue;
                    end
                end
                datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};
                
                tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                tPts = tPts-tPts(1);
                heading = datNow.positionDatMatch.OffsetRotMatch(:,2);
                vR = datNow.positionDatMatch.vRot;
                vF = datNow.positionDatMatch.vF;
                
                [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                darkSkip = darkPer(find(diff(darkPer)>1));
                darkPer(darkSkip(1)+1:end) = [];
                
                CLPer(find(tPts(CLPer)<tPts(CLPer(1))+30)) = [];
                
                flyMov = intersect(find(vR(CLPer)>vRThresh),find(vF(CLPer)>vFThresh));
                [N,angEdge] = histcounts(heading(CLPer(flyMov)),pHistEdges);
                arrowAng = circ_mean((angEdge(1:end-1)+mean(diff(angEdge))/2)',N');
                arrowL = circ_r((angEdge(1:end-1)+mean(diff(angEdge))/2)',N');
                
                if flyIDs(condID) == flyID
                    subplot(6,6,15+3*(condID-1));
                    hold on;
                    rectangle('Position',[-1 -1 2 2],'Curvature',1);
                    axis equal;
                    quiver(0,0,arrowL*sin(arrowAng),arrowL*cos(arrowAng),'Color',pltCol);
                    axis off;
                end
                
                subplot(6,6,21+3*(condID-1));
                hold on;
                rectangle('Position',[-1 -1 2 2],'Curvature',1);
                axis equal;
                quiver(0,0,arrowL*sin(arrowAng),arrowL*cos(arrowAng),'Color',pltCol);
                axis off;
                
            end
        end
    end
end

%% Look at R2 fits for various window sizes

% Set a PVA threshold.
PVAThresh = 0.075;

% Specify a range of times to look at
minTs = [1:15];

% Set the savitzky-golay filter parameters
filtParam1 = 3;
filtParam2 = 11;

% Bin the slopes and fits
RBins = linspace(-1,1,11);
slopeBins = linspace(-4,4,21);

StatPlt = figure('units','normalized','outerposition',[0 0 1 1]);
allRsAcrossFlies = zeros(length(cond),5,2,2,length(minTs));
for minT = minTs
    for condID = 1:length(cond)  
        for flyID = 1:cond{condID}.numFlies
            for hc = 1:2
                if hc == 1
                    trialType = 'All';
                    trialStart = 1;
                    trialStop = 2;
                    pltCol = 'b';
                else
                    trialType = 'All_30C';
                    trialStart = 2;
                    trialStop = 3;
                    pltCol = 'r';
                end

                for visID = 1:2
                    allRs = [];
                    allSlopes = [];
                    if visID == 1
                        visCond = 'Dark';
                    else
                        visCond = 'CL';
                        % Skip trials where the fly tracking failed
                        if condID == 1 & hc == 1 & flyID == 6 & trialID == 2
                            continue;
                        end
                        if condID == 2 & hc == 2 & flyID == 7 & trialID == 1
                            continue;
                        end
                    end
                    chunkedDat.Fly = flyID;
                    chunkedDat.FlyPositionUnwrapNorm_All = {};
                    chunkedDat.BumpPositionUnwrapNorm_All = {};
                    boutCount = 0;

                    for trialID = trialStart:trialStop

                        if trialID >= length(cond{condID}.allFlyData{flyID}.(trialType))
                            continue;
                        end
                            
                        datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};

                        [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                        darkEnd = find(diff(darkPer)>1);
                        if ~isempty(darkEnd)
                            darkPer(darkEnd+1:end) = [];
                        end

                        tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                        heading = datNow.positionDatMatch.OffsetRotMatch(:,2);

                        DF = sgolayfilt(datNow.ROIaveMax, filtParam1, filtParam2, [], 2);
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
                                    % Get correlation coefficient
                                    [R,P] = corrcoef( headingNow(range), PVANow(range));
                                    allRs = [allRs R(1,2)];
                                end
                            end
                        end
                    end
                    figure(StatPlt);
                    subplot(2,4,visID+2*(condID-1)+4*(hc-1));
                    hold on;
                    scatter(minT,median(allRs(find(~isnan(allRs)))),20,pltCol,'filled');
                    title(['median R for ' visCond '-' cond{condID}.name]);
                    
                    allRsAcrossFlies(condID,flyID,hc,visID,minT) = median(allRs(find(~isnan(allRs))));
                end
            end
        end
    end
end

for condID = 1:length(cond)  
    for flyID = 1:cond{condID}.numFlies
        for hc = 1:2
            if hc == 1
                pltCol = 'b';
            else
                pltCol = 'r';
            end
            for visID = 1:2
                subplot(2,4,visID+2*(condID-1)+4*(hc-1));
                hold on;
                plot(squeeze(allRsAcrossFlies(condID,flyID,hc,visID,:)),'Color',pltCol)
                xlabel('window length (t)');
                ylim([-1 1]);
                if flyID == 1
                    plot([min(minTs) max(minTs)],[0 0],'Color','k');
                end
            end
        end
    end
end

for sp = 1:8
    subplot(2,4,sp);
    line([4 4],[-1 1],'Color','k');
end

%% Look at stats for different sliding window lengths

% Set PVA threshold and savitzky-golay filter parameters
PVAThresh = 0.075;
filtParam1 = 3;
filtParam2 = 11;

RBins = linspace(-1,1,11);
slopeBins = linspace(-4,4,21);

for minT = [3 5 7]
    StatPlt = figure('units','normalized','outerposition',[0 0 1 1]);
    for condID = 1:length(cond)  
        for flyID = 1:cond{condID}.numFlies
            for hc = 1:2
                if hc == 1
                    trialType = 'All';
                    trialStart = 1;
                    trialStop = 2;
                    pltCol = 'b';
                else
                    trialType = 'All_30C';
                    trialStart = 2;
                    trialStop = 3;
                    pltCol = 'r';
                end

                for visID = 1:2
                    allRs = [];
                    allSlopes = [];
                    if visID == 1
                        visCond = 'Dark';
                    else
                        visCond = 'CL';
                    end
                    chunkedDat.Fly = flyID;
                    chunkedDat.FlyPositionUnwrapNorm_All = {};
                    chunkedDat.BumpPositionUnwrapNorm_All = {};
                    boutCount = 0;

                    for trialID = trialStart:trialStop
                        % Skip trials where the fly tracking failed
                        if condID == 1 & hc == 1 & flyID == 6 & trialID == 2
                            continue;
                        end
                        if condID == 2 & hc == 2 & flyID == 7 & trialID == 1
                            continue;
                        end

                        datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};

                        [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                        darkEnd = find(diff(darkPer)>1);
                        if ~isempty(darkEnd)
                            darkPer(darkEnd+1:end) = [];
                        end

                        tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                        heading = datNow.positionDatMatch.OffsetRotMatch(:,2);

                        DF = sgolayfilt(datNow.ROIaveMax, filtParam1, filtParam2, [], 2);
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
                                end
                            end
                        end
                    end

                    figure(StatPlt);
                    subplot(3,4,visID+2*(condID-1));
                    hold on;
                    scatter(flyID+10*(hc-1),median(allRs(find(~isnan(allRs)))),20,pltCol,'filled');
                    title(['R for ' visCond '-' cond{condID}.name]);
                    ylim([0 1]);
                    if flyID == 1 & condID == 1 & hc ==1 & visID == 1
                        text(1,1.1,strcat('window = ',num2str(minT),' s'),'FontSize',14);
                    end
                    xlim([0 18]);

                    subplot(3,4,4+visID+2*(condID-1));
                    hold on;
                    scatter(flyID+10*(hc-1),median(allSlopes(find(~isnan(allSlopes)))),20,pltCol,'filled');
                    title(['slope for ' visCond '-' cond{condID}.name]);
                    ylim([0 4.5]);
                    plot([1 18],[1 1],'k');
                    xlim([0 18]);

                    subplot(3,4,8+visID+2*(condID-1));
                    hold on;
                    scatter(flyID+10*(hc-1),prctile(allSlopes(find(~isnan(allSlopes))),75)-prctile(allSlopes(find(~isnan(allSlopes))),25),20,pltCol,'filled');
                    title(['percentile diff for ' visCond '-' cond{condID}.name]);
                    ylim([0 7]);
                    xlim([0 18]);
                end
            end
        end
    end
end

%% Plot the bump width and amplitude as a function of vR

% Set the PVA threshold
PVAThresh = 0.075;

% Specify the bins for the rotational velocity
vRBins = linspace(0,1.5*pi,13);
vRCents = vRBins;
vRCents(end) = [];
vRCents = vRCents + 0.5*mean(diff(vRCents));

% Set the colormaps for the hot and cold periods
blues = colormap(brewermap(10, '*Blues'));
reds = colormap(brewermap(10, '*Reds'));

% Specify the visual condition to consider
visCond = 'Dark';

AmpsAndWids = figure('units','normalized','outerposition',[0 0 1 1]);
for hc = 1:2
    if hc == 1
        trialName = 'All';
        flyCols = blues;
    else
        trialName = 'All_30C';
        flyCols = reds;
    end
    % Get the bump amplitude and width
    WidAmpDat = BumpAmpAndWidth(cond,trialName,PVAThresh);
    
    % Bin the amp and wids
    for condID = 1:length(cond)
        
        for flyID = 1:cond{condID}.numFlies
            
            allAmps = [];
            allWids = [];
            
            vRIDs = discretize(abs(WidAmpDat{condID}.(visCond){flyID}.allvR),vRBins);
            
            for ID = 1:max(vRIDs)
                pltRng = find(vRIDs == ID);
                
                allAmps = vertcat(allAmps, mean(WidAmpDat{condID}.(visCond){flyID}.allAmp(pltRng))-1);
                wids2plt = WidAmpDat{condID}.(visCond){flyID}.allWid(pltRng);
                wids2plt = wids2plt(~isnan(wids2plt));
                allWids = vertcat(allWids,mean(wids2plt));
            end
            
            subplot(4,4,condID+4*(hc-1))
            hold on;
            plot(vRCents(1:ID),allAmps,'Color',flyCols(flyID,:));
            xlim([0 pi]);
            xticks([0 pi/4 pi/2 3*pi/4 pi]);
            ylim([0 1.75]);
            title(strcat(cond{condID}.name,'-',strrep(trialName,'_','-')));
            if condID == 1
                ylabel('mean(max DF/F)');
            end
            if hc == 2
                xlabel('vR (rad/s)');
            end
            set(gca,'FontSize',10);

            subplot(4,4,2+condID+4*(hc-1))
            hold on;
            plot(vRCents(1:ID),allWids,'Color',flyCols(flyID,:));
            xlim([0 pi]);
            xticks([0 pi/4 pi/2 3*pi/4 pi]);
            yticks([0 pi/4 pi/2 3*pi/4 pi 1.25*pi]);
            ylim([0 1.25*pi]);
            title(strcat(cond{condID}.name,'-',strrep(trialName,'_','-')));
            if condID == 1
                ylabel('mean(FWHM)');
            end
            if hc == 2
                xlabel('vR (rad/s)');
            end
            set(gca,'FontSize',10);
        end
    end
end

%% Plot the bump width and amplitude as a function of vR - means - Figure S6 E,F

% Set the PVA threshold
PVAThresh = 0.075;

% Choose the visual condition
visCond = 'Dark';

AmpsAndWids = figure('units','normalized','outerposition',[0 0 1 1]);

% Create arrays to hold the mean values
allAmps = zeros(2,2,10);
allWids = zeros(2,2,10);

for hc = 1:2
    if hc == 1
        trialName = 'All';
    else
        trialName = 'All_30C';
    end
    WidAmpDat = BumpAmpAndWidth(cond,trialName,PVAThresh);
    
    for condID = 1:length(cond)
        for flyID = 1:cond{condID}.numFlies
            allAmps(condID,hc,flyID) = mean(WidAmpDat{condID}.(visCond){flyID}.allAmp-1);
            allWidsNow = WidAmpDat{condID}.(visCond){flyID}.allWid;
            allWids(condID,hc,flyID) = mean(allWidsNow(~isnan(allWidsNow)));
        end
    end
end

for condID = 1:length(cond)
    for flyID = 1:cond{condID}.numFlies
        subplot(2,3,1)
        hold on;
        plot([1,2]+2*(condID-1),[allAmps(condID,1,flyID),allAmps(condID,2,flyID)],'Color',[0.5 0.5 0.5]);
        scatter(1+2*(condID-1),allAmps(condID,1,flyID),100,'b','filled');
        scatter(2+2*(condID-1),allAmps(condID,2,flyID),100,'r','filled');
        xlim([0.5 4.5]);
        xticks([1 2 3 4]);
        xticklabels({'RT','30C','RT','30C'});
        ylim([0 1.75]);
        text(2*(condID-1)+1,-0.25,cond{condID}.name);
        if condID == 1
            ylabel('mean(max DF/F)');
        end
        set(gca,'FontSize',10);
        alpha(0.4);

        subplot(2,3,2)
        hold on;
        plot([1,2]+2*(condID-1),[allWids(condID,1,flyID),allWids(condID,2,flyID)],'Color',[0.5 0.5 0.5]);
        scatter(1+2*(condID-1),allWids(condID,1,flyID),100,'b','filled');
        scatter(2+2*(condID-1),allWids(condID,2,flyID),100,'r','filled');
        xlim([0.5 4.5]);
        xticks([1 2 3 4]);
        xticklabels({'RT','30C','RT','30C'});
        yticks([0 pi/4 pi/2 3*pi/4 pi]);
        ylim([0 pi]);
        text(2*(condID-1)+1,-0.5,cond{condID}.name);
        if condID == 1
            ylabel('mean(FWHM)');
        end
        set(gca,'FontSize',10);
        alpha(0.4);
    end
end

% Calculate the statistics
[h,p] = ttest(squeeze(allAmps(1,1,:)),squeeze(allAmps(1,2,:)))
[h,p] = ttest(squeeze(allAmps(2,1,:)),squeeze(allAmps(2,2,:)))
[h,p] = ttest(squeeze(allWids(1,1,:)),squeeze(allWids(1,2,:)))
[h,p] = ttest(squeeze(allWids(2,1,:)),squeeze(allWids(2,2,:)))

%% Run stats for the median slope and slopes percentile range

% Set the PVA threshold
PVAThresh = 0.075;

% Set the time window size 
minT = 5;

% Set the Savitzky-Golay window parameters
filtParam1 = 3;
filtParam2 = 11;

% Set movement thresholds
vFThresh = 0.1;
vRThresh = pi/10;

SummaryFig = figure('units','normalized','outerposition',[0 0 1 1]); 

% Bin the values
RBins = linspace(-1,1,11);
slopeBins = linspace(-4,4,21);
allPerc = zeros(2,7,2,2);
allMedSlopes = zeros(2,7,2,2);
for condID = 1:length(cond)  
    for flyID = 1:cond{condID}.numFlies
        for hc = 1:2
            if hc == 1
                trialType = 'All';
                trialStart = 1;
                trialStop = 2;
                pltCol = 'b';
            else
                trialType = 'All_30C';
                trialStart = 2;
                trialStop = 3;
                pltCol = 'r';
            end
            
            for visID = 1:2
                allRs = [];
                allSlopes = [];
                if visID == 1
                    visCond = 'Dark';
                else
                    visCond = 'CL';
                end
                
                chunkedDat.Fly = flyID;
                chunkedDat.FlyPositionUnwrapNorm_All = {};
                chunkedDat.BumpPositionUnwrapNorm_All = {};
                boutCount = 0;

                for trialID = trialStart:trialStop

                    % Skip trials where the fly tracking failed
                    if trialID > length(cond{condID}.allFlyData{flyID}.(trialType))
                        continue;
                    end
                    if visID == 1
                        if condID == 1 & hc == 1 & flyID == 6 & trialID == 2
                            continue;
                        end
                        if condID == 2 & hc == 2 & flyID == 7 & trialID == 1
                            continue;
                        end
                    end
                    datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};

                    [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                    darkEnd = find(diff(darkPer)>1);
                    if ~isempty(darkEnd)
                        darkPer(darkEnd+1:end) = [];
                    end

                    tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                    heading = datNow.positionDatMatch.OffsetRotMatch(:,2);

                    DF = sgolayfilt(datNow.ROIaveMax, filtParam1, filtParam2, [], 2);
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
                            end
                        end
                    end
                end
                
                allPerc(condID,flyID,hc,visID) = prctile(allSlopes,75)-prctile(allSlopes,25);
                allMedSlopes(condID,flyID,hc,visID) = median(allSlopes);
            end
        end
    end
end

% Plot the data
for condID = 1:length(cond)
    for flyID = 1:cond{condID}.numFlies
        for visID = 1:2
        
            subplot(2,4,4*visID-3);
            hold on;
            scatter([1:2]+2*(condID-1),allPerc(condID,flyID,:,visID),50,'k','filled');
            plot([1 2]+2*(condID-1),squeeze(allPerc(condID,flyID,:,visID)),'k');
            alpha(0.5);
            ylim([-1 9]);
            xlim([0.5 4.5]);
            ylabel('percentile diff');

            subplot(2,4,4*visID-2);
            hold on;
            scatter([1:2]+2*(condID-1),allMedSlopes(condID,flyID,:,visID),50,'k','filled');
            plot([1 2]+2*(condID-1),squeeze(allMedSlopes(condID,flyID,:,visID)),'k');
            alpha(0.5);
            ylim([-1 5]);
            xlim([0.5 4.5]);
            ylabel('median slope');
            plot([0.5 4.5],[1 1],'k');

            if flyID == 1
                for pltNow = [1:8]
                    if pltNow < 3
                        visCond = 'dark';
                    else
                        visCond = 'CL';
                    end
                    subplot(2,4,pltNow);
                    text(1+2*(condID-1),-3,[visCond '-' cond{condID}.name]);
                    xticks([1 2 3 4]);
                    xticklabels({'RT','30C','RT','30C'});
                end
            end
        end
    end
end

% Calculate the stats
[h,p] = ttest(squeeze(allPerc(1,:,1,1)),squeeze(allPerc(1,:,2,1)))
[h,p] = ttest(squeeze(allPerc(2,:,1,1)),squeeze(allPerc(2,:,2,1)))
[h,p] = ttest(squeeze(allPerc(1,:,1,2)),squeeze(allPerc(1,:,2,2)))
[h,p] = ttest(squeeze(allPerc(2,:,1,2)),squeeze(allPerc(2,:,2,2)))

[h,p] = ttest(squeeze(allMedSlopes(1,:,1,1)),squeeze(allMedSlopes(1,:,2,1)))
[h,p] = ttest(squeeze(allMedSlopes(2,:,1,1)),squeeze(allMedSlopes(2,:,2,1)))
[h,p] = ttest(squeeze(allMedSlopes(1,:,1,2)),squeeze(allMedSlopes(1,:,2,2)))
[h,p] = ttest(squeeze(allMedSlopes(2,:,1,2)),squeeze(allMedSlopes(2,:,2,2)))

%% Replot the window and slope example for the summary figure - Figure 3 N,O

% Specify the flys and conditions to look at
condIDs = [1,2];
flyIDs = [5,1];
visCond = 'CL';

% Select the number of point to sample
numPtsToSamp = 10;

% Set a PVA threshold
PVAThresh = 0.075;

% Select the window length
minT = 5;

% Set filtering parameters
filtParam1 = 3;
filtParam2 = 11;

% Create the figure
SummaryFig = figure('units','normalized','outerposition',[0 0 1 1]);

for condID = condIDs  
    for flyID = flyIDs(condID)
        for hc = 1:2
            if hc == 1
                trialType = 'All';
                trialStart = 1;
                trialStop = 2;
                pltCol = [0 0 1];
            else
                trialType = 'All_30C';
                trialStart = 2;
                trialStop = 3;
                pltCol = [1 0 0];
            end
            
            for visID = 1:2
                allRs = [];
                allSlopes = [];
                if visID == 1
                    visCond = 'Dark';
                else
                    visCond = 'CL';
                end
                
                subplot(2,4,2*(condID-1) + hc + [0 4]);
                hold on;
                plot([-1.5*pi 1.5*pi],[-1.5*pi 1.5*pi],'Color','k');
                plot([-1.5*pi 1.5*pi],[0 0],'Color','k');
                plot([0 0],[-1.5*pi 1.5*pi],'Color','k');
                axis equal;
                axis([-1.5*pi 1.5*pi -1.5*pi 1.5*pi]);
                xticks([-1.5*pi -pi -0.5*pi 0 0.5*pi pi 1.5*pi]);
                yticks([-1.5*pi -pi -0.5*pi 0 0.5*pi pi 1.5*pi]);
                xlabel('heading (rad)');
                ylabel('PVA (rad)');
                title(cond{condID}.name)
                
                for trialID = trialStart:trialStop

                    datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};

                    [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                    darkEnd = find(diff(darkPer)>1);
                    if ~isempty(darkEnd)
                        darkPer(darkEnd+1:end) = [];
                    end

                    tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                    heading = datNow.positionDatMatch.OffsetRotMatch(:,2);

                    DF = sgolayfilt(datNow.ROIaveMax, filtParam1, filtParam2, [], 2);
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
                        
                        sampTs = round((length(per{perNow})-wind+1)*rand(numPtsToSamp,1));
                        
                        for tStep = 1:length(per{perNow})-wind+1
                            range = [tStep:tStep+wind-1];
                            if ~isempty(find(PVAStrenNow(range)<PVAThresh))
                                continue;
                            else
                                CorrGain= robustfit( headingNow(range), PVANow(range));
                                
                                if flyID == flyIDs(condID) & visID == 2 & ismember(tStep,sampTs)
                                    
                                    xFit = linspace(...
                                        min(headingNow(range)),...
                                        max(headingNow(range)),...
                                        3);
                                    yFit = CorrGain(1)+CorrGain(2)*xFit;
                                    
                                    subplot(2,4,2*(condID-1) + hc + [0 4]);
                                    hold on;
                                    scatter(headingNow(range)-headingNow(range(1)),...
                                        PVANow(range)-PVANow(range(1)),5,pltCol,'filled');
                                    alpha(0.2);
                                    plot(xFit-headingNow(range(1)),...
                                        yFit-PVANow(range(1)),...
                                        'Color',pltCol);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Plot the average mean PVA strength - Figure 3 L,M

% Set movement thresholds
vRThresh = pi/10;
vFThresh = 0.1;

PVAStrenSummary = figure('units','normalized','outerposition',[0 0 0.25 0.5]);

% Calculate the mean PVA strenth when the fly is moving
allStrens = zeros(2,2,10,2);
for condID = 1:length(cond)
    subplot(1,3,condID);
    hold on;
    for flyID = 1:cond{condID}.numFlies

        trialNames = fields(cond{condID}.allFlyData{flyID});

        for trialType = 2:3
            for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialNames{trialType}))

                % Skip the trials where the behavior wasn't tracked
                if condID == 1 & trialType == 3 & flyID == 6 & trialID == 2
                    continue;
                end
                if condID == 2 & trialType == 2 & flyID == 7 & trialID == 1
                    continue;
                end
                
                % Get the data for this trial
                datNow = cond{condID}.allFlyData{flyID}.(trialNames{trialType}){trialID};

                % Calculate the PVA strength
                [angs, PVAPlt, PVAStren] = PVA(datNow.ROIaveMax-min(min(datNow.ROIaveMax)));

                % Pull the behavior
                vR = datNow.positionDatMatch.vRot;
                vR = sgolayfilt(vR,3,11);
                vF = datNow.positionDatMatch.vF;
                vF = sgolayfilt(vF,3,11);
                flyAct = union(find(abs(vR)>vRThresh),find(vF>vFThresh));

                tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                tPts = tPts-tPts(1);

                if length(strsplit(trialNames{trialType},'_'))==2
                    if trialID >=2
                        scatter(flyID,mean(PVAStren(flyAct)),150,'r','filled');
                        allStrens(condID,2,flyID,trialID-1) = mean(PVAStren(flyAct));
                    end
                else
                    scatter(flyID,mean(PVAStren(flyAct)),150,'b','filled')
                    allStrens(condID,1,flyID,trialID) = mean(PVAStren(flyAct));
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
    subplot(1,3,condID);
    title(cond{condID}.name)
end


% Split the data by condition
allCtrlPT = [];
allCtrl30 = [];
allD7PT = [];
allD730 = [];
for flyID = 1:10
    allCtrlPT = [allCtrlPT allStrens(1,1,flyID,1)];
    allCtrlPT = [allCtrlPT allStrens(1,1,flyID,2)];
    allCtrl30 = [allCtrl30 allStrens(1,2,flyID,1)];
    allCtrl30 = [allCtrl30 allStrens(1,2,flyID,2)];

    allD7PT = [allD7PT allStrens(2,1,flyID,1)];
    allD7PT = [allD7PT allStrens(2,1,flyID,2)];
    allD730 = [allD730 allStrens(2,2,flyID,1)];
    allD730 = [allD730 allStrens(2,2,flyID,2)];
end

% Calculate the stats
[h,p] = ttest(allCtrlPT(find(allCtrlPT>0)),allCtrl30(find(allCtrlPT>0)));
p
[h,p] = ttest(allD7PT(find(allD7PT>0)),allD730(find(allD730>0)));
p