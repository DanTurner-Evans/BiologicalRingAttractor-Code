%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\Shi\';

%% Load the data

cond{1}.name = 'Empty Gal4';
cond{1}.dirs{1} = strcat(dataDir,'TuBu-Ctrl\20190813');
cond{1}.dirs{2} = strcat(dataDir,'TuBu-Ctrl\20190815');
cond{1}.dirs{3} = strcat(dataDir,'TuBu-Ctrl\20190816');
cond{1}.dirs{4} = strcat(dataDir,'TuBu-Ctrl\20191125');
cond{1}.dirs{5} = strcat(dataDir,'TuBu-Ctrl\20191126');
cond{1}.dirs{6} = strcat(dataDir,'TuBu-Ctrl\20191203');
cond{1}.dirs{7} = strcat(dataDir,'TuBu-Ctrl\20191204');


cond{2}.name = '76B06';
cond{2}.dirs{1} = strcat(dataDir,'TuBu\20190813');
cond{2}.dirs{2} = strcat(dataDir,'TuBu\20190815');
cond{2}.dirs{3} = strcat(dataDir,'TuBu\20190816');
cond{2}.dirs{4} = strcat(dataDir,'TuBu\20191125');
cond{2}.dirs{5} = strcat(dataDir,'TuBu\20191126');
cond{2}.dirs{6} = strcat(dataDir,'TuBu\20191204');

cond = FlyDatLoad(1,cond);

%% Look at the offset distribution in the dark and in the initial closed loop period

trialName = 'All_30C';
PVAThresh = 0.075;
stripeMult = 360/240;
offsetBins = linspace(-pi,pi,17);
offsetCents = offsetBins;
offsetCents(end) = [];
offsetCents = offsetCents + 0.5*mean(diff(offsetCents));

offsetStats = figure('units','normalized','outerposition',[0 0 1 1]);
for condID = 1:length(cond)
    offsetDist = figure('units','normalized','outerposition',[0 0 1 1]);
    for flyID = 1:cond{condID}.numFlies
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.(trialName))
            datNow = cond{condID}.allFlyData{flyID}.(trialName){trialID};
            
            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            tPts = tPts - tPts(1);
            heading = pi/180*datNow.positionDatMatch.PosRotMatch;
            stripePos = datNow.positionDatMatch.OffsetRotMatch(:,2);
            
            DF = datNow.ROIaveMax-1;
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
            
            figure(offsetDist);
            % Plot the distributions
            subplot(length(cond{condID}.allFlyData{flyID}.(trialName)),cond{condID}.numFlies,flyID+cond{condID}.numFlies*(trialID-1));
            hold on;
            h1 = histogram(darkOffset,offsetBins);
            h1.FaceColor = 'k';
            h1.FaceAlpha = 0.1;
            
            h2 = histogram(CLOffset,offsetBins);
            h2.FaceColor = 'b';
            h2.FaceAlpha = 0.1;
            
            xlim([-pi pi]);
            ylim([0 600]);
            if flyID == 1
                if trialID == 1
                    text(-pi,700,cond{condID}.name,'FontSize',14);
                end
            end
            
            % Fit a von Mises function to each
%             F = @(A,kappa,mu,x)A*exp(-kappa*cos(x-mu))/(2*pi*besselj(0,kappa));
%             
            [Ndark,edges] = histcounts(darkOffset,offsetBins);
%             testFitDark = fit(offsetCents', Ndark', F, 'StartPoint',[250 1 0],'Lower',[0 0 -pi],'Upper',[1000 20 pi]);
%             fitDark = testFitDark.A*exp(-testFitDark.kappa*cos(offsetCents-testFitDark.mu))/(2*pi*besselj(0,testFitDark.kappa));
%             plot(offsetCents,fitDark,'k');
%             
            [NCL,edges] = histcounts(CLOffset,offsetBins);
%             testFitCL = fit(offsetCents', NCL', F, 'StartPoint',[250 1 0],'Lower',[0 0 -pi],'Upper',[1000 20 pi]);
%             fitCL = testFitCL.A*exp(-testFitCL.kappa*cos(offsetCents-testFitCL.mu))/(2*pi*besselj(0,testFitCL.kappa));
%             plot(offsetCents,fitCL,'b');
            
            figure(offsetStats);
            subplot(2,1,1);
            hold on;
%             scatter(flyID+10*(condID-1)-0.2,testFitDark.kappa,60,'k','filled');
            scatter(flyID+10*(condID-1)-0.2,circ_r(offsetCents',Ndark'),60,'k','filled');
%             scatter(flyID+10*(condID-1)+0.2,testFitCL.kappa,60,'b','filled');
            scatter(flyID+10*(condID-1)+0.2,circ_r(offsetCents',NCL'),60,'b','filled');
            alpha(0.4);
            xlim([0 21]);
            ylim([0 1]);
%             ylim([0 4]);

            [S s] = circ_var(offsetCents',Ndark');
            cond{condID}.allFlyData{flyID}.(trialName){trialID}.darkPVA = S;
            [S s] = circ_var(offsetCents',NCL');
            cond{condID}.allFlyData{flyID}.(trialName){trialID}.CLPVA = S;
%             cond{condID}.allFlyData{flyID}.(trialName){trialID}.darkPVA = circ_r(offsetCents',Ndark');
%             cond{condID}.allFlyData{flyID}.(trialName){trialID}.CLPVA = circ_r(offsetCents',NCL');
        end
    end     
end

%% Make a summary figure for the paper - Figure 6 G-J

%% Diagram the PVA - Figure 1 I

% Load the average, registered stack
%load('C:\Users\turnerevansd\Documents\RawAnalysis\Shi\EmptyGal4\20190813\Fly4_4-6days_6fx25957_ShixEmptyGal4_30C_All_00002_RegDat.mat');
tifName = 'Fly4_4-6days_6fx25957_ShixEmptyGal4_30C_All_00002.tif';
allPathnameNow = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\'; % Replace this with your dir
[stackMaxIntNow, stackMean] = ImDatLoadBigtiff(tifName,allPathnameNow,0);
frame2plt = 1;
meanRegStack = squeeze(stackMaxIntNow(:,:,frame2plt));%-mean(stackMaxIntNow,3));

paperFig = figure('units','normalized','outerposition',[0 0 1 1]);

% Rotate it
rotMeanStack = imrotate(meanRegStack,-70);

% Gaussian filter it
gaussianSize = [10 10];
gaussianSigma = 4;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
filtRotMeanStack = imfilter(rotMeanStack,Gxy,'replicate','conv');

% Plot it
xoffset = 60;
yoffset = 50;
pltWid = 200;
cMin = 0;
cMax = 120;
cMap = brewermap(64, 'Blues');

subplot(2,5,1);
hold on;
imagesc(imcrop(filtRotMeanStack,[xoffset yoffset pltWid pltWid]));
caxis([cMin cMax]);
colormap(cMap);
axis equal;
axis off;

%Load and plot the ROIs
load('C:\Users\turnerevansd\Documents\RawAnalysis\Shi\EmptyGal4\20190813\Fly4_4-6days_6fx25957_ShixEmptyGal4_30C_ReferenceROIs.mat');
for roiPlt = 1:2
    posNow = position{roiPlt};
    posNow = vertcat(posNow,posNow(1,:));
    patch('XData',posNow(:,1)-xoffset,'YData',posNow(:,2)-yoffset,'FaceColor','none');
end

subplot(2,5,2);
for i=1:64
    patch('XData',[0 0.1 0.1 0],'YData',[(i-1)/64 (i-1)/64 i/64 i/64],'FaceColor',cMap(i,:),'EdgeColor','none');
end
text(0.1,0,num2str(cMin));
text(0.1,1,num2str(cMax));
text(0.1,0.5,'F');
xlim([0 2]);
ylim([-0.5 1.5]);
axis off;

subplot(2,5,3);
hold on;
angs = linspace(-pi,pi,17);
innerD = 0.2;
outerD = 1;
cMin = 0;
cMax = 1.5;
pltData = cond{1}.allFlyData{2}.All_30C{2}.ROIaveMax(:,frame2plt)-1;
PVAVal = circ_mean(angs(1:end-1)'+mean(diff(angs)),pltData);
pltDataNorm = (pltData-cMin)./(cMax-cMin);

for i=1:16
    col = max(1,ceil(64*pltDataNorm(i)));
    angSmooth = linspace(angs(i),angs(i+1),5);
    patch('XData',-[innerD*cos(angSmooth) fliplr(outerD*cos(angSmooth))],...
        'YData',[innerD*sin(angSmooth) fliplr(outerD*sin(angSmooth))],...
        'FaceColor',cMap(col,:),'EdgeColor','k');
    plot([0 -1.2*pltDataNorm(i)*cos(0.5*(angs(i)+angs(i+1)))],[0 1.2*pltDataNorm(i)*sin(0.5*(angs(i)+angs(i+1)))],'Color','k');
end
axis equal;
axis off;
line([0 -1.2*cos(PVAVal)],[0 1.2*sin(PVAVal)],'Color',[0.25 0 0.5]);

subplot(2,5,4);
pltData = cond{1}.allFlyData{2}.All_30C{2}.ROIaveMax(:,frame2plt)-1;
pltDataNorm = (pltData-cMin)./(cMax-cMin);
angs = linspace(-pi,pi,17);
innerD = 0.2;
outerD = 1;
for i=1:16
    col = max(1,ceil(64*pltDataNorm(i)));
    patch('XData',[0 0.1 0.1 0],'YData',[(i-1)/17 (i-1)/17 i/17 i/17],'FaceColor',cMap(col,:),'EdgeColor','none');
end
axis equal;
axis off;
line([-0.02 0.12],[0.5*PVAVal/pi+0.5 0.5*PVAVal/pi+0.5],'Color',[0.25 0 0.5]);

subplot(2,5,5);
for i=1:64
    patch('XData',[0 0.1 0.1 0],'YData',[(i-1)/64 (i-1)/64 i/64 i/64],'FaceColor',cMap(i,:),'EdgeColor','none');
end
text(0.1,0,num2str(cMin));
text(0.1,1,num2str(cMax));
text(0.1,0.5,'DF/F');
xlim([0 2]);
ylim([-0.5 1.5]);
axis off;

%% Make a summary figure for the paper - plot the room temperature trials - Figure S9
paperFig = figure('units','normalized','outerposition',[0 0 1 1]);

trialName = 'All';

% Now plot example of the bump over time
flies = [3,2];
trials = [2,2];
for condID = 1:length(cond)
    for flyID = flies(condID)
        for trialID = trials(condID)
            datNow = cond{condID}.allFlyData{flyID}.(trialName){trialID};

            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            tPts = tPts - tPts(1);
            heading = pi/180*datNow.positionDatMatch.PosRotMatch;
            stripePos = 360/240*datNow.positionDatMatch.OffsetRotMatch(:,2);

            headingPlt = heading;
            stripePosPlt = UnWrap(stripePos,2,0);
%             for tStep = 2:length(tPts)
%                 if abs(headingPlt(tStep)-headingPlt(tStep-1))> pi
%                     headingPlt(tStep-1) = NaN;
%                 end
%                 if abs(stripePosPlt(tStep)-stripePosPlt(tStep-1))> pi
%                     stripePosPlt(tStep-1) = NaN;
%                 end
%             end

            DF = datNow.ROIaveMax-1;
            [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));
            angs = angs-0.5*mean(diff(angs));
            PVAUW = UnWrap(PVAPlt,1.5,1);
            PVAPlt(find(PVAStren<PVAThresh)) = NaN;
            for tStep = 2:length(tPts)
                if abs(PVAPlt(tStep)-PVAPlt(tStep-1))> pi
                    PVAPlt(tStep-1) = NaN;
                end
            end

            [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
            darkJumps = find(diff(darkPer)>1);
            if ~isempty(darkJumps);
                darkPer(darkJumps(1)+1:end) = [];
            end
            % Find when the stripe jumps
            stripeJumps = find(diff(datNow.positionDatMatch.Trans)~=0)+1;
            stripeJumps(end) = [];

            subplot(5,5,10*(condID-1)+[1:3]);
            hold on;
            imagesc(tPts,angs,DF);
            plot(tPts,PVAPlt,'Color',[0.25 0 0.5]);
            caxis([0 1.25]);
            ylim([-pi pi]);
            yticks([-pi 0 pi]);
            ylabel('EB pos. (rad)');
            xlim([tPts(1) tPts(stripeJumps(1))]);
            xlims = get(gca,'XLim');
            text(xlims(1),1.5*pi,cond{condID}.name,'FontSize',14);
            plot([tPts(darkPer(end)) tPts(darkPer(end))],[-pi pi],'--k');


            PVAUW(find(PVAStren<PVAThresh)) = NaN;
            subplot(5,5,10*(condID-1)+[6:8]);
            hold on;
            plot(tPts(darkPer),stripePosPlt(darkPer),'k');
            plot(tPts(CLPer),stripePosPlt(CLPer),'b');
            plot(tPts,PVAUW-PVAUW(darkPer(end))+stripePosPlt(darkPer(end)),'Color',[0.25 0 0.5]);
            xlabel('time (s)');
%             ylim([-pi pi]);
            ylabel('pos. (rad)');
            xlim([tPts(1) tPts(stripeJumps(1))]);
            yticks(linspace(-20*pi,20*pi,21));
            yl = get(gca,'ylim');
            plot([tPts(darkPer(end)) tPts(darkPer(end))],[yl(1) yl(2)],'--k');
            if condID == 1
                legend({'heading (dark)','stripe position (CL)'});
                legend('boxoff');
            end
        end
    end
    colormap(brewermap(64, 'Blues'));
end

%% Plot the bump statistics - bump width, bump amplitude, bump stability

PVAThresh = 0.075;

headingScale = 360/240;

minT = 10;
filtParam1 = 3;
filtParam2 = 11;

vFThresh = 0.1;
vRThresh = pi/10;

SummaryFig = figure('units','normalized','outerposition',[0 0 1 1]);
condIDs = [1,2];
flyIDs = [2,1];
trialIDs_Cold = [1,2];
trialIDs_Hot = [1,3];
visCond = 'CL';%CL

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
        
        datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};
                
        tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
        tPts = tPts-tPts(1);
        heading = headingScale*datNow.positionDatMatch.OffsetRotMatch(:,2);

        vR = datNow.positionDatMatch.vRot;
        vF = datNow.positionDatMatch.vF;
        flyStop = intersect(find(vR<vRThresh),find(vF<vFThresh));
        
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
        
        DF = sgolayfilt(datNow.ROIaveMax-1, filtParam1, filtParam2, [], 2);

        [angs, PVAPlt, PVAStren] = PVA(DF-min(min(DF)));
        
        PVAPaused = PVAPlt;
        for tPt = 1:length(flyStop)
            if PVAStren(flyStop(tPt))<PVAThresh & flyStop(tPt)>1
                PVAPaused(flyStop(tPt)) = PVAPaused(flyStop(tPt)-1);
            end
        end
        
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

                for trialID = trialStart:trialStop

                    datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};

                    [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                    darkEnd = find(diff(darkPer)>1);
                    if ~isempty(darkEnd)
                        darkPer(darkEnd+1:end) = [];
                    end

                    tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                    heading = headingScale*datNow.positionDatMatch.OffsetRotMatch(:,2);

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
                                
                                if flyID == flyIDs(condID) & visID == 2 & (mod(tStep,30) == 0) & hc == 1
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
                if flyID == flyIDs(condID) & visID == 2
                    subplot(6,6,12+3*condID-2+[0 6]);
                    plot([-1.5*pi 1.5*pi],[-1.5*pi 1.5*pi],'Color','k');
                    plot([-1.5*pi 1.5*pi],[0 0],'Color','k');
                    plot([0 0],[-1.5*pi 1.5*pi],'Color','k');
                    yticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi]);
                    xlim([-3*pi 3*pi]);
                    xticks([-2*pi -pi 0 pi 2*pi]);
                    xlabel('heading (rad)');
                    ylabel('PVA (rad)');
                    axis equal;
                end
                    
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
                alpha(0.5);
                ylim([0 1.5]);
                xlim([0.5 4.5]);
                ylabel('percentile diff');
                
                subplot(6,6,28+6*(visID-1));
                hold on;
                scatter(hc+2*(condID-1),median(allSlopes),50,pltCol,'filled');
                alpha(0.5);
                ylim([-0.5 2.5]);
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
                heading = headingScale*datNow.positionDatMatch.OffsetRotMatch(:,2);
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

visCond = 'Dark';
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
            allWidsNow = WidAmpDat{condID}.(visCond){flyID}.allWid;
            allWids(condID,hc,flyID) = mean(allWidsNow(~isnan(allWidsNow)));
        end
    end
end

for condID = 1:length(cond)
    for flyID = 1:cond{condID}.numFlies
        subplot(6,6,[26 27 32 33])
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

visCond = 'CL';
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
            allWidsNow = WidAmpDat{condID}.(visCond){flyID}.allWid;
            allWids(condID,hc,flyID) = mean(allWidsNow(~isnan(allWidsNow)));
        end
    end
end

for condID = 1:length(cond)
    for flyID = 1:cond{condID}.numFlies
        subplot(6,6,[29 30 35 36])
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

%% Plot the bump statistics - bump width, bump amplitude, bump stability - ctrl at RT only

PVAThresh = 0.075;

headingScale = 360/240;

minT = 10;
filtParam1 = 3;
filtParam2 = 11;

vFThresh = 0.1;
vRThresh = pi/10;

SummaryFig = figure('units','normalized','outerposition',[0 0 1 1]);

% First, chunk up the data by bouts where the bump can be tracked
RBins = linspace(-1,1,11);
slopeBins = linspace(-4,4,21);

for condID = 1
    for flyID = 1:cond{condID}.numFlies
        allPerc = zeros(2,2);
        allMedSlopes = zeros(2,2);
        for hc = 1
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

                for trialID = trialStart:trialStop

                    datNow = cond{condID}.allFlyData{flyID}.(trialType){trialID};

                    [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
                    darkEnd = find(diff(darkPer)>1);
                    if ~isempty(darkEnd)
                        darkPer(darkEnd+1:end) = [];
                    end

                    tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
                    heading = headingScale*datNow.positionDatMatch.OffsetRotMatch(:,2);

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
                subplot(2,3,1);
                hold on;
                scatter(visID+2*(condID-1),prctile(allSlopes,75)-prctile(allSlopes,25),50,pltCol,'filled');
                alpha(0.5);
                ylim([0 1]);
                xlim([0.5 2.5]);
                ylabel('percentile diff');
                
                subplot(2,3,2);
                hold on;
                scatter(visID+2*(condID-1),median(allSlopes),50,pltCol,'filled');
                alpha(0.5);
                ylim([0 1.5]);
                xlim([0.5 2.5]);
                ylabel('median slope');
                
                allPerc(hc,visID) = prctile(allSlopes,75)-prctile(allSlopes,25);
                allMedSlopes(hc,visID) = median(allSlopes);
            end
            subplot(2,3,1);
            plot([1 2],allPerc(1,:),'k');
            xticks([1 2]);
            xticklabels({'dark','CL'});
            subplot(2,3,2);
            plot([1 2],allMedSlopes(1,:),'k');
            plot([0.5 2.5],[1 1],'k');
            xticks([1 2]);
            xticklabels({'dark','CL'});
        end
    end
end

visCond = 'Dark';
allWids = zeros(2,2,10);
for hc = 1
    if hc == 1
        trialName = 'All';
    else
        trialName = 'All_30C';
    end
    WidAmpDat = BumpAmpAndWidth(cond,trialName,PVAThresh);
    
    for condID = 1:length(cond)
        for flyID = 1:cond{condID}.numFlies
            allWidsNow = WidAmpDat{condID}.(visCond){flyID}.allWid;
            allWids(condID,1,flyID) = mean(allWidsNow(~isnan(allWidsNow)));
        end
    end
end

for condID = 1
    for flyID = 1:cond{condID}.numFlies
        subplot(2,3,3)
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
end

visCond = 'CL';
for hc = 1
    if hc == 1
        trialName = 'All';
    else
        trialName = 'All_30C';
    end
    WidAmpDat = BumpAmpAndWidth(cond,trialName,PVAThresh);
    
    for condID = 1:length(cond)
        for flyID = 1:cond{condID}.numFlies
            allWidsNow = WidAmpDat{condID}.(visCond){flyID}.allWid;
            allWids(condID,2,flyID) = mean(allWidsNow(~isnan(allWidsNow)));
        end
    end
end

for condID = 1
    for flyID = 1:cond{condID}.numFlies
        subplot(2,3,3)
        hold on;
        scatter(2,allWids(condID,2,flyID),50,'b','filled');
        plot([1 2],allWids(condID,:,flyID),'Color','k');
        xlim([0.5 2.5]);
        xticks([1 2]);
        xticklabels({'dark','CL'});
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