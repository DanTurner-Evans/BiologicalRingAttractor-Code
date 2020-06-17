%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\TwoColor\';

%% Clear out old data and load the new data
cond{1}.name = 'EB: PEN2 - Green, PEG - Red';
cond{1}.dirs{1} = strcat(dataDir,'RedPEGsGreenPEN2s\EB\20180314');
cond{1}.dirs{2} = strcat(dataDir,'RedPEGsGreenPEN2s\EB\20180319');
cond{1}.dirs{3} = strcat(dataDir,'RedPEGsGreenPEN2s\EB\20180321');

cond{2}.name = 'PB: PEN2 - Green, PEG - Red';
cond{2}.dirs{1} = strcat(dataDir,'RedPEGsGreenPEN2s\PB\20180313');
cond{2}.dirs{2} = strcat(dataDir,'RedPEGsGreenPEN2s\PB\20180314');

cond{3}.name = 'EB: PEN2 - Red, PEG - Green';
cond{3}.dirs{1} = strcat(dataDir,'GreenPEGsRedPEN2s\EB\20191202');


cond = FlyDatLoad(2,cond);

%% Specify the PB glomeruli of interest
RLPB = [2:9];
RRPB = [10:17];
GLPB = [1:8];
GRPB = [11:18];

%% Sort the activity by the visual conditions

allAct = {};
% Look at the EB and PB
for condID = 1:length(cond)
   allAct{condID}.name = cond{condID}.name;
   
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
       
      % Sort the data across trials
      for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)
         
         [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch);
         
         % Extract the behavioral parameters
         vR = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.vRot;
         vF = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.vF;
         tPts = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.OffsetRotMatch(:,1);
         stripePos = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.OffsetRotMatch(:,2);
         heading = cond{condID}.allFlyData{flyID}.All{trialID}.positionDatMatch.PosRotMatch;
         
         % Get rid of the last time point so that I have velocities
         % throughout
         if CWPer(end) > length(vR)
             CWPer(end) = [];
         elseif CCWPer(end) > length(vR)
             CCWPer(end) = [];
         elseif darkPer(end) > length(vR)
             darkPer(find(darkPer > OLPer(end))) = [];
         end
         
         % Sort the data by period and color
         for periodID = 1:4
             
             for colorID = 1:2
                 if periodID == 1
                     perNow = darkPer;
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type = 'dark';
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color = [0 0 0];
                 elseif periodID == 2
                     perNow = CLPer;
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type = 'CL';
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color = [0 0 1];
                 elseif periodID == 3
                     perNow = CWPer;
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type = 'CW';
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color = [0.5 0 1];
                 elseif periodID == 4
                     perNow = CCWPer;
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type = 'CCW';
                     allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color = [0 0.5 1];
                 end
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vR = vR(perNow);
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vF = vF(perNow);
                 
                 % Find the time points, stripe position, and heading
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.tPts = tPts(perNow);
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.stripePos = stripePos(perNow);
                 allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.heading = heading(perNow);
             end
             % Find the activity
             if contains(cond{condID}.name,'EB: PEN2 - Green, PEG - Red') || contains(allAct{condID}.name,'EB: PEN2 - Red, PEG - Green')
                 allAct{condID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.act = ...
                     cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax(:,perNow);
                 allAct{condID}.fly{flyID}.color{2}.period{periodID}.trial{trialID}.act = ...
                     cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax(:,perNow);
             elseif contains(cond{condID}.name,'PB: PEN2 - Green, PEG - Red')
                 for RL = 1:2
                     if RL == 1
                         RROI = RLPB;
                         GROI = GLPB;
                     else
                         RROI = RRPB;
                         GROI = GRPB;
                     end
                     allAct{condID}.fly{flyID}.color{1}.period{periodID}.trial{trialID}.PBSide{RL}.act =...
                         cond{condID}.allFlyData{flyID}.All{trialID}.RROIaveMax(RROI,perNow);
                     allAct{condID}.fly{flyID}.color{2}.period{periodID}.trial{trialID}.PBSide{RL}.act =...
                         cond{condID}.allFlyData{flyID}.All{trialID}.GROIaveMax(GROI,perNow);
                 end
             end
         end
      end
   end
end

%% Plot example trials - Figure 7 G,H and Figure S11 P,Q
PVAThresh = 0.1;

% Specify the EB example
EBFly = 5;
EBTrial = 2;

% Specify the PB example
PBFly = 3;
PBTrial = 3;

angs = linspace(-pi,pi,17);
angs(end) = [];
angs = angs+0.5*mean(diff(angs));

for condID = 1:2
    
    actEx = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Specify EB and PB specifics
    if condID == 1
        flyID = EBFly;
        trialID = EBTrial;
        yLimits = [-pi pi];
    else
        flyID = PBFly;
        trialID = PBTrial;
        angs = [1:18];
        yLimits = [1 18];
    end
    
    datNow = cond{condID}.allFlyData{flyID}.All{trialID};
    
    tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
    tPts = tPts - tPts(1);
    heading = datNow.positionDatMatch.OffsetRotMatch(:,2);

    [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
    darkJump = find(diff(darkPer)>1);
    if ~isempty(darkJump)
        darkPer(darkJump(1)+1:end) = [];
    end
    
    if condID == 1
        stripePos = -stripePos;
        heading = -heading;
    end

    subplot(4,3,[1:2])
    if condID == 1
        RAct = datNow.RROIaveMax-1;
        [angs, PVAPlt_R, PVAStren_R] = PVA(RAct-min(min(RAct)));
        angs = angs-0.5*mean(diff(angs));
        PVAPlt_R(find(PVAStren_R<PVAThresh)) = NaN;
    else
        RActLPB = datNow.RROIaveMax(RLPB,:)-1;
        [angsPB, PVAPlt_RL, PVAStren_RL] = PVA(RActLPB-min(min(RActLPB)));
        PVAPlt_RL(find(PVAStren_RL<PVAThresh)) = NaN;
        RActRPB = datNow.RROIaveMax(RRPB,:)-1;
        [angsPB, PVAPlt_RR, PVAStren_RR] = PVA(RActRPB-min(min(RActRPB)));
        PVAPlt_RR(find(PVAStren_RR<PVAThresh)) = NaN;
        RAct = zeros([20, length(RActLPB)]);
        RAct(RLPB,:) = RActLPB;
        RAct(RRPB+2,:) = RActRPB;
    end
    RIm = zeros([size(RAct) 3]);
    RIm(:,:,1) = RAct./1.5;
    RIm(:,:,3) = RAct./1.5;
    image(tPts,angs,RIm);
    xlim([tPts(darkPer(1)) tPts(CLPer(end))]);
    if condID == 1
        ylabel('EB position');
    else
        ylabel('R PB        L PB');
    end
    title('P-EG act.');
    
    subplot(4,3,3);
    mgtaCbar = zeros(64,2,3);
    mgtaCbar(:,1,1) = linspace(0,1,64);
    mgtaCbar(:,2,1) = linspace(0,1,64);
    mgtaCbar(:,1,3) = linspace(0,1,64);
    mgtaCbar(:,2,3) = linspace(0,1,64);
    image([1:2],[1:64],flipud(mgtaCbar));
    xlim([1 20]);
    text(3,1,'1.5');
    text(3,64,'0');
    text(3,32,'DF/F');
    axis off;

    subplot(4,3,[4:5])
    if condID == 1
        GAct = datNow.GROIaveMax-1;
        [angs, PVAPlt_G, PVAStren_G] = PVA(GAct-min(min(GAct)));
        angs = angs-0.5*mean(diff(angs));
        PVAPlt_G(find(PVAStren_G<PVAThresh)) = NaN;
    else
        GActLPB = datNow.GROIaveMax(GLPB,:)-1;
        [angsPB, PVAPlt_GL, PVAStren_GL] = PVA(GActLPB-min(min(GActLPB)));
        PVAPlt_GL(find(PVAStren_GL<PVAThresh)) = NaN;
        GActRPB = datNow.GROIaveMax(GRPB,:)-1;
        [angsPB, PVAPlt_GR, PVAStren_GR] = PVA(GActRPB-min(min(GActRPB)));
        PVAPlt_GR(find(PVAStren_GR<PVAThresh)) = NaN;
        GAct = zeros([20, length(GActLPB)]);
        GAct(GLPB,:) = GActLPB;
        GAct(GRPB+2,:) = GActRPB;
    end
    GIm = zeros([size(GAct) 3]);
    GIm(:,:,2) = GAct./1.5;
    image(tPts,angs,GIm);
    xlim([tPts(darkPer(1)) tPts(CLPer(end))]);
    if condID == 1
        ylabel('EB position');
    else
        ylabel('R PB        L PB');
    end
    title('P-EN2 act.');

    subplot(4,3,6);
    gCbar = zeros(64,2,3);
    gCbar(:,1,2) = linspace(0,1,64);
    gCbar(:,2,2) = linspace(0,1,64);
    image([1:2],[1:64],flipud(gCbar));
    xlim([1 20]);
    text(3,1,'1.5');
    text(3,64,'0');
    text(3,32,'DF/F');
    axis off;

    subplot(4,3,[7:8])
    overlayIm = zeros([size(RAct) 3]);
    overlayIm(:,:,1) = RAct;
    overlayIm(:,:,3) = RAct;
    overlayIm(:,:,2) = GAct;
    image(tPts,angs,overlayIm);
    xlim([tPts(1) tPts(end)]);
    if condID == 1
        ylabel('EB position');
    else
        ylabel('R PB        L PB');
    end
    
    headingPlt = heading;
    for tPt = 2:length(tPts)
        if abs(headingPlt(tPt) - headingPlt(tPt-1)) > pi
            headingPlt(tPt-1) = NaN;
        end
        if condID == 1
            if abs(PVAPlt_R(tPt) - PVAPlt_R(tPt-1)) > pi
                PVAPlt_R(tPt-1) = NaN;
            end
            if abs(PVAPlt_G(tPt) - PVAPlt_G(tPt-1)) > pi
                PVAPlt_G(tPt-1) = NaN;
            end
        else
            if abs(PVAPlt_RR(tPt) - PVAPlt_RR(tPt-1)) > pi
                PVAPlt_RR(tPt-1) = NaN;
            end
            if abs(PVAPlt_RL(tPt) - PVAPlt_RL(tPt-1)) > pi
                PVAPlt_RL(tPt-1) = NaN;
            end
            if abs(PVAPlt_GR(tPt) - PVAPlt_GR(tPt-1)) > pi
                PVAPlt_GR(tPt-1) = NaN;
            end 
            if abs(PVAPlt_GL(tPt) - PVAPlt_GL(tPt-1)) > pi
                PVAPlt_GL(tPt-1) = NaN;
            end 
        end
        
    end
    subplot(4,3,[10:11])
    hold on;
    plot(tPts(darkPer),headingPlt(darkPer),'Color','k');
    plot(tPts(CLPer),headingPlt(CLPer),'Color','b');
    if condID == 1
        plot(tPts,-PVAPlt_R,'m');
        plot(tPts,-PVAPlt_G,'g');
        legend({'dark','closed-loop stripe','P-EG PVA','E-PG PVA'});
    else
        plot(tPts,-PVAPlt_RR,'m');
        plot(tPts,-PVAPlt_RL,'m');
        plot(tPts,-PVAPlt_GR,'g');
        plot(tPts,-PVAPlt_GL,'g');
        legend({'dark','closed-loop stripe','P-EG PVA (R)','P-EG PVA (L)','E-PG PVA (R)','E-PG PVA (L)'});
    end
    ylim([-pi pi]);
    xlim([tPts(darkPer(1)) tPts(CLPer(end))]);
    xlabel('time (s)');
    ylabel('position (rad)');
    
%         set(actEx,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%         print(actEx,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\EMPaper2ColorCode\',...
%         allAct{regionID}.name,'_',...
%         allAct{regionID}.fly{flyID}.color{1}.period{periodID}.type,'_Ex'),'-dpdf');
end

%% Run a regression with the peak maximum against the forward and rotational velocities - for each trial
numPts = 41;

% Look at the EB and PB
for condID = 1:length(cond)
   reg = figure('units','normalized','outerposition',[0 0 1 1]);
   
   numFlies = length(allAct{condID}.fly);
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
       
      % Sort the data across trials
      for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)
          
          numPeriods = length(allAct{condID}.fly{flyID}.color{1}.period);
          
          % Sort the data by period
          for periodID = 1:4
              
              % Sort the data by color
              for colorID = 1:2
                  
                  % Find the time step
                  tStep = mean(diff(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.tPts));
                  
                  % Get the variables of choice
                  vR = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vR;
                  vF = allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vF;
                  
                  % Get the maximum activity
                  if contains(allAct{condID}.name,'EB: PEN2 - Green, PEG - Red') || contains(allAct{condID}.name,'EB: PEN2 - Red, PEG - Green')
                      maxAct = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.act);
                      
                      autoC = zeros(numPts,1);
                      for lag = 1:numPts
                          lagCCs = corrcoef(...
                              abs(vR(ceil(numPts/2):end-floor(numPts/2)-1)),...
                              maxAct(numPts-lag+1:end-lag));
                          autoC(lag) = lagCCs(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vRAutoC = autoC;
                      
                      subplot(2*numFlies,2*numPeriods,4*numPeriods*(flyID-1)+2*periodID-2+colorID);
                      hold on;
                      plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoC,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('|vR| autocorrelation');
                          text(-2*tStep*numPts/2,1.1,strcat('fly #',num2str(flyID)));
                          if flyID == 1
                              text(-tStep*numPts/2,1.25,allAct{condID}.name);
                          end
                      end
                      if flyID == 1
                          if colorID == 1
                              title(strcat(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type,'-red'),...
                                  'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                          else
                              title(strcat(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type,'-green'),...
                                  'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                          end
                      end
                      
                      autoC = zeros(numPts,1);
                      for lag = 1:numPts
                          lagCCs = corrcoef(...
                              vF(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxAct(numPts-lag+1:end-lag));
                          autoC(lag) = lagCCs(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vFAutoC = autoC;
                      
                      subplot(2*numFlies,2*numPeriods,4*numPeriods*(flyID-1)+2*periodID-2+colorID+2*numPeriods);
                      hold on;
                      plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoC,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('vF autocorrelation');
                      end
                      if flyID == cond{condID}.numFlies
                          xlabel('time delay (s)');
                      end
                      
                      
                  elseif contains(allAct{condID}.name,'PB: PEN2 - Green, PEG - Red')
                      maxActL = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.act);
                      maxActR = max(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.act);
                      
                      vRPos = zeros(size(vR));
                      vRPos(find(vR>0)) = vR(find(vR > 0));
                      autoCL = zeros(numPts,1);
                      autoCR = zeros(numPts,1);
                      for lag = 1:numPts
                          
                          lagCCsL = corrcoef(...
                              vRPos(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActL(numPts-lag+1:end-lag));
                          autoCL(lag) = lagCCsL(2,1);
                          lagCCsR = corrcoef(...
                              vRPos(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActR(numPts-lag+1:end-lag));
                          autoCR(lag) = lagCCsR(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRPosAutoC = autoCL;
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRPosAutoC = autoCR;
                      
                      subplot(3*numFlies,2*numPeriods,6*numPeriods*(flyID-1)+2*periodID-2+colorID);
                      hold on;
                      Lplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCL,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      Rplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCR,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                          'LineWidth',2,'LineStyle','--');
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('vR > 0 autocorrelation');
                            text(-2*tStep*numPts/2,1.1,strcat('fly #',num2str(flyID)));
                            if flyID == 1
                                text(-tStep*numPts/2,1.25,allAct{condID}.name);
                                legend([Lplt Rplt],{'L PB','R PB'});
                            end
                      end
                      if flyID == 1
                          if colorID == 1
                              title(strcat(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type,'-red'),...
                                  'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                          else
                              title(strcat(allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.type,'-green'),...
                                  'Color',allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                          end
                      end
                      
                      vRNeg = zeros(size(vR));
                      vRNeg(find(vR < 0)) = -vR(find(vR < 0));
                      autoCL = zeros(numPts,1);
                      autoCR = zeros(numPts,1);
                      for lag = 1:numPts
                          
                          lagCCsL = corrcoef(...
                              vRNeg(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActL(numPts-lag+1:end-lag));
                          autoCL(lag) = lagCCsL(2,1);
                          lagCCsR = corrcoef(...
                              vRNeg(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActR(numPts-lag+1:end-lag));
                          autoCR(lag) = lagCCsR(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRNegAutoC = autoCL;
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRNegAutoC = autoCR;
                      
                      subplot(3*numFlies,2*numPeriods,6*numPeriods*(flyID-1)+2*periodID-2+colorID+2*numPeriods);
                      hold on;
                      Lplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCL,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      Rplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCR,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                          'LineWidth',2,'LineStyle','--');
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('vR < 0 autocorrelation');
                      end
                      
                      autoCL = zeros(numPts,1);
                      autoCR = zeros(numPts,1);
                      for lag = 1:numPts
                          lagCCsL = corrcoef(...
                              vF(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActL(numPts-lag+1:end-lag));
                          autoCL(lag) = lagCCsL(2,1);
                          lagCCsR = corrcoef(...
                              vF(ceil(numPts/2):end-floor(numPts/2)-1),...
                              maxActR(numPts-lag+1:end-lag));
                          autoCR(lag) = lagCCsR(2,1);
                      end
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vFAutoC = autoCL;
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vFAutoC = autoCR;
                      
                      subplot(3*numFlies,2*numPeriods,6*numPeriods*(flyID-1)+2*periodID-2+colorID+4*numPeriods);
                      hold on;
                      Lplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCL,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      Rplt = plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),autoCR,'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                          'LineWidth',2,'LineStyle','--');
                      line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                      line([0 0],[-0.5 1], 'Color', 'k','LineStyle','--');
                      ylim([-0.5 1]);
                      xlim([-tStep*numPts/2 tStep*numPts/2]);
                      if periodID == 1 & colorID == 1
                          ylabel('vF autocorrelation');
                      end
                      if flyID == cond{condID}.numFlies
                          xlabel('time delay (s)');
                      end
                  end
                  
              end
          end
      end
   end
end
          
%% Run a regression with the peak maximum against the forward and rotational velocities - averages

% Look at the EB and PB
for condID = 1:length(cond)
   reg = figure('units','normalized','outerposition',[0 0 1 1]);
   
   numFlies = length(allAct{condID}.fly);
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
      % Sort the data by period
      for periodID = 1:4

          % Sort the data by color
          for colorID = 1:2
              if contains(allAct{condID}.name,'EB: PEN2 - Green, PEG - Red') || contains(allAct{condID}.name,'EB: PEN2 - Red, PEG - Green')
                  vRAll = [];
                  vFAll = [];
              elseif contains(allAct{condID}.name,'PB: PEN2 - Green, PEG - Red')
                  vRPosLAll = [];
                  vRNegLAll = [];
                  vFLAll = [];
                  vRPosRAll = [];
                  vRNegRAll = [];
                  vFRAll = [];
              end
              % Sort the data across trials
              for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                  % Group across trials
                  if contains(allAct{condID}.name,'EB: PEN2 - Green, PEG - Red') || contains(allAct{condID}.name,'EB: PEN2 - Red, PEG - Green')
                      if condID == 3 && flyID == 4 && trialID == 4
                          break;
                      end
                      vRAll = horzcat(vRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vRAutoC);
                      vFAll = horzcat(vFAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vFAutoC);
                  elseif contains(allAct{condID}.name,'PB: PEN2 - Green, PEG - Red') 
                      vRPosLAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRPosAutoC);
                      vRPosRAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRPosAutoC);

                      vRNegLAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRNegAutoC);
                      vRNegRAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRNegAutoC);

                      vFLAll = horzcat(vFLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vFAutoC);
                      vFRAll = horzcat(vFRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vFAutoC);
                  end
              end
              if contains(allAct{condID}.name,'EB: PEN2 - Green, PEG - Red') || contains(allAct{condID}.name,'EB: PEN2 - Red, PEG - Green')
                    subplot(2,2*numFlies,2*(flyID-1)+colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRAll,2),'Color',...
                    allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    if colorID == 1
                        ylabel('|vR| autocorrelation');
                        if flyID == 1
                            text(-tStep*numPts/2,0.55,allAct{condID}.name);
                            if periodID == 4
                               legend({'dark','CL','CW','CCW'});
                            end
                        end
                        title(strcat('fly #',num2str(flyID),'-red'));
                    else
                        title(strcat('fly #',num2str(flyID),'-green'));
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 0.5]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end
                    
                    subplot(2,2*numFlies,2*(flyID-1)+colorID+2*numFlies);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                    line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                    ylim([-0.25 0.5]);
                    xlim([-tStep*numPts/2 tStep*numPts/2]);
                    if colorID == 1 & flyID == 1
                        ylabel('vF autocorrelation');
                    end
                    xlabel('time delay (s)');
              elseif contains(allAct{condID}.name,'PB: PEN2 - Green, PEG - Red')
                    subplot(3,2*numFlies,2*(flyID-1)+colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    if colorID == 1
                        ylabel('vR > 0 autocorrelation');
                        if flyID == 1
                            text(-tStep*numPts/2,0.6,allAct{condID}.name);
                        end
                        title(strcat('fly #',num2str(flyID),'-red'));
                    else
                        title(strcat('fly #',num2str(flyID),'-green'));
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 0.5]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end
                    
                    subplot(3,2*numFlies,2*(flyID-1)+colorID+2*numFlies);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    if colorID == 1
                        ylabel('vR < 0 autocorrelation');
                        title(strcat('fly #',num2str(flyID),'-red'));
                    else
                        title(strcat('fly #',num2str(flyID),'-green'));
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 0.5]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end

                    subplot(3,2*numFlies,2*(flyID-1)+colorID+4*numFlies);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFRAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFLAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                    line([0 0],[-0.25 0.5], 'Color', 'k','LineStyle','--');
                    ylim([-0.25 0.5]);
                    xlim([-tStep*numPts/2 tStep*numPts/2]);
                    if colorID == 1 & flyID == 1
                        ylabel('vF autocorrelation');
                    end
                    xlabel('time delay (s)');
              end
          end
      end
   end   
%    set(reg,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%     print(reg,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\RSS2191G20739\',...
%     allAct{condID}.name,'\',allAct{condID}.name,'_MeanAutoCC'),'-dpdf');
end

%% Plot the regressions nicely - Figure S11 A-D
% For each condition, pick and plot an example fly, and the correlation
% offset and delay for each condition

EBExFly = 5;
PBExFly = 3;

% Look at the EB and PB
for condID = 1:length(cond)
   reg = figure('units','normalized','outerposition',[0 0 1 1]);
   
   numFlies = length(allAct{condID}.fly);
   % Step through the flies
   for flyID = 1:cond{condID}.numFlies
      % Sort the data by period
      for periodID = 1:4

          % Sort the data by color
          for colorID = 1:2
              if contains(allAct{condID}.name,'EB: PEN2 - Green, PEG - Red') || contains(allAct{condID}.name,'EB: PEN2 - Red, PEG - Green')
                  vRAll = [];
                  vFAll = [];
              elseif contains(allAct{condID}.name,'PB: PEN2 - Green, PEG - Red')
                  vRPosLAll = [];
                  vRNegLAll = [];
                  vFLAll = [];
                  vRPosRAll = [];
                  vRNegRAll = [];
                  vFRAll = [];
              end
              % Sort the data across trials
              for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

                  % Group across trials
                  if contains(allAct{condID}.name,'EB: PEN2 - Green, PEG - Red') || contains(allAct{condID}.name,'EB: PEN2 - Red, PEG - Green')
                      if condID == 3 && flyID == 4 && trialID == 4
                          break;
                      end
                      vRAll = horzcat(vRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vRAutoC);
                      vFAll = horzcat(vFAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.vFAutoC);
                  elseif contains(allAct{condID}.name,'PB: PEN2 - Green, PEG - Red') 
                      vRPosLAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRPosAutoC);
                      vRPosRAll = horzcat(vRPosLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRPosAutoC);

                      vRNegLAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vRNegAutoC);
                      vRNegRAll = horzcat(vRNegLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vRNegAutoC);

                      vFLAll = horzcat(vFLAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{1}.vFAutoC);
                      vFRAll = horzcat(vFRAll,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.trial{trialID}.PBSide{2}.vFAutoC);
                  end
              end
              if contains(allAct{condID}.name,'EB: PEN2 - Green, PEG - Red') || contains(allAct{condID}.name,'EB: PEN2 - Red, PEG - Green')
                  if flyID == EBExFly
                      subplot(2,4,colorID);
                      hold on;
                      plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRAll,2),'Color',...
                          allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      ax = gca;
                      if colorID == 1
                          ylabel('|vR| autocorrelation');
                          text(-tStep*numPts/2-1,1.15,strcat(allAct{condID}.name,'-Summary'),'FontSize',14);
                          if periodID == 4
                             legend({'dark','CL','CW','CCW'});
                             legend('boxoff');
                          end
                          title('red channel');
                          ax.XColor = [0.5 0 0.5];
                          ax.YColor = [0.5 0 0.5];
                      else
                          title('green channel');
                          ax.XColor = [0 0.5 0];
                          ax.YColor = [0 0.5 0];
                      end
                      if periodID == 4
                          line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                          line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                          ylim([-0.25 1]);
                          xlim([-tStep*numPts/2 tStep*numPts/2]);
                      end
                                          
                      subplot(2,4,4+colorID);
                      hold on;
                      plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                      if periodID == 4
                          line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                          line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                          ylim([-0.25 1]);
                          xlim([-tStep*numPts/2 tStep*numPts/2]);
                      end
                      ax = gca;
                      if colorID == 1
                          ylabel('vF autocorrelation');
                          ax.XColor = [0.5 0 0.5];
                          ax.YColor = [0.5 0 0.5];
                      else
                          ax.XColor = [0 0.5 0];
                          ax.YColor = [0 0.5 0];
                      end
                      xlabel('time delay (s)');
                  end
                           
                  for vVar = 1:2
                      if vVar ==1 
                          maxCC = max(mean(vRAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRAll,2)==maxCC)-1)*tStep;
                          ylab = 'vR CC';
                          if colorID == 1
                              pastMax1 = maxCC;
                          end
                      else
                          maxCC = max(mean(vFAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vFAll,2)==maxCC)-1)*tStep;
                          ylab = 'vF CC';
                          if colorID == 1
                              pastMax2 = maxCC;
                          end
                      end
                      
                          
                      subplot(2,4,2+vVar);
                      hold on;
%                       scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
%                           40,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color, 'filled');
                      if colorID == 1
                        scatter(0.25*(colorID-1.5)+periodID,maxCC,...
                              40,[0.5 0 0.5],'filled');
                      else
                          scatter(0.25*(colorID-1.5)+periodID,maxCC,...
                              40,[0 0.5 0],'filled');
                          if vVar == 1
                              plot([0.25*(colorID-2.5)+periodID 0.25*(colorID-1.5)+periodID],[pastMax1 maxCC],'k');
                          else
                              plot([0.25*(colorID-2.5)+periodID 0.25*(colorID-1.5)+periodID],[pastMax2 maxCC],'k');
                          end
                      end
        %                   alpha(0.5);
                      ylim([-0 1]);
                      ylabel(ylab);
                      set(gca,'XTick',[1:4],'XTickLabels',{'dark','CL','CW','CCW'});
                      line([0 5],[0 0],'Color','k','LineStyle','--');

                  end
                  
              elseif contains(allAct{condID}.name,'PB: PEN2 - Green, PEG - Red')
                  if flyID == PBExFly
                    subplot(3,6,colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRPosLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    ax = gca;
                    if colorID == 1
                        ylabel('vR > 0 autocorrelation');
                        text(-tStep*numPts/2-1,1.15,strcat(allAct{condID}.name,'-Summary'));
                        title('red channel');
                        ax.XColor = [0.5 0 0.5];
                        ax.YColor = [0.5 0 0.5];
                    else
                        title('green channel');
                        ax.XColor = [0 0.5 0];
                        ax.YColor = [0 0.5 0];
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 1]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                        if colorID == 1
                            legend({'dark (R)','dark (L)','CL (R)','CL(L)','CW (R)','CW (L)','CCW (R)','CCW (L)'});
                            legend('boxoff');
                        end
                    end
                    
                    subplot(3,6,6+colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegRAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vRNegLAll,2),'Color',...
                        allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    ax = gca;
                    if colorID == 1
                        ylabel('vR < 0 autocorrelation');
                        ax.XColor = [0.5 0 0.5];
                        ax.YColor = [0.5 0 0.5];
                    else
                        ax.XColor = [0 0.5 0.5];
                        ax.YColor = [0 0.5 0.5];
                    end
                    if periodID == 4
                        line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                        line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                        ylim([-0.25 1]);
                        xlim([-tStep*numPts/2 tStep*numPts/2]);
                    end

                    subplot(3,6,12+colorID);
                    hold on;
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFRAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color);
                    plot(-tStep*numPts/2+linspace(0,(numPts-1)*tStep,numPts),mean(vFLAll,2),'Color',...
                      allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color,...
                      'LineWidth',2,'LineStyle','--');
                    line([-tStep*numPts/2 tStep*numPts/2],[0 0], 'Color', 'k','LineStyle','--');
                    line([0 0],[-0.25 1], 'Color', 'k','LineStyle','--');
                    ylim([-0.25 1]);
                    xlim([-tStep*numPts/2 tStep*numPts/2]);
                    ax = gca;
                    if colorID == 1
                        ylabel('vF autocorrelation');
                        ax.XColor = [0.5 0 0.5];
                        ax.YColor = [0.5 0 0.5];
                    else
                        ax.XColor = [0 0.5 0];
                        ax.YColor = [0 0.5 0];
                    end
                    xlabel('time delay (s)');
                  end
                  
                  for vVar = 1:6
                      if vVar == 1 
                          maxCC = max(mean(vRPosLAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRPosLAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 2
                          maxCC = max(mean(vRNegLAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRNegLAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 3
                          maxCC = max(mean(vFLAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vFLAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 4
                          maxCC = max(mean(vRPosRAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRPosRAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 5
                          maxCC = max(mean(vRNegRAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vRNegRAll,2)==maxCC)-1)*tStep;
                      elseif vVar == 6
                          maxCC = max(mean(vFRAll,2));
                          CClag = -tStep*numPts/2+(find(mean(vFRAll,2)==maxCC)-1)*tStep;
                      end
                      
                          
                      subplot(3,6,3+6*mod(vVar-1,3)+floor((vVar-1)/3));
                      hold on;
%                       scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
%                           30,allAct{condID}.fly{flyID}.color{colorID}.period{periodID}.color, 'filled');
                      if colorID == 1
                        scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
                              30,[0.5 0 0.5],'filled');
                      else
                          scatter(0.25*(colorID-1.5)+periodID+0.05*(flyID-ceil(cond{condID}.numFlies/2)),maxCC,...
                              30,[0 0.5 0],'filled');
                      end
        %                   alpha(0.5);
                      ylim([-0.25 1]);
                      set(gca,'XTick',[1:4],'XTickLabels',{'dark','CL','CW','CCW'});
                      line([0 5],[0 0],'Color','k','LineStyle','--');
                      if vVar < 4
                          title('left PB');
                      else
                          title('right PB');
                      end
                  end
              end
          end
      end
   end
%    reg.Renderer='Painters';
%    set(reg,'PaperPositionMode','manual','PaperOrientation','landscape','PaperUnits','inches','PaperPosition',[0 0 11 8.5]);
%    print(reg,strcat('C:\Users\turnerevansd\Documents\RawAnalysis\EMPaper2ColorCode\',...
%    allAct{condID}.name,'_CCSummary'),'-dpdf');
end

%% In the EB, look at the PVA diff as a function of velocity - Figure 7 I,J
PVAThresh = 0.1;

sgolayOrder = 3;
sgolayFrames = 11;

vREdges = linspace(-pi/2,pi/2,9);
vRCents = vREdges;
vRCents(end) = [];
vRCents = vRCents+0.5*(mean(diff(vRCents)));

condID = 1;

figure;

WidAmpDat = BumpAmpAndWidth2Cols(cond,'All',PVAThresh);
for condID = [1,3]
    for flyID = 1:cond{condID}.numFlies

        allDarkvR = [];
        allCLvR = [];
        allCWvR = [];
        allCCWvR = [];
        allDarkDiff = [];
        allCLDiff = [];
        allCWDiff = [];
        allCCWDiff = [];
        for trialID = 1:length(cond{condID}.allFlyData{flyID}.All)

            datNow = cond{condID}.allFlyData{flyID}.All{trialID};

            tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
            tPts = tPts - tPts(1);
            vR = datNow.positionDatMatch.vRot;
            vR = sgolayfilt(vR,sgolayOrder,sgolayFrames);

            [darkPer,OLPer,CLPer,CWPer,CCWPer] = SortVis(datNow.positionDatMatch);
            darkJump = find(diff(darkPer)>1);
            if ~isempty(darkJump)
                darkPer(darkJump(1)+1:end) = [];
            end

            RAct = datNow.RROIaveMax-1;
            [angs, PVAPlt_R, PVAStren_R] = PVA(RAct-min(min(RAct)));
            GAct = datNow.GROIaveMax-1;
            [angs, PVAPlt_G, PVAStren_G] = PVA(GAct-min(min(GAct)));

            for tPt = 1:length(darkPer)
                if PVAStren_R(darkPer(tPt)) > PVAThresh & PVAStren_G(darkPer(tPt)) > PVAThresh
                    PVADiff = min(abs(PVAPlt_R(darkPer(tPt))-PVAPlt_G(darkPer(tPt))),...
                        2*pi-abs(PVAPlt_R(darkPer(tPt))-PVAPlt_G(darkPer(tPt))));
                    allDarkvR = [allDarkvR vR(darkPer(tPt))];
                    allDarkDiff = [allDarkDiff PVADiff];
                end
            end

            for tPt = 1:length(CLPer)
                if PVAStren_R(CLPer(tPt)) > PVAThresh & PVAStren_G(CLPer(tPt)) > PVAThresh
                    PVADiff = min(abs(PVAPlt_R(CLPer(tPt))-PVAPlt_G(CLPer(tPt))),...
                        2*pi-abs(PVAPlt_R(CLPer(tPt))-PVAPlt_G(CLPer(tPt))));
                    allCLvR = [allCLvR vR(CLPer(tPt))];
                    allCLDiff = [allCLDiff PVADiff];
                end
            end

            for tPt = 1:length(CWPer)
                if PVAStren_R(CWPer(tPt)) > PVAThresh & PVAStren_G(CWPer(tPt)) > PVAThresh
                    PVADiff = min(abs(PVAPlt_R(CWPer(tPt))-PVAPlt_G(CWPer(tPt))),...
                        2*pi-abs(PVAPlt_R(CWPer(tPt))-PVAPlt_G(CWPer(tPt))));
                    allCWvR = [allCWvR vR(CWPer(tPt))];
                    allCWDiff = [allCWDiff PVADiff];
                end
            end

            for tPt = 1:length(CCWPer)-1
                if PVAStren_R(CCWPer(tPt)) > PVAThresh & PVAStren_G(CCWPer(tPt)) > PVAThresh
                    PVADiff = min(abs(PVAPlt_R(CCWPer(tPt))-PVAPlt_G(CCWPer(tPt))),...
                        2*pi-abs(PVAPlt_R(CCWPer(tPt))-PVAPlt_G(CCWPer(tPt))));
                    allCCWvR = [allCCWvR vR(CCWPer(tPt))];
                    allCCWDiff = [allCCWDiff PVADiff];
                end
            end
        end


        darkMeans = zeros(length(vRCents),1);
        dark25 = zeros(1,length(vRCents));
        dark75 = zeros(1,length(vRCents));
        vRDisc = discretize(allDarkvR,vREdges);
        for vRBin = 1:length(vRCents)
            rng = find(vRDisc == vRBin);
            darkMeans(vRBin) = mean(allDarkDiff(rng));
            dark25(vRBin) = prctile(allDarkDiff(rng),25);
            dark75(vRBin) = prctile(allDarkDiff(rng),75);
        end

        CLMeans = zeros(length(vRCents),1);
        CL25 = zeros(1,length(vRCents));
        CL75 = zeros(1,length(vRCents));
        vRDisc = discretize(allCLvR,vREdges);
        for vRBin = 1:length(vRCents)
            rng = find(vRDisc == vRBin);
            CLMeans(vRBin) = mean(allCLDiff(rng));
            CL25(vRBin) = prctile(allCLDiff(rng),25);
            CL75(vRBin) = prctile(allCLDiff(rng),75);
        end

        if condID == 1
            subplot(4,5,flyID);
            hold on;
            plot(vRCents,darkMeans,'k');
            patch('XData',[vRCents fliplr(vRCents)],'YData',[dark25 fliplr(dark75)],'EdgeColor','none','FaceColor','k','FaceAlpha',0.1);
            plot(vRCents,CLMeans,'b');
            patch('XData',[vRCents fliplr(vRCents)],'YData',[CL25 fliplr(CL75)],'EdgeColor','none','FaceColor','b','FaceAlpha',0.1);
            xlim([-pi pi]);
            ylim([0 3*pi/8]);
            xlabel('vR (rad/s)');
            xticks([-pi -pi/2 0 pi/2 pi]);
            ylabel('PVA difference');
            yticks([0 pi/8 pi/4 3*pi/8]);

            if flyID == EBFly
                subplot(4,5,16);
                hold on;
                plot(vRCents,darkMeans,'k');
                patch('XData',[vRCents fliplr(vRCents)],'YData',[dark25 fliplr(dark75)],'EdgeColor','none','FaceColor','k','FaceAlpha',0.1);
                plot(vRCents,CLMeans,'b');
                patch('XData',[vRCents fliplr(vRCents)],'YData',[CL25 fliplr(CL75)],'EdgeColor','none','FaceColor','b','FaceAlpha',0.1);
                xlim([-pi pi]);
                xticks([-pi -pi/2 0 pi/2 pi]);
                xlabel('vR (rad/s)');
                ylim([0 pi/4]);
                yticks(linspace(0,pi/4,5));
                ylabel('PVA difference');
            end

            subplot(4,5,17);
            hold on;
            scatter(1,mean(allDarkDiff),40,'k','filled');
            scatter(2,mean(allCLDiff),40,'b','filled');
            line([1 2],[mean(allDarkDiff) mean(allCLDiff)],'Color','k');
            scatter(3,mean(allCWDiff),40,'c','filled');
            line([2 3],[mean(allCLDiff) mean(allCWDiff)],'Color','k');
            scatter(4,mean(allCCWDiff),40,'c','filled');
            line([3 4],[mean(allCWDiff) mean(allCCWDiff)],'Color','k');
            alpha(0.5);
            xlim([0.5 4.5]);
            xticks([1 2 3 4]);
            xticklabels({'dark', 'CL', 'CW', 'CCW'});
            ylim([0 pi/4]);
            yticks(linspace(0,pi/4,5));
            ylabel('mean PVA difference');
        end

        visCond = 'Dark';
        allWids_R_Dark = [];
        allWids_G_Dark = [];

        vRIDs_Dark = discretize(WidAmpDat{condID}.(visCond){flyID}.allvR,vREdges);

        for ID = 1:max(vRIDs_Dark)
            pltRng = find(vRIDs_Dark == ID);

            wids2plt = WidAmpDat{condID}.(visCond){flyID}.allWid_R(pltRng);
            wids2plt = wids2plt(~isnan(wids2plt));
            allWids_R_Dark = vertcat(allWids_R_Dark,mean(wids2plt));

            wids2plt = WidAmpDat{condID}.(visCond){flyID}.allWid_G(pltRng);
            wids2plt = wids2plt(~isnan(wids2plt));
            allWids_G_Dark = vertcat(allWids_G_Dark,mean(wids2plt));
        end

        visCond = 'CL';
        allWids_R_CL = [];
        allWids_G_CL = [];

        vRIDs_CL = discretize(WidAmpDat{condID}.(visCond){flyID}.allvR,vREdges);

        for ID = 1:max(vRIDs_CL)
            pltRng = find(vRIDs_CL == ID);

            wids2plt = WidAmpDat{condID}.(visCond){flyID}.allWid_R(pltRng);
            wids2plt = wids2plt(~isnan(wids2plt));
            allWids_R_CL = vertcat(allWids_R_CL,mean(wids2plt));

            wids2plt = WidAmpDat{condID}.(visCond){flyID}.allWid_G(pltRng);
            wids2plt = wids2plt(~isnan(wids2plt));
            allWids_G_CL = vertcat(allWids_G_CL,mean(wids2plt));
        end

        if condID == 1
            subplot(4,5,5+flyID)
            hold on;
            plot(vRCents,allWids_R_Dark,'Color','m');
            plot(vRCents,allWids_R_CL,'Color','m');
            plot(vRCents,allWids_G_Dark,'Color','g');
            plot(vRCents,allWids_G_CL,'Color','g');
            xlim([0 pi]);
            xticks([0 pi/4 pi/2 3*pi/4 pi]);
            yticks([0 pi/4 pi/2 3*pi/4 pi]);
            ylim([0 pi]);
            title(cond{condID}.name);
            if condID == 1
                ylabel('mean(FWHM)');
            end
            xlabel('vR (rad/s)');
            set(gca,'FontSize',10);

            subplot(4,5,10+flyID);
            hold on;
            plot(vRCents,darkMeans./(0.5*(allWids_R_Dark+allWids_G_Dark)),'k');
            patch('XData',[vRCents fliplr(vRCents)],...
                'YData',[dark25./(0.5*(allWids_R_Dark'+allWids_G_Dark')) fliplr(dark75)./(0.5*(allWids_R_Dark'+allWids_G_Dark'))],...
                'EdgeColor','none','FaceColor','k','FaceAlpha',0.1);
            plot(vRCents,CLMeans./(0.5*(allWids_R_CL+allWids_G_CL)),'b');
            patch('XData',[vRCents fliplr(vRCents)],...
                'YData',[CL25./(0.5*(allWids_R_CL'+allWids_G_CL')) fliplr(CL75)./(0.5*(allWids_R_CL'+allWids_G_CL'))],...
                'EdgeColor','none','FaceColor','b','FaceAlpha',0.1);
            xlim([-pi pi]);
            ylim([0 0.5]);
            xlabel('vR (rad/s)');
            xticks([-pi -pi/2 0 pi/2 pi]);
            ylabel('PVA difference (% of bump width)');

            if flyID == EBFly
                subplot(4,5,19);
                hold on;
                plot(vRCents,darkMeans./(0.5*(allWids_R_Dark+allWids_G_Dark)),'k');
                patch('XData',[vRCents fliplr(vRCents)],...
                    'YData',[dark25./(0.5*(allWids_R_Dark'+allWids_G_Dark')) fliplr(dark75)./(0.5*(allWids_R_Dark'+allWids_G_Dark'))],...
                    'EdgeColor','none','FaceColor','k','FaceAlpha',0.1);
                plot(vRCents,CLMeans./(0.5*(allWids_R_CL+allWids_G_CL)),'b');
                patch('XData',[vRCents fliplr(vRCents)],...
                    'YData',[CL25./(0.5*(allWids_R_CL'+allWids_G_CL')) fliplr(CL75)./(0.5*(allWids_R_CL'+allWids_G_CL'))],...
                    'EdgeColor','none','FaceColor','b','FaceAlpha',0.1);
                xlim([-pi pi]);
                xticks([-pi -pi/2 0 pi/2 pi]);
                xlabel('vR (rad/s)');
                ylim([0 0.4]);
                ylabel('PVA difference (% of bump width)');
            end
        end

        visCond = 'Dark';
        darkWids_R = WidAmpDat{condID}.(visCond){flyID}.allWid_R;
        darkWids_R = darkWids_R(~isnan(darkWids_R));
        darkWids_G = WidAmpDat{condID}.(visCond){flyID}.allWid_G;
        darkWids_G = darkWids_G(~isnan(darkWids_G));

        visCond = 'CL';
        CLWids_R = WidAmpDat{condID}.(visCond){flyID}.allWid_R;
        CLWids_R = CLWids_R(~isnan(CLWids_R));
        CLWids_G = WidAmpDat{condID}.(visCond){flyID}.allWid_G;
        CLWids_G = CLWids_G(~isnan(CLWids_G));

        subplot(4,5,20);
        hold on;
        scatter(1,mean(allDarkDiff)/(0.5*mean(darkWids_R)+0.5*mean(darkWids_G)),40,'k','filled');
        scatter(2,mean(allCLDiff)/(0.5*mean(CLWids_R)+0.5*mean(CLWids_G)),40,'b','filled');
        if condID == 1
            pltCol = 'g';
        else
            pltCol = 'm';
        end
        line([1 2],[mean(allDarkDiff)/(0.5*mean(darkWids_R)+0.5*mean(darkWids_G)) mean(allCLDiff)/(0.5*mean(CLWids_R)+0.5*mean(CLWids_G))],'Color',pltCol);
        alpha(0.5);
        xlim([0.5 2.5]);
        xticks([1 2]);
        xticklabels({'dark', 'CL'});
        ylim([0 0.6]);
        ylabel('mean PVA difference');
    end
end