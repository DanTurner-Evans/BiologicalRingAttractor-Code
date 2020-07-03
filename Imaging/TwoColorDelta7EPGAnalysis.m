%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\TwoColor\';

%% Load the data
cond{1}.name = 'Green - E-PG, Red - D7';
cond{1}.dirs{1} = strcat(dataDir,'GreenEPGsRedDelta7s\20170717');

cond{2}.name = 'Red - E-PG, Green - D7';
cond{2}.dirs{1} = strcat(dataDir,'RedEPGsGreenDelta7s\20170706');
cond{2}.dirs{2} = strcat(dataDir,'RedEPGsGreenDelta7s\20170707');

cond = FlyDatLoad(2,cond);


%% Specify parameters to be used in multiple analyses
glomShift = 3; % How many glomeruli to shift the data for display
vRMin = 0; % min rotational velocity bin
vRMax = 720; % max rotational velocity bin
vRSpan = 60; % span of rotational velocities

% The glomeruli where the neurons of interest arborize
redSpan = [2:17];
greenSpan = [2:17];

% Color maps
blues = brewermap(64, 'Blues');
greens(:,1) = blues(:,1);
greens(:,2) = blues(:,3);
greens(:,3) = blues(:,2);
magentas = blues;
magentas(:,1) = blues(:,3);
magentas(:,2) = blues(:,1);
magentas(:,3) = blues(:,3);
colorAssign = linspace(0,0.5,64);

% Align the bumps
[LPB1.dark.dataR, RPB1.dark.dataR, LPB1.dark.dataG, RPB1.dark.dataG] = ...
    BumpAlignOnePeak(cond{1}.allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'Dark',glomShift);
[LPB2.dark.dataR, RPB2.dark.dataR, LPB2.dark.dataG, RPB2.dark.dataG] = ...
    BumpAlignOnePeak(cond{2}.allFlyData,vRMin,vRMax,vRSpan,redSpan(1:8),greenSpan(1:8),'Dark',glomShift);

% Find the number of velocity bins
vRBinNum = round((vRMax-vRMin)/vRSpan);

% Specify the number of ROIs and get the correspondin angles
num_ROIs = 8;
angsraw = (1:num_ROIs)*2*pi/num_ROIs-pi;
angsraw = angsraw';


%% Align the peaks and plot the normalized mean values as sorted by velocity - Figure 4 I 

% Step through the conditions
for condNow = 1:2
    
    if condNow == 1
        LPB = LPB1;
        RPB = RPB1;
    elseif condNow == 2
        LPB = LPB2;
        RPB = RPB2;
    end
    
    PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

    % Set the range of values to plot
    plotRange = 180/vRSpan;
    
    % Step through the flies
    for flyID = 1:cond{condNow}.numFlies

        % Start by looking at the periods where the fly is stopped
        LPBDataRAllStop = LPB.dark.dataR{flyID}.Stop;
        RPBDataRAllStop = RPB.dark.dataR{flyID}.Stop;
        LPBDataGAllStop = LPB.dark.dataG{flyID}.Stop;
        RPBDataGAllStop = RPB.dark.dataG{flyID}.Stop;

        % Find the mean signal and normalize it
        % For the right and left red channel data
        LRStop = mean(LPBDataRAllStop,2);
        LRStop = (LRStop-min(LRStop(redSpan(1:8))))./(max(LRStop)-min(LRStop(redSpan(1:8))));
        RRStop = mean(RPBDataRAllStop,2);
        RRStop = (RRStop-min(RRStop(greenSpan(1:8))))./(max(RRStop)-min(RRStop(greenSpan(1:8))));
        RStop = vertcat(LRStop,RRStop);

        % For the right and left green channel data
        LGStop = mean(LPBDataGAllStop,2);
        LGStop = (LGStop-min(LGStop(greenSpan(1:8))))./(max(LGStop)-min(LGStop(greenSpan(1:8))));
        RGStop = mean(RPBDataGAllStop,2);
        RGStop = (RGStop-min(RGStop(redSpan(1:8))))./(max(RGStop)-min(RGStop(redSpan(1:8))));
        GStop = vertcat(LGStop,RGStop);

        % Find the peak locations and the difference between peaks
        LRPk = find(LRStop == max(LRStop));
        LGPk = find(LGStop == max(LGStop));
        pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
        RRPk = find(RRStop == max(RRStop));
        RGPk = find(RGStop == max(RGStop));
        pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));

        % Plot the actiivty profiles
        subplot(2*cond{condNow}.numFlies,2*(plotRange+1),1+4*(plotRange+1)*(flyID-1))
        LBinR = mean(LPBDataRAllStop,2);
        RBinR = mean(RPBDataRAllStop,2);
        RBinR = circshift(RBinR,-1);
        actVals = vertcat(LBinR,RBinR)'-1;
        colorIm = zeros(2,18,3);
        for colorStep = 1:18
            colorNow = find(actVals(colorStep)<colorAssign)';
            colorNow(end+1) = 64;
            colorIm(1,colorStep,:) = magentas(colorNow(1),:);
        end
        LBinG =  mean(LPBDataGAllStop,2);
        RBinG = mean(RPBDataGAllStop,2);
        RBinG = circshift(RBinG,-1);
        actVals = vertcat(LBinG,RBinG)'-1;
        for colorStep = 1:18
            colorNow = find(actVals(colorStep)<colorAssign)';
            colorNow(end+1) = 64;
            colorIm(2,colorStep,:) = greens(colorNow(1),:);
        end
        image(colorIm);
        line([circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
            circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
            [0.5 1.5],'Color','k','LineWidth',2);
        line([circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
            circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
            [1.5 2.5],'Color','k','LineWidth',2);
%         text(2,0,num2str(pkDiffL));
%         text(11,0,num2str(pkDiffR));
        axis off;
        title(strcat('stopped, Fly ',num2str(flyID)));
        
        subplot(2*cond{condNow}.numFlies,2*(plotRange+1),2+4*(plotRange+1)*(flyID-1))
        RBinR = mean(RPBDataRAllStop,2);
        actVals = RBinR(2:9)-1;
        colorIm = zeros(2,8,3);
        for colorStep = 1:8
            colorNow = find(actVals(colorStep)<colorAssign)';
            colorNow(end+1) = 64;
            colorIm(1,colorStep,:) = magentas(colorNow(1),:);
        end
        RBinG = mean(RPBDataGAllStop,2);
        actVals = RBinG(2:9)-1;
        for colorStep = 1:8
            colorNow = find(actVals(colorStep)<colorAssign)';
            colorNow(end+1) = 64;
            colorIm(2,colorStep,:) = greens(colorNow(1),:);
        end
        image(colorIm);
        line([circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
            circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
            [0.5 1.5],'Color','k','LineWidth',2);
        line([circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
            circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
            [1.5 2.5],'Color','k','LineWidth',2);
%         text(2,0,num2str(pkDiffL));
%         text(11,0,num2str(pkDiffR));
        axis off;
        title(strcat('stopped, Fly ',num2str(flyID)));


        % Now, do the same as above for the CW turns
        for binID=1:plotRange
            LPBDataRAllCW = LPB.dark.dataR{flyID}.CW{binID};
            RPBDataRAllCW = RPB.dark.dataR{flyID}.CW{binID};
            LPBDataGAllCW = LPB.dark.dataG{flyID}.CW{binID};
            RPBDataGAllCW = RPB.dark.dataG{flyID}.CW{binID};

            if ~isempty(LPBDataRAllCW)
                LRCW = mean(LPBDataRAllCW,2);
                LRCW = (LRCW-min(LRCW(redSpan(1:8))))./(max(LRCW)-min(LRCW(redSpan(1:8))));
                RRCW = mean(RPBDataRAllCW,2);
                RRCW = (RRCW-min(RRCW(greenSpan(1:8))))./(max(RRCW)-min(RRCW(greenSpan(1:8))));
                RCW = vertcat(LRCW,RRCW);
                LGCW = mean(LPBDataGAllCW,2);
                LGCW = (LGCW-min(LGCW(greenSpan(1:8))))./(max(LGCW)-min(LGCW(greenSpan(1:8))));
                RGCW = mean(RPBDataGAllCW,2);
                RGCW = (RGCW-min(RGCW(redSpan(1:8))))./(max(RGCW)-min(RGCW(redSpan(1:8))));
                GCW = vertcat(LGCW,RGCW);

                LRPk = find(LRCW == max(LRCW));
                LGPk = find(LGCW == max(LGCW));
                pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
                RRPk = find(RRCW == max(RRCW));
                RGPk = find(RGCW == max(RGCW));
                pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));


                subplot(2*cond{condNow}.numFlies,2*(plotRange+1),1+2*binID+4*(plotRange+1)*(flyID-1))
                LBinR = mean(LPBDataRAllCW,2);
                actVals = LBinR(2:9)-1;
                colorIm = zeros(2,8,3);
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(1,colorStep,:) = magentas(colorNow(1),:);
                end
                LBinG = mean(LPBDataGAllCW,2);
                actVals = LBinG(2:9)-1;
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(2,colorStep,:) = greens(colorNow(1),:);
                end
                image(colorIm);
                line([circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [0.5 1.5],'Color','k','LineWidth',2);
                line([circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [1.5 2.5],'Color','k','LineWidth',2);
%                 text(2,0,num2str(pkDiffL));
%                 text(11,0,num2str(pkDiffR));
                axis off;
                title(strcat(num2str(vRSpan*(binID-1)),'-',num2str(vRSpan*binID),' deg/s CW'));
                
                subplot(2*cond{condNow}.numFlies,2*(plotRange+1),2+2*binID+4*(plotRange+1)*(flyID-1))
                RBinR = mean(RPBDataRAllCW,2);
                actVals = RBinR(2:9)-1;
                colorIm = zeros(2,8,3);
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(1,colorStep,:) = magentas(colorNow(1),:);
                end
                RBinG = mean(RPBDataGAllCW,2);
                actVals = RBinG(2:9)-1;
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(2,colorStep,:) = greens(colorNow(1),:);
                end
                image(colorIm);
                line([circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [0.5 1.5],'Color','k','LineWidth',2);
                line([circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [1.5 2.5],'Color','k','LineWidth',2);
%                 text(2,0,num2str(pkDiffL));
%                 text(11,0,num2str(pkDiffR));
                axis off;
                title(strcat(num2str(vRSpan*(binID-1)),'-',num2str(vRSpan*binID),' deg/s CW'));
            end

            % Now, do as above for the CCW turns
            LPBDataRAllCCW = LPB.dark.dataR{flyID}.CCW{binID};
            RPBDataRAllCCW = RPB.dark.dataR{flyID}.CCW{binID};
            LPBDataGAllCCW = LPB.dark.dataG{flyID}.CCW{binID};
            RPBDataGAllCCW = RPB.dark.dataG{flyID}.CCW{binID};

            if ~isempty(LPBDataRAllCCW)
                LRCCW = mean(LPBDataRAllCCW,2);
                LRCCW = (LRCCW-min(LRCCW(redSpan(1:8))))./(max(LRCCW)-min(LRCCW(redSpan(1:8))));
                RRCCW = mean(RPBDataRAllCCW,2);
                RRCCW = (RRCCW-min(RRCCW(greenSpan(1:8))))./(max(RRCCW)-min(RRCCW(greenSpan(1:8))));
                RCCW = vertcat(LRCCW,RRCCW);
                LGCCW = mean(LPBDataGAllCCW,2);
                LGCCW = (LGCCW-min(LGCCW(greenSpan(1:8))))./(max(LGCCW)-min(LGCCW(greenSpan(1:8))));
                RGCCW = mean(RPBDataGAllCCW,2);
                RGCCW = (RGCCW-min(RGCCW(redSpan(1:8))))./(max(RGCCW)-min(RGCCW(redSpan(1:8))));
                GCCW = vertcat(LGCCW,RGCCW);

                LRPk = find(LRCCW == max(LRCCW));
                LGPk = find(LGCCW == max(LGCCW));
                pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
                RRPk = find(RRCCW == max(RRCCW));
                RGPk = find(RGCCW == max(RGCCW));
                pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));

                subplot(2*cond{condNow}.numFlies,2*(plotRange+1),1+2*binID+4*(plotRange+1)*(flyID-1)+2*(plotRange+1))
                LBinR = mean(LPBDataRAllCCW,2);
                actVals = LBinR(2:9)-1;
                colorIm = zeros(2,8,3);
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(1,colorStep,:) = magentas(colorNow(1),:);
                end
                LBinG = mean(LPBDataGAllCCW,2);
                actVals = LBinG(2:9)-1;
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(2,colorStep,:) = greens(colorNow(1),:);
                end
                image(colorIm);
                line([circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,LBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [0.5 1.5],'Color','k','LineWidth',2);
                line([circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,LBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [1.5 2.5],'Color','k','LineWidth',2);
%                 text(2,0,num2str(pkDiffL));
%                 text(11,0,num2str(pkDiffR));
                axis off;
                title(strcat(num2str(vRSpan*(binID-1)),'-',num2str(vRSpan*binID),' deg/s CCW'));
                
                subplot(2*cond{condNow}.numFlies,2*(plotRange+1),2+2*binID+4*(plotRange+1)*(flyID-1)+2*(plotRange+1))
                RBinR = mean(RPBDataRAllCCW,2);
                actVals = RBinR(2:9)-1;
                colorIm = zeros(2,8,3);
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(1,colorStep,:) = magentas(colorNow(1),:);
                end
                RBinG = mean(RPBDataGAllCCW,2);
                actVals = RBinG(2:9)-1;
                for colorStep = 1:8
                    colorNow = find(actVals(colorStep)<colorAssign)';
                    colorNow(end+1) = 64;
                    colorIm(2,colorStep,:) = greens(colorNow(1),:);
                end
                image(colorIm);
                line([circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,RBinR(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [0.5 1.5],'Color','k','LineWidth',2);
                line([circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2,...
                    circ_mean(angsraw,RBinG(2:9))*num_ROIs/(2*pi)+num_ROIs/2],...
                    [1.5 2.5],'Color','k','LineWidth',2);
%                 text(2,0,num2str(pkDiffL));
%                 text(11,0,num2str(pkDiffR));
                axis off;
                title(strcat(num2str(vRSpan*(binID-1)),'-',num2str(vRSpan*binID),' deg/s CCW'));
            end
        end
    end

    figure(PBFig);
    [ax,h] = suplabel(cond{condNow}.name,'t',[0.1 0.1 0.85 0.85]);
end

%% Make colorbars - Figure 4 I

blues = brewermap(64, 'Blues');
greens(:,1) = blues(:,1);
greens(:,2) = blues(:,3);
greens(:,3) = blues(:,2);
magentas = blues;
magentas(:,1) = blues(:,3);
magentas(:,2) = blues(:,1);
magentas(:,3) = blues(:,3);

cBar = figure;
subplot(2,10,1);
hold on;
for colNow = 1:length(colorAssign)
    patch('XData',[0 10 10 0],'YData',[colNow colNow colNow+1 colNow+1],'FaceColor',magentas(colNow,:),'EdgeColor','none');
end
text(11,0,'0');
text(11,length(colorAssign),num2str(max(colorAssign)));
axis off;

subplot(2,10,2);
hold on;
for colNow = 1:length(colorAssign)
    patch('XData',[0 10 10 0],'YData',[colNow colNow colNow+1 colNow+1],'FaceColor',greens(colNow,:),'EdgeColor','none');
end
text(11,0,'0');
text(11,length(colorAssign),num2str(max(colorAssign)));
axis off;


%% Run linear regressions with activity peak, rotational velocity, and forward velocity

% Specify the glomeruli to look at
PBRange = [2:9];

% Pick the number of lead/lag time points to consider
numtPts = 51;

% Step through the conditions
for condNow = 1:2
    
    % Initialize a figure
    figure;
    
    % Step through the flies
    for flyID = 1:cond{condNow}.numFlies
        
        % Initialize arrays to hold the regression coefficients
        lrCoefR = zeros(numtPts,3);
        lrCoefG = zeros(numtPts,3);
        
        % Initialize arrays to hold the behavioral and imaging data
        vRot = [];
        vF = [];
        RMaxVals = [];
        GMaxVals = [];
        
        % Step through the trials, concatentating the data together
        for trial = 1:length(cond{condNow}.allFlyData{flyID}.Dark)
           vRot = vertcat(vRot, ...
               cond{condNow}.allFlyData{flyID}.Dark{trial}.positionDatMatch.vRot);
           vF = vertcat(vF, ...
               cond{condNow}.allFlyData{flyID}.Dark{trial}.positionDatMatch.vF);
           RMaxVals = vertcat(RMaxVals,...
               max(cond{condNow}.allFlyData{flyID}.Dark{trial}.RROIaveMax(PBRange,2:end))');
           GMaxVals = vertcat(GMaxVals,...
               max(cond{condNow}.allFlyData{flyID}.Dark{trial}.GROIaveMax(PBRange,2:end))');
        end
        
        % Sort the velocity
        vCW = zeros(size(vRot));
        vCCW = zeros(size(vRot));
        vPos = find(vRot > 0);
        vNeg = find(vRot < 0);
        vCW(vPos) = vRot(vPos);
        vCCW(vNeg) = abs(vRot(vNeg));
        vPred = horzcat(zscore(vCW),zscore(vCCW),zscore(vF)); 
        
        % Calculate the regression coefficients
        for tStep = 1:numtPts
            lrCoefR(tStep,:) = ...
                regress(...
                zscore(RMaxVals(tStep:end-numtPts+tStep)),...
                vPred(floor(numtPts/2):end-ceil(numtPts/2),:));
            lrCoefG(tStep,:) = ...
                regress(...
                zscore(GMaxVals(tStep:end-numtPts+tStep)),...
                vPred(floor(numtPts/2):end-ceil(numtPts/2),:));
        end
        
        % Plot!
        subplot(2,cond{condNow}.numFlies,flyID);
        hold on;
        plot(lrCoefR(:,1),'b');
        plot(lrCoefR(:,2),'r');
        plot(lrCoefR(:,3),'m');
        
        subplot(2,cond{condNow}.numFlies,flyID+cond{condNow}.numFlies);
        hold on;
        plot(lrCoefG(:,1),'b');
        plot(lrCoefG(:,2),'r');
        plot(lrCoefG(:,3),'m');
        
    end
    
end

%% Look at relationship between red and green peak values 

% Specify the glomeruli to look at
PBRange = [2:9];

% Step through the conditions
for condNow = 1:2
    
    % Initialize the figure
    figure;
    
    % Step through the flies
    for flyID = 1:cond{condNow}.numFlies
        
        % Look at the maximum values
        RMaxVals = [];
        GMaxVals = [];
        
        for trial = 1:length(cond{condNow}.allFlyData{flyID}.Dark)
           RMaxVals = vertcat(RMaxVals,...
               max(cond{condNow}.allFlyData{flyID}.Dark{trial}.RROIaveMax(PBRange,:))');
           GMaxVals = vertcat(GMaxVals,...
               max(cond{condNow}.allFlyData{flyID}.Dark{trial}.GROIaveMax(PBRange,:))');
        end
        
        subplot(1,cond{condNow}.numFlies,flyID);
        scatter(RMaxVals,GMaxVals,20,'filled')
        alpha(0.1);
        
    end
    
end


%% Load an example tiff - Figure 4 G
pathName = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\TwoColor\GreenEPGsRedDelta7s\20170717\';
tiffFilename = 'Fly2_3-4day_6fx60D05_jRGCx55G08_Dark_00001.tif';

behavFilename = 'Fly2_3-4day_6fx60D05_jRGCx55G08_Dark_01.TXT';
positionDat = VRDatLoad(behavFilename,pathName,0);

% Load the imaging stack and get the behavioral data
[RstackMaxInt, GstackMaxInt, RstackMean, GstackMean] = ImDatLoadBigtiff2Color(tiffFilename,pathName);

%% Preprocess the imaging stack - Figure 4 G
% Subtract the background
RMean = mean(RstackMaxInt,3);
GMean = mean(GstackMaxInt,3);
for frame = 1:size(RstackMaxInt,3)
    RstackSub(:,:,frame) = RstackMaxInt(:,:,frame) - RMean;
    GstackSub(:,:,frame) = GstackMaxInt(:,:,frame) - GMean;
end

% Gaussian filter the stacks
gaussianSize = [7 7];
gaussianSigma = 3;
Gxy = fspecial('gaussian',gaussianSize,gaussianSigma);
RstackXYfiltMax = double(zeros(size(RstackMaxInt)));
GstackXYfiltMax = double(zeros(size(GstackMaxInt)));

h = waitbar(0.0,'Gaussian filtering stack...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Gaussian filtering TIFF stack...');
for i = 1:size(RstackMaxInt,3)
    if mod(i,100)==0
        waitbar(i/length(RstackMaxInt),h,['Filtering frame# ' num2str(i) ' out of ' num2str(length(RstackMaxInt))]);
    end
    RstackXYfiltMax(:,:,i) = imfilter(RstackMaxInt(:,:,i),Gxy,'replicate','conv');
    GstackXYfiltMax(:,:,i) = imfilter(GstackMaxInt(:,:,i),Gxy,'replicate','conv');
end
delete(h);

%% Show stills of the bump moving in both populations across three time points and the max. int. plot over time - Figure 4 G

% Find the min and max of the imaging stack and the imaging rate
RmaxCa = max(max(max(RstackXYfiltMax)));
GmaxCa = max(max(max(GstackXYfiltMax)));
RminCa = RmaxCa/4;% min(min(min(RstackXYfiltMax)));
GminCa = GmaxCa/6;% min(min(min(GstackXYfiltMax)));
imFreq = 10000/mean(diff(positionDat.tFrameGrab))/num_planes; 

% Pick example time points to plot
exPoint1 = round(imFreq*38.7);
exPoint2 = round(imFreq*42.7);
exPoint3 = round(imFreq*51.7);

% Initialize the figure
TwoPop = figure('Position',[50 50 800 1000]);

% Find the mean red activity
RMean = mean(RstackMaxInt(60:160,:,:),3);
RMean = RMean./max(max(RMean));

% Plot the mean red activity
subplot(6,4,1);
imshow(RMean);

% Find the mean green activity
GMean = mean(GstackMaxInt(60:160,:,:),3);
GMean = GMean./max(max(GMean));

% Plot the mean green activity
subplot(6,4,5);
imshow(GMean);

% Create the colored images from the two stacks for the various time points
% and plot them
curRFrame = 2*(RstackXYfiltMax(60:160,:,exPoint1)-...
    RminCa)./(RmaxCa-RminCa);
curGFrame = 6*(GstackXYfiltMax(60:160,:,exPoint1)-...
    GminCa)./(GmaxCa-GminCa);

RIm = zeros([size(curRFrame) 3]);
RIm(:,:,1) = curRFrame;
RIm(:,:,3) = curRFrame;
GIm = zeros([size(curGFrame) 3]);
GIm(:,:,2) = curGFrame;

subplot(6,4,2);
cla;
imshow(RIm);
axis off;

subplot(6,4,6);
cla;
imshow(GIm);
axis off;

subplot(6,4,10);
cla;
imshow(GIm+RIm);
axis off;

colorBarIm = zeros(2,50,3);
colorBarIm(1,:,1) = linspace(1,50,50)./50;
colorBarIm(1,:,3) = linspace(1,50,50)./50;
colorBarIm(2,:,2) = linspace(1,50,50)./50;

subplot(6,4,9);
image(1:50,1:2,colorBarIm)
set(gca,'XTick',[1 50],'XTickLabels',...
    {num2str(GminCa) num2str((GmaxCa-GminCa)/6+GminCa)});
text(-5,0.25,num2str(RminCa));
text(45,0.25,num2str((RmaxCa-RminCa)/2+RminCa));

curRFrame = 2*(RstackXYfiltMax(60:160,:,exPoint2)-...
    RminCa)./(RmaxCa-RminCa);
curGFrame = 6*(GstackXYfiltMax(60:160,:,exPoint2)-...
    GminCa)./(GmaxCa-GminCa);

RIm = zeros([size(curRFrame) 3]);
RIm(:,:,1) = curRFrame;
RIm(:,:,3) = curRFrame;
GIm = zeros([size(curGFrame) 3]);
GIm(:,:,2) = curGFrame;

subplot(6,4,3);
cla;
imshow(RIm);
axis off;

subplot(6,4,7);
cla;
imshow(GIm);
axis off;

subplot(6,4,11);
cla;
imshow(GIm+RIm);
axis off;

curRFrame = 2*(RstackXYfiltMax(60:160,:,exPoint3)-...
    RminCa)./(RmaxCa-RminCa);
curGFrame = 6*(GstackXYfiltMax(60:160,:,exPoint3)-...
    GminCa)./(GmaxCa-GminCa);

RIm = zeros([size(curRFrame) 3]);
RIm(:,:,1) = curRFrame;
RIm(:,:,3) = curRFrame;
GIm = zeros([size(curGFrame) 3]);
GIm(:,:,2) = curGFrame;

subplot(6,4,4);
cla;
imshow(RIm);
axis off;

subplot(6,4,8);
cla;
imshow(GIm);
axis off;

subplot(6,4,12);
cla;
imshow(GIm+RIm);
axis off;

fmStart = round(imFreq*35);
fmEnd = round(imFreq*55);

tAll = cond{1}.allFlyData{2}.Dark{1}.positionDatMatch.OffsetRotMatch(fmStart:fmEnd,1);
t1 = cond{1}.allFlyData{2}.Dark{1}.positionDatMatch.OffsetRotMatch(exPoint1,1);
t2 = cond{1}.allFlyData{2}.Dark{1}.positionDatMatch.OffsetRotMatch(exPoint2,1);
t3 = cond{1}.allFlyData{2}.Dark{1}.positionDatMatch.OffsetRotMatch(exPoint3,1);
RData = cond{1}.allFlyData{2}.Dark{1}.RROIaveMax(2:17,fmStart:fmEnd);
GData = cond{1}.allFlyData{2}.Dark{1}.GROIaveMax(2:17,fmStart:fmEnd);

RIm = zeros([size(RData) 3]);
GIm = zeros([size(GData) 3]);

RMax = max(max(RData));
RMin = min(min(RData));
GMax = max(max(GData));
GMin = min(min(GData));

RIm(:,:,1) = (RData-RMin)./(RMax-RMin);
RIm(:,:,3) = (RData-RMin)./(RMax-RMin);
GIm(:,:,2) = (GData-GMin)./(GMax-GMin);

subplot(6,4,[14:16]);
image(tAll,2:17,RIm);
set(gca,'XTick',[]);
ylabel('PB glomerulus');
title('\Delta7 neurons');
line([t1 t1],[1 18],'Color','w','LineStyle','--');
line([t2 t2],[1 18],'Color','w','LineStyle','--');
line([t3 t3],[1 18],'Color','w','LineStyle','--');
line([tAll(1) tAll(end)],[9.5 9.5],'Color','w');
set(gca,'XColor','w');

subplot(6,4,[18:20]);
image(tAll,2:17,GIm);
ylabel('PB glomerulus');
title('compass neurons');
line([t1 t1],[1 18],'Color','w','LineStyle','--');
line([t2 t2],[1 18],'Color','w','LineStyle','--');
line([t3 t3],[1 18],'Color','w','LineStyle','--');
line([tAll(1) tAll(end)],[9.5 9.5],'Color','w');


subplot(6,4,[22:24]);
image(tAll,2:17,RIm+GIm);
xlabel('time (s)')
ylabel('PB glomerulus');
line([t1 t1],[1 18],'Color','w','LineStyle','--');
line([t2 t2],[1 18],'Color','w','LineStyle','--');
line([t3 t3],[1 18],'Color','w','LineStyle','--');
line([tAll(1) tAll(end)],[9.5 9.5],'Color','w');


%% Plot the offsets as a function of velocity - Figure 4 J, Figure S6 C

% Initialize the figure
PBFig = figure('units','normalized','outerposition',[0 0 1 1]);

% Step through the conditions
for condNow = 1:2
    
    % Specify data and plotting parameters for each condition
    if condNow == 1
        LPB = LPB1;
        RPB = RPB1;
        scatColL = [0.2 0.4 0.8];
        scatColR = [0.4 0.2 0.8];
    elseif condNow == 2
        LPB = LPB2;
        RPB = RPB2;
        scatColL = [0.8 0.4 0.2];
        scatColR = [0.8 0.2 0.4];
    end

    % Specify the plot range
    plotRange = 180/vRSpan;
    
    % Step through the flies
    for flyID = 1:cond{condNow}.numFlies

        % Start by looking at the periods where the fly is stopped
        LPBDataRAllStop = LPB.dark.dataR{flyID}.Stop;
        RPBDataRAllStop = RPB.dark.dataR{flyID}.Stop;
        LPBDataGAllStop = LPB.dark.dataG{flyID}.Stop;
        RPBDataGAllStop = RPB.dark.dataG{flyID}.Stop;

        % Find the mean signal and normalize it
        % For the right and left red channel data
        LRStop = mean(LPBDataRAllStop,2);
        LRStop = (LRStop-min(LRStop(redSpan(1:8))))./(max(LRStop)-min(LRStop(redSpan(1:8))));
        RRStop = mean(RPBDataRAllStop,2);
        RRStop = (RRStop-min(RRStop(greenSpan(1:8))))./(max(RRStop)-min(RRStop(greenSpan(1:8))));
        RStop = vertcat(LRStop,RRStop);

        % For the right and left green channel data
        LGStop = mean(LPBDataGAllStop,2);
        LGStop = (LGStop-min(LGStop(greenSpan(1:8))))./(max(LGStop)-min(LGStop(greenSpan(1:8))));
        RGStop = median(RPBDataGAllStop,2);
        RGStop = (RGStop-min(RGStop(redSpan(1:8))))./(max(RGStop)-min(RGStop(redSpan(1:8))));
        GStop = vertcat(LGStop,RGStop);

        % Find the peak locations and the difference between peaks
        LRPk = find(LRStop == max(LRStop));
        LGPk = find(LGStop == max(LGStop));
        pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
        RRPk = find(RRStop == max(RRStop));
        RGPk = find(RGStop == max(RGStop));
        pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
        
        % Find the PVA difference between peaks
        PVADiffL = abs(circ_mean(angsraw,LGStop(2:9))-circ_mean(angsraw,LRStop(2:9)))...
            *num_ROIs/(2*pi);
        PVADiffL = min(PVADiffL,8-PVADiffL);
        PVADiffR = abs(circ_mean(angsraw,RGStop(2:9))-circ_mean(angsraw,RRStop(2:9)))...
            *num_ROIs/(2*pi);
        PVADiffR = min(PVADiffR,8-PVADiffR);

        % Plot the max difference and PVA difference between peaks
        subplot(2,2,1);
        hold on;
        scatter(0,pkDiffL,20,scatColL,'filled');
        scatter(0,pkDiffR,20,scatColR,'filled');
        xlabel('vR (deg/s)');
        ylabel('peak diff (# of glom.)');
        ylim([0 5]);
        set(gca,'XTick',[-180:30:180])
        
        subplot(2,2,3);
        hold on;
        scatter(0,PVADiffL,20,scatColL,'filled');
        scatter(0,PVADiffR,20,scatColR,'filled');
        xlabel('vR (deg/s)');
        ylabel('PVA diff (# of glom.)');
        ylim([0 5]);
        set(gca,'XTick',[-180:30:180])

        % As above, not for clockwise turns
        for binID=1:plotRange
            LPBDataRAllCW = LPB.dark.dataR{flyID}.CW{binID};
            RPBDataRAllCW = RPB.dark.dataR{flyID}.CW{binID};
            LPBDataGAllCW = LPB.dark.dataG{flyID}.CW{binID};
            RPBDataGAllCW = RPB.dark.dataG{flyID}.CW{binID};

            if ~isempty(LPBDataRAllCW)
                LRCW = mean(LPBDataRAllCW,2);
                LRCW = (LRCW-min(LRCW(redSpan(1:8))))./(max(LRCW)-min(LRCW(redSpan(1:8))));
                RRCW = mean(RPBDataRAllCW,2);
                RRCW = (RRCW-min(RRCW(greenSpan(1:8))))./(max(RRCW)-min(RRCW(greenSpan(1:8))));
                RCW = vertcat(LRCW,RRCW);
                LGCW = mean(LPBDataGAllCW,2);
                LGCW = (LGCW-min(LGCW(greenSpan(1:8))))./(max(LGCW)-min(LGCW(greenSpan(1:8))));
                RGCW = mean(RPBDataGAllCW,2);
                RGCW = (RGCW-min(RGCW(redSpan(1:8))))./(max(RGCW)-min(RGCW(redSpan(1:8))));
                GCW = vertcat(LGCW,RGCW);

                LRPk = find(LRCW == max(LRCW));
                LGPk = find(LGCW == max(LGCW));
                pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
                RRPk = find(RRCW == max(RRCW));
                RGPk = find(RGCW == max(RGCW));
                pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
                
                PVADiffL = abs(circ_mean(angsraw,LGCW(2:9))-circ_mean(angsraw,LRCW(2:9)))...
                    *num_ROIs/(2*pi);
                PVADiffL = min(PVADiffL,8-PVADiffL);
                PVADiffR = abs(circ_mean(angsraw,RGCW(2:9))-circ_mean(angsraw,RRCW(2:9)))...
                    *num_ROIs/(2*pi);
                PVADiffR = min(PVADiffR,8-PVADiffR);
                
                subplot(2,2,1);
                scatter(vRSpan*binID,pkDiffL,20,scatColL,'filled');
                scatter(vRSpan*binID,pkDiffR,20,scatColR,'filled');
                
                subplot(2,2,3);
                scatter(vRSpan*binID,PVADiffL,20,scatColL,'filled');
                scatter(vRSpan*binID,PVADiffR,20,scatColR,'filled');

            end

            % As above, not for counterclockwise turns
            LPBDataRAllCCW = LPB.dark.dataR{flyID}.CCW{binID};
            RPBDataRAllCCW = RPB.dark.dataR{flyID}.CCW{binID};
            LPBDataGAllCCW = LPB.dark.dataG{flyID}.CCW{binID};
            RPBDataGAllCCW = RPB.dark.dataG{flyID}.CCW{binID};

            if ~isempty(LPBDataRAllCCW)
                LRCCW = mean(LPBDataRAllCCW,2);
                LRCCW = (LRCCW-min(LRCCW(redSpan(1:8))))./(max(LRCCW)-min(LRCCW(redSpan(1:8))));
                RRCCW = mean(RPBDataRAllCCW,2);
                RRCCW = (RRCCW-min(RRCCW(greenSpan(1:8))))./(max(RRCCW)-min(RRCCW(greenSpan(1:8))));
                RCCW = vertcat(LRCCW,RRCCW);
                LGCCW = mean(LPBDataGAllCCW,2);
                LGCCW = (LGCCW-min(LGCCW(greenSpan(1:8))))./(max(LGCCW)-min(LGCCW(greenSpan(1:8))));
                RGCCW = mean(RPBDataGAllCCW,2);
                RGCCW = (RGCCW-min(RGCCW(redSpan(1:8))))./(max(RGCCW)-min(RGCCW(redSpan(1:8))));
                GCCW = vertcat(LGCCW,RGCCW);

                LRPk = find(LRCCW == max(LRCCW));
                LGPk = find(LGCCW == max(LGCCW));
                pkDiffL = min(abs(LGPk-LRPk),9-abs(LGPk-LRPk));
                RRPk = find(RRCCW == max(RRCCW));
                RGPk = find(RGCCW == max(RGCCW));
                pkDiffR = min(abs(RRPk-RGPk),9-abs(RRPk-RGPk));
                
                PVADiffL = abs(circ_mean(angsraw,LGCCW(2:9))-circ_mean(angsraw,LRCCW(2:9)))...
                    *num_ROIs/(2*pi);
                PVADiffL = min(PVADiffL,8-PVADiffL);
                PVADiffR = abs(circ_mean(angsraw,RGCCW(2:9))-circ_mean(angsraw,RRCCW(2:9)))...
                    *num_ROIs/(2*pi);
                PVADiffR = min(PVADiffR,8-PVADiffR);

                subplot(2,2,1);
                scatter(-vRSpan*binID,pkDiffL,20,scatColL,'filled');
                scatter(-vRSpan*binID,pkDiffR,20,scatColR,'filled');
                
                subplot(2,2,3);
                scatter(-vRSpan*binID,PVADiffL,20,scatColL,'filled');
                scatter(-vRSpan*binID,PVADiffR,20,scatColR,'filled');
            end
        end
    end
end

subplot(2,2,1);
alpha(0.5);
subplot(2,2,3);
alpha(0.5);

PBFig.Renderer='Painters';

%% Plot the activity over time for an example - Figure 4 H

% Set the PVA threshold
PVAThresh = 0.075;

% Choose an example fly
condID = 2;
flyID = 2;
trialID = 1;

% Set the DF/F ranges for the plot
minR = 0;
maxR = 0.4;
minG = 0;
maxG = 1.25;

% Set the savitzky-golay filter parameters
sgolayOrder = 3;
sgolayFrames = 11;

% Make the figure
TwoColPlt = figure('units','normalized','outerposition',[0 0 1 1]);

% Load the data
datNow = cond{condID}.allFlyData{flyID}.Dark{trialID};

% Pull out the behavior
tPts = datNow.positionDatMatch.OffsetRotMatch(:,1);
tPts = tPts-tPts(1);
stripePos = datNow.positionDatMatch.OffsetRotMatch(:,2);
heading = pi/360*datNow.positionDatMatch.PosRotMatch;

% Pull out the ca activity
DF_R = datNow.RROIaveMax-1;
DF_R = sgolayfilt(DF_R,sgolayOrder,sgolayFrames,[],2);
DF_G = datNow.GROIaveMax-1;
DF_G = sgolayfilt(DF_G,sgolayOrder,sgolayFrames,[],2);

% Remove the glomeruli where the neurons don't arborize and calculate the
% PVA
if condID == 1
    [angs, PVAPlt_R, PVAStren_R] = PVA(DF_R-min(min(DF_R)));
    PVAPlt_R(find(PVAStren_R<PVAThresh)) = NaN;
    [angs, PVAPlt_G, PVAStren_G] = PVA(DF_G-min(min(DF_G)));
    PVAPlt_G(find(PVAStren_G<PVAThresh)) = NaN;
else
    [angs, PVAPlt_R_L, PVAStren_R_L] = PVA(DF_R(1:9,:)-min(min(DF_R(1:9,:))));
    PVAPlt_R_L(find(PVAStren_R_L<PVAThresh)) = NaN;
    [angs, PVAPlt_G_L, PVAStren_G_L] = PVA(DF_G(1:9,:)-min(min(DF_G(1:9,:))));
    PVAPlt_G_L(find(PVAStren_G_L<PVAThresh)) = NaN;

    [angs, PVAPlt_R_R, PVAStren_R_R] = PVA(DF_R(10:18,:)-min(min(DF_R(10:18,:))));
    PVAPlt_R_R(find(PVAStren_R_R<PVAThresh)) = NaN;
    [angs, PVAPlt_G_R, PVAStren_G_R] = PVA(DF_G(10:18,:)-min(min(DF_G(10:18,:))));
    PVAPlt_R_L(find(PVAStren_R_L<PVAThresh)) = NaN;
end

% Remove jumping artifacts from the plot
for tPt = 2:length(tPts)
    if condID == 1
        if abs(PVAPlt_R(tPt)-PVAPlt_R(tPt-1))>pi
            PVAPlt_R(tPt-1) = NaN;
        end
        if abs(PVAPlt_G(tPt)-PVAPlt_G(tPt-1))>pi
            PVAPlt_G(tPt-1) = NaN;
        end
    elseif condID == 2
        if abs(PVAPlt_R_L(tPt)-PVAPlt_R_L(tPt-1))>pi
            PVAPlt_R_L(tPt-1) = NaN;
        end
        if abs(PVAPlt_G_L(tPt)-PVAPlt_G_L(tPt-1))>pi
            PVAPlt_G_L(tPt-1) = NaN;
        end
        if abs(PVAPlt_R_R(tPt)-PVAPlt_R_R(tPt-1))>pi
            PVAPlt_R_R(tPt-1) = NaN;
        end
        if abs(PVAPlt_G_R(tPt)-PVAPlt_G_R(tPt-1))>pi
            PVAPlt_G_R(tPt-1) = NaN;
        end
    end
end

% Normalize the data according to the plot bounds
RIm = zeros([size(DF_R) 3]);
RIm(:,:,1) = (DF_R-minR)./...
    (maxR-minR);
RIm(:,:,3) = (DF_R-minR)./...
    (maxR-minR);

GIm = zeros([size(DF_G) 3]);
GIm(:,:,2) = (DF_G-minG)./...
    (maxG-minG);

% Plot it!
subplot(2,1,2*(trialID-1)+1);
overlayIm = RIm+GIm;
if condID == 1
    image(tPts,angs,flipud(overlayIm));
else
    image(tPts,[1:18],overlayIm);
end
if trialID == 1
    title(strcat(cond{condID}.name,'-fly #',num2str(flyID)));
end
ylim([1.5 17.5]);

subplot(2,1,2*(trialID-1)+2);
hold on;
plot(tPts,stripePos,'k');
if condID == 1
    plot(tPts,PVAPlt_R,'r');
    plot(tPts,PVAPlt_G,'g');
else
    plot(tPts,-PVAPlt_R_L,'r');
    plot(tPts,-PVAPlt_G_L,'g');
    plot(tPts,-PVAPlt_R_R,'r');
    plot(tPts,-PVAPlt_G_R,'g');
end
xlim([tPts(1) tPts(end)]);
ylim([-pi pi]);

%% Find the mean and sd of the offset across flies

% Initialize an array to hold the offsets
allOffsets = []

% Step through the conditions
for condNow = 1:2
    
    if condNow == 1
        LPB = LPB1;
        RPB = RPB1;
    elseif condNow == 2
        LPB = LPB2;
        RPB = RPB2;
    end
    
    % Step through the flies
    for flyID = 1:cond{condNow}.numFlies

        LPBDataRAll = [LPB.dark.dataR{flyID}.Stop LPB.dark.dataR{flyID}.CW{1} LPB.dark.dataR{flyID}.CCW{1}];
        RPBDataRAll = [RPB.dark.dataR{flyID}.Stop RPB.dark.dataR{flyID}.CW{1} RPB.dark.dataR{flyID}.CCW{1}];
        LPBDataGAll = [LPB.dark.dataG{flyID}.Stop  LPB.dark.dataG{flyID}.CW{1} LPB.dark.dataG{flyID}.CCW{1}];
        RPBDataGAll = [RPB.dark.dataG{flyID}.Stop  RPB.dark.dataG{flyID}.CW{1} RPB.dark.dataG{flyID}.CCW{1}];

        LR = mean(LPBDataRAll,2);
        RR = mean(RPBDataRAll,2);

        LG = mean(LPBDataGAll,2);
        RG = median(RPBDataGAll,2);
        
        PVADiffL = abs(circ_mean(angsraw,LG(2:9))-circ_mean(angsraw,LR(2:9)))...
            *num_ROIs/(2*pi);
        PVADiffL = min(PVADiffL,8-PVADiffL);
        allOffsets = [allOffsets PVADiffL];
        
        PVADiffR = abs(circ_mean(angsraw,RG(2:9))-circ_mean(angsraw,RR(2:9)))...
            *num_ROIs/(2*pi);
        PVADiffR = min(PVADiffR,8-PVADiffR);
        allOffsets = [allOffsets PVADiffR];
    end
end

mean(allOffsets)
std(allOffsets)