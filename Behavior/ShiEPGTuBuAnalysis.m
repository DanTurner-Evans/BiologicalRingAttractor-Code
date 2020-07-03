%% Specify the data directory
% Replace this path with your own path to the data.
dataDir = 'C:\Users\turnerevansd\Documents\TheNeuroanatomicalUltrastructureAndFunctionOfABiologicalRingAttractor\Data\Behavior\';

%% Choose one of the following datasets to load
%% Specify the directories - AOTU
numDirs = 6;
saveName = '76B06';
allPathname{1} = strcat(dataDir,'TuBu\20170818\');
allPathname{2} = strcat(dataDir,'TuBu\20170819\');
allPathname{3} = strcat(dataDir,'TuBu\20170904\');
allPathname{4} = strcat(dataDir,'TuBu\20170906\');
allPathname{5} = strcat(dataDir,'TuBu\20170913\');
allPathname{6} = strcat(dataDir,'TuBu\20180327\');
%% Specify the directories - Empty Gal4
numDirs = 5;
saveName = 'EmptyGal4';
allPathname{1} = strcat(dataDir,'ctrl\20170818\');
allPathname{2} = strcat(dataDir,'ctrl\20170819\');
allPathname{3} = strcat(dataDir,'ctrl\20170905\');
allPathname{4} = strcat(dataDir,'ctrl\20170907\');
allPathname{5} = strcat(dataDir,'ctrl\20180327\');
%% Specify the directories - E-PGs
numDirs = 3;
saveName = 'SS96';
allPathname{1} = strcat(dataDir,'EPG\20180223\');
allPathname{2} = strcat(dataDir,'EPG\20180225\');
allPathname{3} = strcat(dataDir,'EPG\20180226\');

%% Load the positional data
flyID = 0;
for dirNow = 1:numDirs
    dirParts = strsplit(allPathname{dirNow},'\');
    fileNames = dir(allPathname{dirNow});
    for moveID = 3:length(fileNames)
        moveName =fileNames(moveID).name;
        moveNameParts = strsplit(moveName,'_');
        pathNameParts = strsplit(allPathname{dirNow},'\');
        % Find the associated position data 
        if (contains(moveName,'.TXT') & ~contains(moveName,'SYNC'))
            flyName = strcat(pathNameParts{end-1},moveNameParts{1});
            if flyID == 0
                flyID = 1;
                allPosData{flyID}.name = flyName
                allPosData{flyID}.Stripe_30C = {};
                allPosData{flyID}.Stripe = {};
            elseif ~strcmp(allPosData{flyID}.name,flyName)
                flyID = flyID + 1;
                allPosData{flyID}.name = flyName
                allPosData{flyID}.Stripe_30C = {};
                allPosData{flyID}.Stripe = {};
            end
            moveName
            positionDat = VRDatLoad(moveName,allPathname{dirNow},0);
            
            num_planes = 10;
            tSpan = positionDat.t - positionDat.t(1) + positionDat.tVR(1)/10000;
            minFG = ceil(min(find(positionDat.tFrameGrab >= positionDat.tVR(1)))/num_planes);
            maxFG = round(length(positionDat.tFrameGrab)/num_planes);
            positionDat.num_planes = num_planes;
            positionDat.minFG = minFG;
            positionDat.maxFG = maxFG;

            % convert the raw ball position data to position (in case
            % the visual display is not closed loop with a gain of 1)
            [posRot, posFor, posLat] = PositionConverter(positionDat,dirParts{end-1});

            % Match the position data to the framegrab times
            OffsetRotMatch = MatchData(positionDat.t,positionDat);
            OffsetRotMatch(:,2) = MatchData(pi/180*positionDat.OffsetRot,positionDat);
            OffsetForMatch = MatchData(positionDat.OffsetFor,positionDat);
            OffsetLatMatch = MatchData(positionDat.OffsetLat,positionDat);
            PosRotMatch = MatchData(posRot,positionDat);
            PosForMatch = MatchData(posFor,positionDat);
            PosLatMatch = MatchData(posLat,positionDat);
            Closed = MatchData(positionDat.closed,positionDat);
            Direction = MatchData(positionDat.direction,positionDat);
            Trans = MatchData(positionDat.trans,positionDat);
            if isfield(positionDat,'olgain')
                OLGain = MatchData(positionDat.olgain,positionDat);
            end
            if isfield(positionDat,'clgain')
                CLGain = MatchData(positionDat.clgain,positionDat);
            end

            if isfield(positionDat,'tStim')
                tStim = zeros(positionDat.maxFG-positionDat.minFG+1,1);
                for stim = 1:length(positionDat.tStim)
                    tMatch = find(positionDat.tFrameGrab > positionDat.tStim(stim));
                    volNow = ceil(tMatch(1)./positionDat.num_planes) - ...
                        positionDat.minFG+1;
                    tStim(volNow) = 1;
                end
            end


            % Unwrap the rotation and calculate the velocities
            OffsetRotUnwrap = UnWrap(pi/180*PosRotMatch,2,0);
            vRot = diff(OffsetRotUnwrap)./mean(diff(OffsetRotMatch(:,1)));
            vF = sqrt(diff(PosForMatch).^2+diff(PosLatMatch).^2)./mean(diff(OffsetRotMatch(:,1)));

            % populate the positionDatMatch cell
            positionDatMatch = {};
            positionDatMatch.OffsetRotMatch = OffsetRotMatch;
            positionDatMatch.OffsetForMatch = OffsetForMatch;
            positionDatMatch.OffsetLatMatch = OffsetLatMatch;
            positionDatMatch.vRot = vRot;
            positionDatMatch.vF = vF;
            positionDatMatch.minFG = minFG;
            positionDatMatch.maxFG = maxFG;
            positionDatMatch.PosRotMatch = PosRotMatch;
            positionDatMatch.PosForMatch = PosForMatch;
            positionDatMatch.PosLatMatch = PosLatMatch;
            positionDatMatch.Closed = Closed;
            positionDatMatch.Direction = Direction;
            positionDatMatch.Trans = Trans;
            if isfield(positionDat,'olgain')
                positionDatMatch.OLGain = OLGain;
            end
            if isfield(positionDat,'clgain')
                positionDatMatch.CLGain = CLGain;
            end
            if isfield(positionDat,'tStim')
                positionDatMatch.tStim = tStim;
            end
                    
            trialParts = strsplit(moveNameParts{end},'.');
            trialID = str2num(trialParts{1});
            if contains(moveName,'30C')
                allPosData{flyID}.Stripe_30C{trialID}.positionDatMatch = positionDatMatch;
            else
                allPosData{flyID}.Stripe{trialID}.positionDatMatch = positionDatMatch;
            end
        end
    end  
end

%% Plot the stripe position distribution only when the animal is moving - Figure S12 A-D

% Stretch the arena to cover 360 deg
visMult = 360/240;

% Set a threshold on movement
vRThresh = 0.1*pi;
vFThresh = 0.1;

% Specify plotting parameters
pHistEdges = linspace(-pi,pi,32);
innerD = 0.8;
outerD = 1;
       
% Plot the data for each fly
for flyID = 1:length(allPosData)

    % Plot the data
    for hc = 1:2
        if hc==1
            trial = 'Stripe';
        else
            trial = 'Stripe_30C';
        end
        for rt = 1:length(allPosData{flyID}.(trial))

            colorTxt = [0 0 0];
            
            % Find where the stripe appears
            stripeVis = find(diff(allPosData{flyID}.(trial){rt}.positionDatMatch.Direction));
            stripeVis = stripeVis(1)+1;
            if stripeVis > 1000
                stripeVis = 1;
            end

            % Find where the stripe first jumps
            stripeJump = find(diff(allPosData{flyID}.(trial){rt}.positionDatMatch.Trans));
            stripeJump = stripeJump+1;
            if length(stripeJump) == 1
                stripeJump(2) = length(allPosData{flyID}.(trial){rt}.positionDatMatch.Trans)-1;
            end
            
            % Identify the dark and first and second stripe periods
            darkPer = [1:stripeVis-1];
            if ~isempty(stripeJump)
                stripePer1 = [stripeVis:stripeJump(1)-1];
                stripePer2 = [stripeJump(1):stripeJump(2)-1];
            else
                stripePer1 = [stripeVis:length(allPosData{flyID}.(trial){rt}.positionDatMatch.Direction)];
                stripePer2 = [];
            end

            % Only consider periods where the animal is moving
            vRot = allPosData{flyID}.(trial){rt}.positionDatMatch.vRot;
            vF = allPosData{flyID}.(trial){rt}.positionDatMatch.vF;
            rotAn = find(abs(vRot) > vRThresh);
            forAn = find(vF > vFThresh);
            movingAn = union(rotAn,forAn);

            % Plot histograms of the stripe position
            tAll = allPosData{flyID}.(trial){rt}.positionDatMatch.OffsetRotMatch(:,1);
            RotAll = allPosData{flyID}.(trial){rt}.positionDatMatch.OffsetRotMatch(:,2);

            RotNow = visMult*RotAll(intersect(stripePer1,movingAn));
            [N1,angEdge1] = histcounts(RotNow,pHistEdges);
            RotNow = visMult*RotAll(intersect(stripePer2,movingAn));
            [N2,angEdges2] = histcounts(RotNow,pHistEdges);

            arrowAng = circ_mean((angEdge1(1:end-1)+mean(diff(angEdge1))/2)',N1');
            arrowL = innerD*circ_r((angEdge1(1:end-1)+mean(diff(angEdge1))/2)',N1');
            allPosData{flyID}.(trial){rt}.prejump.ang = arrowAng;
            allPosData{flyID}.(trial){rt}.prejump.L = arrowL;
            
            arrowAng = circ_mean((angEdge1(1:end-1)+mean(diff(angEdge1))/2)',N2');
            arrowL = innerD*circ_r((angEdge1(1:end-1)+mean(diff(angEdge1))/2)',N2');
            allPosData{flyID}.(trial){rt}.postjump.ang = arrowAng;
            allPosData{flyID}.(trial){rt}.postjump.L = arrowL;


        end
    end
end

% Create plot
StripePosStats = figure('units','normalized','outerposition',[0 0 0.25 0.5]);

subplot(4,3,[1 2 4 5]);
hold on;
rectangle('Position',[-outerD -outerD 2*outerD 2*outerD],'Curvature',1,'EdgeColor','k');
rectangle('Position',[-0.75*outerD -0.75*outerD 1.5*outerD 1.5*outerD],'Curvature',1,'EdgeColor',[0.25 0.25 0.25]);
rectangle('Position',[-0.5*outerD -0.5*outerD outerD outerD],'Curvature',1,'EdgeColor',[0.25 0.25 0.25]);
rectangle('Position',[-0.25*outerD -0.25*outerD 0.5*outerD 0.5*outerD],'Curvature',1,'EdgeColor',[0.25 0.25 0.25]);

subplot(4,3,[7 8 10 11]);
hold on;
rectangle('Position',[-outerD -outerD 2*outerD 2*outerD],'Curvature',1,'EdgeColor','k');
rectangle('Position',[-0.75*outerD -0.75*outerD 1.5*outerD 1.5*outerD],'Curvature',1,'EdgeColor',[0.25 0.25 0.25]);
rectangle('Position',[-0.5*outerD -0.5*outerD outerD outerD],'Curvature',1,'EdgeColor',[0.25 0.25 0.25]);
rectangle('Position',[-0.25*outerD -0.25*outerD 0.5*outerD 0.5*outerD],'Curvature',1,'EdgeColor',[0.25 0.25 0.25]);

% Specify colors for each fly
colCode = jet(2*length(allPosData));

% Set up a vector to hold the resultant vector lengths and angles
vecLengths = cell(2,1);
angles = cell(2,1);
angDiffs = cell(2,1);

% Plot the data for each fly
for flyID = 1:length(allPosData)
    % Plot the data
    for hc = 1:2
        if hc==1
            trial = 'Stripe';
            trial1 = 1;
        else
            trial = 'Stripe_30C';
            trial1 = 2;
        end
        for rt = trial1:length(allPosData{flyID}.(trial))
            
            subplot(4,3,6*(hc-1) + [1 2 4 5]);
            hold on;
            
            arrowAng1 =  allPosData{flyID}.(trial){rt}.prejump.ang;
            angles{hc} = [angles{hc},arrowAng1];
            arrowL = allPosData{flyID}.(trial){rt}.prejump.L;
            vecLengths{hc} = [vecLengths{hc},arrowL];
            scatter(arrowL*sin(arrowAng1),...
                arrowL*cos(arrowAng1),...
                50,...
                colCode(2*(flyID-1)+1+(rt-trial1),:),...
                'filled');
            scatter(1.05*sin(arrowAng1),...
                1.05*cos(arrowAng1),...
                25,...
                'k',...
                'filled');
            
            arrowAng2 =  allPosData{flyID}.(trial){rt}.postjump.ang;
            angles{hc} = [angles{hc},arrowAng2];
            arrowL = allPosData{flyID}.(trial){rt}.postjump.L;
            vecLengths{hc} = [vecLengths{hc},arrowL];
            scatter(arrowL*sin(arrowAng2),...
                arrowL*cos(arrowAng2),...
                50,...
                colCode(2*(flyID-1)+1+(rt-trial1),:),...
                'filled');
            scatter(1.05*sin(arrowAng2),...
                1.05*cos(arrowAng2),...
                25,...
                'k',...
                'filled');
            
            angDiffs{hc} = [angDiffs{hc},arrowAng1 - arrowAng2];

            if rt == 1 & hc == 1
                title('All flies','FontSize',10,'FontWeight','normal','color','k');
            end
            xlim([-1.05*outerD 1.05*outerD]);
            ylim([-1.05*outerD 1.05*outerD]);
            axis square;
            axis off;
            line([0 0],[0 outerD],'Color','k','LineStyle','--');

        end
    end
end
   

% Create plot
StripePosStats = figure('units','normalized','outerposition',[0 0 0.25 0.5]);

subplot(2,2,1);
hold on;
[counts,bins] = hist(mod(angDiffs{1}-pi,2*pi)-pi,[-pi:pi/4:pi]); %# get counts and bin locations
bar(bins,counts,'b')
ylim([0 10]);
xticks([-pi:pi/2:pi]);
ylabel('counts');
xlabel('ang. pref change');

subplot(2,2,2);
hold on;
[counts,bins] = hist(angles{1},[-pi:pi/4:pi]); %# get counts and bin locations
bar(bins,counts,'b')
ylim([0 15]);
xticks([-pi:pi/2:pi]);
xlabel('mean angle');
ylabel('counts');

subplot(2,2,3);
hold on;
[counts,bins] = hist(mod(angDiffs{2}-pi,2*pi)-pi,[-pi:pi/4:pi]); %# get counts and bin locations
bar(bins,counts,'r')
ylim([0 10]);
xticks([-pi:pi/2:pi]);
ylabel('counts');
xlabel('ang. pref change');

subplot(2,2,4);
hold on;
[counts,bins] = hist(angles{2},[-pi:pi/4:pi]); %# get counts and bin locations
bar(bins,counts,'r')
ylim([0 15]);
xticks([-pi:pi/2:pi]);
xlabel('mean angle');
ylabel('counts');