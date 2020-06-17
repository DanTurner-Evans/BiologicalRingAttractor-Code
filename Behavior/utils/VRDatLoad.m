%% Load the imaging and DAQ data from a 2D VR experiment
function positionDat = VRDatLoad(posFilename,posPathname,recStim)

% Get the DAQ info
SYNCFilename = strcat(posFilename(1:end-4),'_SYNC',posFilename(end-3:end));
SYNCPathname = posPathname;

% Pull out the fly positional information
fileID = fopen(strcat(posPathname,posFilename));
tstamp = fgetl(fileID);
exTypeStr = strsplit(tstamp,'_');
exType = exTypeStr{end}(1:end-4);
formatSpec = '%s %f %s %f %s %f %s %f %s %f %s %f %s %d %s %d %s %d %s %d %s %d %s %f %s %f';
N=400000;
C = textscan(fileID,formatSpec,N,'CommentStyle','Current','Delimiter','\t');
t = C{1,2}; % Time
OffsetRot = C{1,4}; % Stripe rotational offset
OffsetRot = mod(OffsetRot+180, 360)-180;
OffsetFor = C{1,6}; % Stripe forward offset
OffsetLat = C{1,8}; % Stripe lateral offset
dx0 = C{1,10}; % X position of the ball from camera 1 
dx1 = C{1,12}; % X position of the ball from camera 2
dy0 = C{1,14};
dy1 = C{1,16};
closed = C{1,18};
direction = C{1,20};
trans = C{1,22};
olgain = C{1,24};
clgain = C{1,26};
fclose(fileID);

numDatPts = length(OffsetFor);

positionDat = {};
positionDat.t = t(1:numDatPts);
positionDat.OffsetRot = OffsetRot(1:numDatPts);
positionDat.OffsetFor = OffsetFor;
positionDat.OffsetLat = OffsetLat;
positionDat.dx0 = dx0;
positionDat.dx1 = dx1;
positionDat.dy0 = dy0;
positionDat.dy1 = dy1;
positionDat.closed = closed;
positionDat.direction = direction;
positionDat.trans = trans;
positionDat.olgain = olgain;
positionDat.clgain = clgain;
positionDat.exType = exType;

% Load photodiode and frame time stamps
timeStamps = importdata(strcat(SYNCPathname,SYNCFilename));

% Pull out the time stamps for the frame grab signal
clear tFrameGrab;
sampleData = 1;
upperLim = max(timeStamps(:,1));
startFrameGrab = find(timeStamps(:,1) > 0.8*upperLim);
incDat = startFrameGrab(1)-2;
inct = 1;
while (sampleData)
    if (timeStamps(incDat+1,1) > 0.8*upperLim && timeStamps(incDat,1) < 0.8*upperLim)
        tFrameGrab(inct) = incDat+1;
        inct = inct +1;
        incDat = incDat + 1;
    end
    incDat=incDat+1;
    if incDat > length(timeStamps)-1
        break
    end
end

positionDat.tFrameGrab = tFrameGrab;

% Pull out the time stamps for the VR refresh
clear tVR;
sampleData = 1;
upperLim = max(timeStamps(:,2));
offset = round(0.6/(360)*10000);
VRthresh = 0.4;
startVR = find(timeStamps(:,2) > VRthresh*upperLim);
if startVR(1) == 1
    startVR = startVR + 1;
end
incDat = startVR(1)-2;
inct = 1;
while (sampleData)
    if (timeStamps(incDat+1,2) < VRthresh*upperLim && (timeStamps(incDat-1,2) < timeStamps(incDat+1,2) || timeStamps(incDat,2) < timeStamps(incDat+1,2)))
        tVR(inct) = incDat+1;
        inct = inct +1;
        incDat = incDat + offset;
    end
    incDat=incDat+1;
    if incDat > length(timeStamps)-1
        break
    end
end

positionDat.tVR = tVR;

if recStim
    % Pull out the time stamps for the stimulation
    clear tStim;
    sampleData = 1;
    upperLim = max(timeStamps(:,3));
    startFrameGrab = find(timeStamps(:,3) > 0.9*upperLim);
    incDat = startFrameGrab(1)-2;
    inct = 1;
    while (sampleData)
        if (timeStamps(incDat+1,3) > 0.9*upperLim && timeStamps(incDat,3) < 0.9*timeStamps(incDat+1,3))
            tStim(inct) = incDat+1;
            inct = inct +1;
            incDat = incDat + offset;
        end
        incDat=incDat+1;
        if incDat > length(timeStamps)-1
            break
        end
    end

    positionDat.tStim = tStim;
end

end
