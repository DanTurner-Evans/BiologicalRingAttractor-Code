function [dataRLPB, dataRRPB, dataGLPB, dataGRPB] = BumpAlignOnePeak(allFlyData, Min, Max, Span, rangeR, rangeG, trial, glomShift)
    
%function to sort data into bins

%input:
%allFlyData = the data to be sorted
%Min = the minimum velocity bin (deg)
%Max = the maximum velocity bin (deg)
%Span = the span of each bin (deg)
%RG = Switch to tell whether 60D05 is green (0) or red (1)
%LR = Switch to tell whether we're looking at the left PB (0) or right
%PB (1)
%trial = the trial to be considered

%output:
%dataRLPB = the sorted data for the left PB red channel
%dataRRPB = the sorted data for the right PB red channel
%dataGLPB = the sorted data for the left PB green channel
%dataGRPB = the sorted data for the right PB green channel

% Specify the range of rotational velocities to consider
vRMin = Min*pi/180;
vRMax = Max*pi/180;
vRSpan = Span*pi/180;
vRBinNum = round((vRMax-vRMin)/vRSpan);

% Initilaize the sorted output data arrays
dataRLPB = cell(size(allFlyData,2),1);
dataRRPB = cell(size(allFlyData,2),1);
dataGLPB = cell(size(allFlyData,2),1);
dataGRPB = cell(size(allFlyData,2),1);

% Specify which range of ROIs to consider
if rangeR(1) == 2 || rangeR(1) == 10
    RG = 0;
else
    RG = 1;
end

% Sort the data
for flyID = 1:length(allFlyData)
    dataRLPB{flyID}.CW = cell(vRBinNum,1);
    dataRLPB{flyID}.CCW = cell(vRBinNum,1);
    dataRLPB{flyID}.Stop = [];  
    dataRRPB{flyID}.CW = cell(vRBinNum,1);
    dataRRPB{flyID}.CCW = cell(vRBinNum,1);
    dataRRPB{flyID}.Stop = [];  
    
    dataGLPB{flyID}.CW = cell(vRBinNum,1);
    dataGLPB{flyID}.CCW = cell(vRBinNum,1);
    dataGLPB{flyID}.Stop = [];  
    dataGRPB{flyID}.CW = cell(vRBinNum,1);
    dataGRPB{flyID}.CCW = cell(vRBinNum,1);
    dataGRPB{flyID}.Stop = [];  

    for trialID = 1:length(allFlyData{flyID}.(trial))
        vRot = allFlyData{flyID}.(trial){trialID}.positionDatMatch.vRot;
        RSig = allFlyData{flyID}.(trial){trialID}.RROIaveMax;
        GSig = allFlyData{flyID}.(trial){trialID}.GROIaveMax;
        for tStep = 1:length(vRot)
            if RG
                pkMax = max(GSig(rangeG,tStep+1));
                pkShift = -find(GSig(rangeG,tStep+1) == pkMax);
            else
                pkMax = max(RSig(rangeR,tStep+1));
                pkShift = -find(RSig(rangeR,tStep+1) == pkMax);
            end
            RSigNowLPB = zeros(9,1);
            RSigNowRPB = zeros(9,1);
            GSigNowLPB = zeros(9,1);
            GSigNowRPB = zeros(9,1);
            if rangeR(1) < 9
                RSigNowLPB(rangeR) = circshift(RSig(rangeR,tStep+1),pkShift+glomShift);
                RSigNowRPB(rangeG) = circshift(RSig(sort(19-rangeR),tStep+1),pkShift+glomShift);
                GSigNowLPB(rangeG) = circshift(GSig(rangeG,tStep+1),pkShift+glomShift);
                GSigNowRPB(rangeR) = circshift(GSig(sort(19-rangeG),tStep+1),pkShift+glomShift);
            else
                RSigNowLPB(sort(19-rangeR)) = circshift(RSig(sort(19-rangeR),tStep+1),pkShift+glomShift);
                RSigNowRPB(sort(19-rangeG)) = circshift(RSig(rangeR,tStep+1),pkShift+glomShift);
                GSigNowLPB(sort(19-rangeG)) = circshift(GSig(sort(19-rangeG),tStep+1),pkShift+glomShift);
                GSigNowRPB(sort(19-rangeR)) = circshift(GSig(rangeG,tStep+1),pkShift+glomShift);
            end
            
            if vRot(tStep) == 0
                dataRLPB{flyID}.Stop = horzcat(dataRLPB{flyID}.Stop, RSigNowLPB);
                dataRRPB{flyID}.Stop = horzcat(dataRRPB{flyID}.Stop, RSigNowRPB);
                dataGLPB{flyID}.Stop = horzcat(dataGLPB{flyID}.Stop, GSigNowLPB);
                dataGRPB{flyID}.Stop = horzcat(dataGRPB{flyID}.Stop, GSigNowRPB);

            elseif vRot(tStep) > 0
                dataRLPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataRLPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    RSigNowLPB);
                dataRRPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataRRPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    RSigNowRPB);
                dataGLPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataGLPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    GSigNowLPB);
                dataGRPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)} = ...
                    horzcat(dataGRPB{flyID}.CCW{ceil(vRot(tStep)/vRSpan)},...
                    GSigNowRPB);
            else
                dataRLPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataRLPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    RSigNowLPB);
                dataRRPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataRRPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    RSigNowRPB);
                dataGLPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataGLPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    GSigNowLPB);
                dataGRPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)} = ...
                    horzcat(dataGRPB{flyID}.CW{ceil(-vRot(tStep)/vRSpan)},...
                    GSigNowRPB);
            end
        end
    end
end
    