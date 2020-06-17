%% Load the imaging and DAQ data from a 2D VR experiment
function [stackMaxIntRot, stackMeanRot] = ImDatLoadBigtiff(imageFilename,imagePathname,rotAng)

% Get the tiff header info
fullpath = strcat(imagePathname,imageFilename);
reader=ScanImageTiffReader(fullpath);
desc = reader.metadata;
planeLoc = strfind(desc, 'numFramesPerVolume');
discardLoc = strfind(desc, 'numDiscardFlybackFrames');
num_planes = str2num(desc(planeLoc+21:planeLoc+22));
if isempty(num_planes)
    num_planes = 1;
end
num_discards = str2num(desc(discardLoc+26));

% Load the tifStack
vol=ScanImageTiffReader(fullpath).data();
width = size(vol,1);
height = size(vol,2);
num_images = size(vol,3);
numVolumes = floor(num_images/num_planes);

% Find the max intensity projection
stackMaxInt = double(zeros(height,width,numVolumes));
miniStack = double(zeros(height,width,num_planes-num_discards));
h = waitbar(0.0,'Loading TIFF stack (max calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
for incIm = 1:num_images
    if mod(incIm,100)==0
        waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
    end
    if mod(incIm,num_planes) ~= 0 & mod(incIm,num_planes) < num_planes-num_discards+1
        miniStack(:,:,mod(incIm,num_planes)) = vol(:,:,incIm)';
        if mod(incIm,num_planes) == num_planes-num_discards
            stackMaxInt(:,:,ceil((incIm)/num_planes)) = max(miniStack,[],3);
        end
    end
end
delete(h);

% Rotate the stacks so that the EB is aligned
stkMean = mean(stackMaxInt,3);
stkMean = imrotate(stkMean,rotAng);
stackMaxIntRot = zeros([size(stkMean) size(stackMaxInt,3)]);
for frame = 1:size(stackMaxIntRot,3)
    stackMaxIntRot(:,:,frame) = imrotate(stackMaxInt(:,:,frame),rotAng);
end

% Find the mean volume projection
stackMean = double(zeros(height,width,numVolumes));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
for incIm = 1:num_images
    if mod(incIm,100)==0
        waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
    end
    if mod(incIm,num_planes) ~= 0
        stackMean(:,:,ceil(incIm/num_planes)) = double(vol(:,:,incIm)')+ stackMean(:,:,ceil(incIm/num_planes));
    end
end
delete(h);

% Rotate the stacks so that the EB is aligned
stkMean = mean(stackMean,3);
stkMean = imrotate(stkMean,rotAng);
stackMeanRot = zeros([size(stkMean) size(stackMean,3)]);
for frame = 1:size(stackMean,3)
    stackMeanRot(:,:,frame) = imrotate(stackMean(:,:,frame),rotAng);
end

end