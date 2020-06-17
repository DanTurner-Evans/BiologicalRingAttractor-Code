%% Load the imaging and DAQ data from a two color, 2D VR experiment
function [RstackMaxIntRot, GstackMaxIntRot, RstackMeanRot, GstackMeanRot] = ImDatLoadBigtiff2Color(imageFilename,imagePathname,rotAng)

% Get the tiff header info
fullpath = strcat(imagePathname,imageFilename);
reader=ScanImageTiffReader(fullpath);
desc = reader.metadata;
planeLoc = strfind(desc, 'numFramesPerVolume');
discardLoc = strfind(desc, 'numDiscardFlybackFrames');
num_planes = str2num(desc(planeLoc+21:planeLoc+22));
num_discards = str2num(desc(discardLoc+26));

% Load the tifStack
vol=ScanImageTiffReader(fullpath).data();
width = size(vol,1);
height = size(vol,2);
num_images = size(vol,3);
numVolumes = floor(num_images/num_planes/2);

% Find the max intensity projection
RstackMaxInt = double(zeros(height,width,numVolumes));
GstackMaxInt = double(zeros(height,width,numVolumes));
RminiStack = double(zeros(height,width,(num_planes-num_discards)));
GminiStack = double(zeros(height,width,(num_planes-num_discards)));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
for incIm = 1:num_images
    if mod(incIm,100)==0
        waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
    end
    if mod(ceil(incIm/2),num_planes) ~= 0
        if mod(incIm,2)
            RminiStack(:,:,mod(ceil(incIm/2),num_planes)) = vol(:,:,incIm)';
        else
            GminiStack(:,:,mod(ceil(incIm/2),num_planes)) = vol(:,:,incIm)';
        end
        if (mod(incIm/2,num_planes) == num_planes-num_discards)
            RstackMaxInt(:,:,ceil(incIm/(2*num_planes))) = max(RminiStack,[],3);
            GstackMaxInt(:,:,ceil(incIm/(2*num_planes))) = max(GminiStack,[],3);
        end
    end
end
delete(h);

% Rotate the stacks so that the EB is aligned
RMean = mean(RstackMaxInt,3);
RMean = imrotate(RMean,rotAng);
GMean = mean(GstackMaxInt,3);
GMean = imrotate(GMean,rotAng);
RstackMaxIntRot = zeros([size(RMean) size(RstackMaxInt,3)]);
GstackMaxIntRot = zeros([size(GMean) size(GstackMaxInt,3)]);
for frame = 1:size(RstackMaxIntRot,3)
    RstackMaxIntRot(:,:,frame) = imrotate(RstackMaxInt(:,:,frame),rotAng);
    GstackMaxIntRot(:,:,frame) = imrotate(GstackMaxInt(:,:,frame),rotAng);
end

% Find the mean volume projection
RstackMean = double(zeros(height,width,numVolumes));
GstackMean = double(zeros(height,width,numVolumes));
h = waitbar(0.0,'Loading TIFF stack (mean calc)...');
set(h,'Position',[50 50 360 72]);
set(h,'Name','Loading TIFF stack...');
for incIm = 1:num_images
    if mod(incIm,100)==0
        waitbar(incIm/num_images,h,['Loading frame# ' num2str(incIm) ' out of ' num2str(num_images)]);
    end
    if mod(ceil(incIm/2),num_planes) ~= 0
        if mod(incIm,2)
            RstackMean(:,:,ceil(incIm/(2*num_planes))) = double(vol(:,:,incIm)')+ RstackMean(:,:,ceil(incIm/(2*num_planes)));
        else
            GstackMean(:,:,ceil(incIm/(2*num_planes))) = double(vol(:,:,incIm)')+ GstackMean(:,:,ceil(incIm/(2*num_planes)));
        end
    end
end
delete(h);

% Rotate the stacks so that the EB is aligned
RMean = mean(RstackMean,3);
RMean = imrotate(RMean,rotAng);
GMean = mean(GstackMean,3);
GMean = imrotate(GMean,rotAng);
RstackMeanRot = zeros([size(RMean) size(RstackMean,3)]);
GstackMeanRot = zeros([size(GMean) size(GstackMean,3)]);
for frame = 1:size(RstackMeanRot,3)
    RstackMeanRot(:,:,frame) = imrotate(RstackMean(:,:,frame),rotAng);
    GstackMeanRot(:,:,frame) = imrotate(GstackMean(:,:,frame),rotAng);
end

end