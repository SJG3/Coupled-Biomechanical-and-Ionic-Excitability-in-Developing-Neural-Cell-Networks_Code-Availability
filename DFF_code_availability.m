
clear all, close all, clc; 
% this code is used to load a Tif or ND2 image, find all cell center using
% difference of gausian approach with otsu thresholding and then extract
% the flourescent intensity and compute delta F/F. These can be used for
% downstream analaysis like correlation analysis 


%% Load Tiff File                            
% Use uigetfile to select the .tif file
[file, path] = uigetfile('*.tif', 'Select a TIF file');

% Check if the user clicked on Cancel
if isequal(file, 0)
    disp('File selection canceled.');
    return;
end
tic;
% Combine the path and file name to get the full file path
fullFilePath = fullfile(path, file);

% Get information about the TIFF file
info = imfinfo(fullFilePath);
numFrames = numel(info);

% Preallocate the 3D matrix to store the image stack (X, Y, T)
IMG = zeros(info(1).Height, info(1).Width, numFrames, 'uint16');

% Read each frame and store it in the 3D matrix
for t = 1:numFrames
    IMG(:,:,t) = imread(fullFilePath, 'Index', t, 'Info', info);
end
toc;
IMG = im2double(IMG);


% Check if the image was successfully loaded
if isempty(IMG)
    error('Error loading the TIF image file.');
else
    disp('TIF image loaded successfully.');
end

IMGmean = mean(IMG,3);
%% Load ND2 File

addpath('/Users/kos-mos/Documents/bfmatlab');

[file, path] = uigetfile('*.nd2', 'Select a TIF file');

% Check if the user clicked on Cancel
if isequal(file, 0)
    disp('File selection canceled.');
    return;
end
tic;
% Combine the path and file name to get the full file path
fullFilePath = fullfile(path, file);

data = bfopen(fullFilePath);

% Get information about the TIFF file
% info = imfinfo(fullFilePath);
numFrames = size(data{1,1},1);

% Preallocate the 3D matrix to store the image stack (X, Y, T)
IMG = zeros(size(data{1,1}{1},1), size(data{1,1}{1},2), numFrames, 'uint16');

% Read each frame and store it in the 3D matrix
for t = 1:numFrames
    IMG(:,:,t) = data{1,1}{t};
end
toc;
IMG = im2double(IMG);


% Check if the image was successfully loaded
if isempty(IMG)
    error('Error loading the TIF image file.');
else
    disp('image loaded successfully.');
end

IMGmean = mean(IMG,3);
%% Parameters
val = randi(100);

FPS = 20; %capture rate 

expectedRadiusPx = 7; 
window_size = 15; %300;
baseline_percentile = 75;

%% Time Corrections                     
TimeTotal = 1/FPS * t; 
TimeVals = seconds(0:TimeTotal/(t-1):TimeTotal);

XTimes = seconds(TimeVals); 

% Initialize the variable 'TimeIn' to 0 initially
TimeIn = "seconds";

% Create a cell array of options
options = {'Seconds', 'Minutes', 'Hours'};

% Display a menu to the user and get their choice
choice = menu('Select a Time Display:', options);

% Check the user's choice and update the 'TimeIn' variable accordingly
if choice == 1
    TimeIn = "seconds";
    XTimes = seconds(TimeVals);
elseif choice == 2
    TimeIn = "minutes";
    XTimes = minutes(TimeVals);
elseif choice == 3
    TimeIn = "hours";
    XTimes = hours(TimeVals);
end

% Display the selected option
fprintf('You selected time to display in: %s\n', options{choice});
%% Define the standard deviations for the two Gaussian filters (large and small)
sigma_large = 6;
sigma_small = 5;
%% DoG automated cell finder using Otsu on image            
clear center_x center_y
% Assuming you have the IMGmean variable (2D matrix, X, Y) already loaded

% Apply Gaussian blurring using fspecial and imfilter functions
large_blur = imgaussfilt(IMGmean, sigma_large);
small_blur = imgaussfilt(IMGmean, sigma_small);

% Compute the Difference of Gaussians (DoG) image
DoG_image = mat2gray(small_blur - large_blur);

thresh_Factor = 1.0; % play around with this until you get value that gives multiple dots corresponding with cells

Objs_image = (DoG_image) > graythresh(DoG_image) * thresh_Factor;
imagesc(Objs_image); 

BW_props = regionprops(Objs_image,'centroid','circularity','Area');
        for i = 1:length(BW_props)
            if BW_props(i).Circularity >=0.25 && BW_props(i).Area > 10
                center_y(i,1) = BW_props(i).Centroid(2);
                center_x(i,1) = BW_props(i).Centroid(1);
            else
                center_x(i,1) = NaN;
                center_y(i,1) = NaN;
            end
        end
center_x = rmmissing(center_x);
center_y = rmmissing(center_y);
%% Display the original image and the identified centers
figure (2);
% subplot(2,2,1);
%     imagesc(Objs_image); 
%     hold on;
%     plot(center_x, center_y, 'm.', 'MarkerSize', 10, 'LineWidth', 2);  axis image;
%     title('Identified Centers of Objects using DoG');
%     hold off;
    
% subplot(2,2,2);
    imagesc(imadjust(IMGmean)); 
%     caxis([0, 0.002])
    hold on;
    plot(center_x, center_y, 'm.', 'MarkerSize', 10, 'LineWidth', 2);  axis image;
    title( size(center_x,1) + " Identified Centers of Objects using DoG");
    hold off;
    
% subplot(2,2,3);
%     imagesc(bwlabel(Objs_image)); 
% %     caxis([])
%     hold on;
%     plot(center_x, center_y, 'm.', 'MarkerSize', 10, 'LineWidth', 2);  axis image;
%     title('Identified Centers of Objects using DoG');
%     hold off;
%% Extract Fluorecent Intensity (ZB) 

% [rawFluos,roiBW] = extractFluo(IMG,center_x,center_y,expectedRadiusPx);
[rawFluos,roiBW] = extractFluo(IMG,center_x,center_y,expectedRadiusPx);
ROIsMap = sum(cat(3,roiBW{:}),3); %in case you want to see all ROIs together
%% Calculate dF/F (ZB)                  
% Normalize the rawFluo matrix to the range [0, 1]
% normalized_rawFluo = (rawFluos - min(rawFluos(:))) / (max(rawFluos(:)) - min(rawFluos(:)));
tic;
normalized_rawFluo = normalize(rawFluos,1,"norm");
%Call the calculateDFF function
DFF = calculateDFF(normalized_rawFluo', window_size, baseline_percentile);
toc;
%% PLOT DFFs

% Generate rwb colormap  
    % Define the red, white, and blue colors
    redColor = [1, 0, 0];
    whiteColor = [0.8, 0.8, 0.8];
    blueColor = [0, 0, 1];

    % Create the colormap
    numRows = 100;
    positions = [0, 0.5, 1];
    colors = [redColor; whiteColor; blueColor];
    rwbColormap = interp1(positions, colors, linspace(0, 1, numRows));

% Plot Kymograph    
    cval = max(DFF(:))*1.01; %Or make 1 to have max and min be 1
    figure(val); imagesc( XTimes, 1:size(DFF,1), DFF); colormap (rwbColormap); caxis([-cval, cval]);
    colorbar;
    xlabel("Time ("+TimeIn+")");
    ylabel('Index #');
  
% Plot Traces Overlapping
    figure(val+1); plot(XTimes, DFF');
    xlabel("Time ("+TimeIn+")");
    ylabel('Delta F/F Normalized Intensity');
    ax = gca; % Get the current Axes object
    ax.ColorOrder = turbo(size(DFF,1));

% Plot Traces on own row with an offset in the y-dimension
    figure(val+2);
    offset = 1;
    hold on;
%     colormap('turbo');
    for i = 1:size(DFF, 1)
        y_offset = (i - 1) * offset; % Calculate the offset for the current row
        plot(XTimes , normalize(DFF(i, :),"range") + y_offset); % Plot with offset
%         plot(XTimes , DFF(i, :) + y_offset); % Plot with offset

    end
    xlabel("Time ("+TimeIn+")");
    ax = gca; % Get the current Axes object
    ax.ColorOrder = turbo(size(DFF,1));
%     ax.Color = [0.6 0.6 0.6];
    ylim([0,max(DFF(i, :))+y_offset]);

    
% Plot Map of Cells and their index
    figure(val+3);
%     imagesc(bwlabel(Objs_image)); colormap ([0.1 0.1 0.1; turbo]);
    imagesc(imadjust(IMGmean)); colormap ([0.1 0.1 0.1; bone]); axis image;
%     caxis([])
    Plt_colors = turbo(size(center_x,1));
    
    hold on; 
    colororder(Plt_colors);
    scatter(center_x, center_y, 50, Plt_colors, 'LineWidth', 2); axis image;
    hold off;
    title('Identified Centers of Objects using DoG');

%% SAVE Variables

% Create the new folder
mkdir(path,"AnalysisOf_"+file);


% Save each variable as a separate .mat file
save(fullfile(path,"AnalysisOf_"+file , "DFF_"+file+"_"+date+".mat"), 'DFF');
save(fullfile(path,"AnalysisOf_"+file , "centers_XY_"+file+"_"+date+".mat"), 'center_x','center_y');

% Display a message to confirm the saving
disp('Variables saved successfully.');
 

%% FUNCTION Create ROI Map

function ROImap = createROImap(center_y, center_x, expectedRadiusPx, IMG)
    % Initialize the ROImap as a binary image with zeros
    ROImap = zeros(size(IMG,1), size(IMG,2));

    % Number of objects
    numObjects = numel(center_y);

    % Create circular ROIs at each object center
    for objIdx = 1:numObjects
        % Get the coordinates of the current object center
        y_center = center_y(objIdx);
        x_center = center_x(objIdx);

        % Create a circular mask around the current object center
        [x, y] = meshgrid(1:size(IMG,2), 1:size(IMG,1));
        mask = (x - x_center).^2 + (y - y_center).^2 <= expectedRadiusPx^2;

        % Set the circular ROI to 1 in the ROImap
        ROImap(mask) = 1;
    end
end
%% FUNCTION: Extract raw fluorescent intensity
function rawFluo = extractFluos(IMG, center_y, center_x, expectedRadiusPx)
    % Get the size of the image sequence (X, Y, T)
    [X, Y, T] = size(IMG);

    % Number of objects
    numObjects = numel(center_y);

    % Initialize the rawFluo matrix to store the fluorescent intensity values
    rawFluo = zeros(numObjects, T);

    % Create circular masks at each object center and extract fluorescent intensity
    for objIdx = 1:numObjects
        % Get the coordinates of the current object center
        y_center = center_y(objIdx);
        x_center = center_x(objIdx);

        % Create a circular mask around the current object center
        [x, y] = meshgrid(1:Y, 1:X);
        mask = (x - x_center).^2 + (y - y_center).^2 <= expectedRadiusPx^2;

        % Extract fluorescent intensity within the circular mask for each frame
        for t = 1:T
            frame = IMG(:, :, t);
            rawFluo(objIdx, t) = mean(frame(mask));
        end
    end
end
%% FUNCTION: Fluo extraction function ZacBowen-- I left in a lot of commented out code in case
% we go back and want to define cell-by-cell custom ROIs based on fluo.
% Right now all ROIs are the same size and shape.
function [rawFluo,roiBW] = extractFluo(IMG,xc,yc,expectedNeuronRadiusPix)


% narginchk(4,7)
% if ~exist('winSizeSeconds','var'); winSizeSeconds = 10; end % default to 10 second window
% if ~exist('percentBaselineSub','var'); percentBaselineSub = 50; end % default to 50% baseline subtraction

numNeurons = length(xc);

%% PREALLOCATE
traceExtract_start = tic;

% pcimg = cell (numNeurons , 13 , 360 );
% imgCrop = cell ( numNeurons , 27, 27 );
% imgCropNorm = cell (numNeurons, 28 , 28 );
% roiBoundaries = cell ( numNeurons , 360 , 3 );
% smRoiBoundaries = cell ( numNeurons , 360 , 3 );
% ROIOut = zeros ( 360 , 2 );

% 2-1-2023 SJG commented out and it seems to allow code to run tell Anna
ignor = 1;
if ignor == 0; 
    ROIxvOut = cell ( numNeurons , 360 , 1 );
    ROIyvOut = cell ( numNeurons , 360 , 1 );
    roiBW = cell ( numNeurons , size(IMG,1) , size(IMG,2));
end 
%% PREALLOCATE FLUO AND NPFLUO MATRICES
rawFluo = zeros( size(IMG,3) , numNeurons );

%% FIND THE BOUNDARIES OF CLICKED NEURONS
roiBounds = [deg2rad(1:360)' repmat(expectedNeuronRadiusPix,360,1)];
for pp = 1:numNeurons
    
%     % Find fluorescent ring using local peak finding on each neuron
%     xpt=xc(pp);
%     ypt=yc(pp);
%     imgCrop{pp} = imcrop(meanIMG,[xpt-15 ypt-15 31 31]);
%     imgCropNorm{pp} = (imgCrop{pp} - min(imgCrop{pp}(:)))  ./ (max(imgCrop{pp}(:)) - min(imgCrop{pp}(:)));
%     pcimg{pp} = imgpolarcoord (imgCropNorm{pp} );  % this comes from Matlab Central Download
%     
%     RingPks = zeros(size(pcimg{pp},2),1); %reset vals to zero
%     
%     tmpNeuron = pcimg{pp}; %iterate this for each selected neuron
%     
%     for cc = 1:size(tmpNeuron,2) % for every direction - find the inner part of the ring
%         
%         pkTmp = find(diff(tmpNeuron(:,cc)) == min(diff(tmpNeuron(:,cc)))); % DW07122015_changed to make this more robust - seems to be working right now - continue testing
%         
%         if ~isempty(pkTmp)
%             if length(pkTmp) > 1 %more than one pixel identified - grab the first one
%                 if pkTmp(1) < expectedNeuronRadiusPix && pkTmp(1) > 2
%                     RingPks(cc) = pkTmp(1);
%                 else
%                     RingPks(cc) = expectedNeuronRadiusPix;
%                 end
%             else
%                 if pkTmp < expectedNeuronRadiusPix && pkTmp(1) > 2
%                     RingPks(cc) = pkTmp;
%                 else
%                     RingPks(cc) = expectedNeuronRadiusPix;
%                 end
%             end
%         elseif cc == 1 % if it's the first direction and no peaks are found
%             RingPks(cc) = expectedNeuronRadiusPix;  %made this dependent on mag factor DW_02022015
%         else
%             RingPks(cc) = RingPks(cc-1);
%         end
%         ROIOut(cc,:) = [ deg2rad(cc)  RingPks(cc) ];
% 
%     end
%     roiBoundaries{pp} = [ ROIOut(:,1) ROIOut(:,2)]; % [PolarCoords (0-2Pi)     OuterRing]
%     smRoiBoundaries{pp} = [ ROIOut(:,1) smooth(ROIOut(:,2),10)]; % [PolarCoords (0-2pi)     OuterRing]
    
    % CREATE MASKS FOR ALL CLICKED ROIS, THEN SHOW THEM -- DW 11232015
    % renamed variable for consistency
%     ROIxvOut{pp} =  xpt + smRoiBoundaries{pp}(:,2) .* (cos(smRoiBoundaries{pp}(:,1))) ;
%     ROIyvOut{pp} =  ypt + smRoiBoundaries{pp}(:,2) .* (sin(smRoiBoundaries{pp}(:,1))) ;
    ROIxvOut{pp} =  xc(pp) + roiBounds(:,2) .* (cos(roiBounds(:,1))) ;
    ROIyvOut{pp} =  yc(pp) + roiBounds(:,2) .* (sin(roiBounds(:,1))) ;
    roiBW{pp} = poly2mask( ROIxvOut{pp} , ROIyvOut{pp} , size(IMG,1) , size(IMG,2));
end

%DW 11232015 - adjusted for inclusion of neuropil correction
% correct for overlapping ROIs (exclude from both)
disp('Adjusting ROI masks for overlap....');
tStartROICorr = tic;
AllMasksTMP =  sum ( cat ( 3 , roiBW{:} ) , 3 ); % first term of cat (i.e., '3') points to element-wise alignement/stacking of arrays
[oLapRoiY, oLapRoiX] = find( AllMasksTMP > 1 );
for ii = 1:numNeurons
    for yy = 1:length(oLapRoiX)
        roiBW{ii}(oLapRoiY(yy),oLapRoiX(yy)) = 0;
    end
end

 ROImap =  sum ( cat ( 3 , roiBW{:} ) , 3 ); % first term of cat (i.e., '3') points to element-wise alignement/stacking of arrays
tElapsedROICorr = toc(tStartROICorr);
disp(['    Time elapsed for ROI mask Correction was ', num2str(tElapsedROICorr/60),' minutes']);

%% Extract fluorescence traces for each neuron
% Loop through each neuron to get somatic fluo and neuropil fluo
for nn = 1:numNeurons
    
    [r,c]=find(roiBW{nn}~=0);
    tmpPixels = NaN(length(r),size(IMG,3));
    for i = 1:length(r)
        tmpPixels(i,:) = IMG(r(i),c(i),:);
    end
    rawFluo(:,nn) = nanmean(tmpPixels,1);
    
end

traceExtract_finish = toc(traceExtract_start);
fprintf('Trace extraction took %.1f minutes\n',traceExtract_finish/60)


end
%% FUNCTION: Sliding window dF/F calculation ZacBowen
% Inputs: F       - raw fluorecense
%         winsize - baseline window size
%         percent - lower percent of baseline values to average
% Output: dFF     - relative change in fluorescense in percent
function dFF = calculateDFF(F,winsize,percent)
nRois = size(F,1);
nFrames = size(F,2);
dFF = zeros(size(F));
%M = zeros(1, nrois);
%SD = zeros(1, nrois);
for j = 1 : nRois
  
    for k = 1 : nFrames
            lWin = max(1, k-winsize);
            rWin = min(k+winsize, nFrames);
            tWin = sort(F(j,(lWin:rWin)),'ascend');
%             percentWin = floor(percent/100*(rWin-lWin));
%             F0 = mean(tWin(1:percentWin));
            F0 = mean(rmoutliers(tWin,"percentiles",[0 percent]));
            
%             nDFF(j,k) = DFF(j,k)- twinm;   
            dFF(j,k) = 1 * (F(j,k)- F0) / F0; %SJG looks like ZB is multiplying 100, which is used in peak finding in other section of code, to correct can divide values by 100 as well to get true normalized DFF   
    end
    
    %twin = sort(nDFF(j,:));
    %pcentwin = floor(percent/100*(nframes));
    %M(j) = mean(twin(1:pcentwin));
    %SD(j) = std(twin(1:pcentwin));
    
end

end

%%

function [im_ch1, im_ch2, im_ch3] = nd2read(filename, varargin)
tic
finfo = nd2finfo(filename);
disp(['analyzing file structure used ', sprintf('%0.2f', toc), ' seconds'])

im_ch1 = zeros(finfo.img_width, finfo.img_height, 'uint16');
im_ch2 = zeros(finfo.img_width, finfo.img_height, 'uint16');
im_ch3 = zeros(finfo.img_width, finfo.img_height, 'uint16');
if finfo.ch_count == 4
  im_ch4 = zeros(finfo.img_width, finfo.img_height, 'uint16');
end

fid = fopen(filename, 'r');
fseek(fid, finfo.file_structure(strncmp('ImageDataSeq', ...
  {finfo.file_structure(:).nameAttribute}, 12)).dataStartPos, 'bof');

tic
% Image extracted from ND2 has image width defined by its first dimension.
if finfo.padding_style == 1
  if finfo.ch_count == 4
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch3(:, ii) = temp(1, :);
        im_ch1(:, ii) = temp(2, :);
        im_ch2(:, ii) = temp(3, :);
        im_ch4(:, ii) = temp(4, :);
        fseek(fid, 2, 'cof');
    end
  else
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch3(:, ii) = temp(1, :);
        im_ch1(:, ii) = temp(2, :);
        im_ch2(:, ii) = temp(3, :);
        fseek(fid, 2, 'cof');
    end
  end
else
  if finfo.ch_count == 4
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch1(:, ii) = temp(1, :);
        im_ch2(:, ii) = temp(2, :);
        im_ch3(:, ii) = temp(3, :);
        im_ch4(:, ii) = temp(4, :);
    end
  else
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch1(:, ii) = temp(1, :);
        im_ch2(:, ii) = temp(2, :);
        im_ch3(:, ii) = temp(3, :);
    end 
  end
end

fclose(fid);

im_ch1 = permute(im_ch1, [2 1]);
im_ch2 = permute(im_ch2, [2 1]);
im_ch3 = permute(im_ch3, [2 1]);
if finfo.ch_count == 4
    im_ch4 = permute(im_ch4, [2 1]);
end
if any(strcmpi(varargin, 'use_ch4'))
  im_ch3 = im_ch4;
end
  

disp(['reading complete image data used ', sprintf('%0.2f', toc), ' seconds'])
end
