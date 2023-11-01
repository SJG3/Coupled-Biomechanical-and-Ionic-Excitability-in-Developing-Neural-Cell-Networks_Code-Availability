%% SJG OF Lite            
clear all; clc;

% GOOD TO TEST ON
% /Users/sylvesterjgiii/Documents/LOSERT LAB/matlab of code lenny/FlowClusterTracking-master/HEK
% /Volumes/DS/M_P Data Dump/Losert_Lab/Microscopes/Multiscale/2021_2_12_ReN-lifeact_inGF/rescan_4spf_34msexp/2048HQ_10%_12/small_sec_for_OF_testing
% /Volumes/DS/M_P Data Dump/Losert_Lab/Microscopes/Multiscale/2021_3_3_ReN-lifeact_DiffD16/2048HQ_15%_0interval_6/testing/feature
% /Volumes/DS/M_P Data Dump/Losert_Lab/Microscopes/Multiscale/2021_2_12_ReN-lifeact_inGF/rescan_4spf_34msexp/2048HQ_10%_12/small_sec_for_OF_testing
% /Volumes/DS/M_P Data Dump/Losert_Lab/Microscopes/Multiscale/2021_3_5_ReN-lifeact_sparse_DiffD0/15%_34ms_57frames_5/features
% /Volumes/DS/M_P Data Dump/Losert_Lab/Microscopes/Multiscale/2021_4_7_ReNlifeact_puro_DiffD0_sparse_seed2daysprior_0.3e5_+1umCalbryte/561_att2V_15%laser-60ms_140int_inTyNa_1000frames_4x4_3/488-34ms_20s_10%_0.25att_1/512crop
%% (troubleshooting) create synthetic data
xstep = @(n) n .^ 1.3;
ystep = @(n) sin( n ); 
for i = 1:50
    % Create a logical image of a circle with specified
    % diameter, center, and image size.
    % First create the image.
    imageSizeX = 256;
    imageSizeY = 256;
    [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
    
    % Next create the circle in the image.
    centerX = 20 +  xstep(i); % Wherever you want.
    centerY = (imageSizeY / 1/4); %+ 50*ystep(1/5 * i); 
    radius = 15;
    circlePixels(:,:,i) = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
     
    circlePixels_noise(:,:,i) = (rand(imageSizeX)*0.1) + imgaussfilt(double(circlePixels(:,:,i)),5);

end 


% two balls runnnig away from one another
circlePixels2_noise = [flipud(fliplr(circlePixels_noise(:,:,1:50))),circlePixels_noise(:,:,1:50)];
% two balls running the same direction
circlePixels3_noise = [circlePixels_noise(:,:,1:50) , circlePixels_noise(:,:,1:50)]; 

for i = 1:50
    figure(365842); imagesc(circlePixels2_noise(:,:,i)); axis image;
    pause(0.001);
end
%% (troubleshooting) plot velocity of synthetic data
xy = [0 0; 0 0];
for i = 1:50
xy(i,:) = [i ,  (xstep(i)- xstep(i-1))/1 ];
end
figure(674);
plot(xy(:,1),xy(:,2)); xlabel ("frame"); ylabel("velocity (px/frame)");
%% (troubleshooting) load synthetic data for OF analysis
IMG = im2double(circlePixels2_noise);
IMGog = im2double(circlePixels2_noise);

%% START: LOAD Image                                                        

LOADER_STYLE = 1; 

% TifLink for loading
if LOADER_STYLE == 1 
  [fileX,pathX] = uigetfile('*.tif');
    disp("Opening File:" + fileX);
    im_path = fullfile(pathX,fileX);
    info = imfinfo(im_path);
    FileTif= im_path;

    %save the Tif into an array
    tic
    TifLink = Tiff(FileTif, 'r');
    for frame = 1: numel(info)
        TifLink.setDirectory(frame);
        IMGog(:,:,frame) = imresize((TifLink.read()),1);
    end
    TifLink.close();
         IMG= im2double(IMGog);
    disp(":::Completed TifLink Loading "+ numel(info)+" frames of the "+info(1).Width+"x"+info(1).Height+" file in " + secs2hms(toc));
end 

% ImRead for loading
if LOADER_STYLE == 2 
[fileX,pathX] = uigetfile({'*.gif;*.tif'});
    disp("Opening File:" + fileX);
    im_path = fullfile(pathX,fileX);
    info = imfinfo(im_path);
    FileTif= im_path;
    
    tic;
    for frame = 1: numel(info)
        IMGog(:,:,frame) = (imread(im_path,frame));
    end
    IMG= im2double(IMGog);
%     IMG = im2double(uint8(IMG));
    disp(":::Completed ImRead Loading "+ numel(info)+" frames of the "+info(1).Width+"x"+info(1).Height+" file in " + secs2hms(toc));
end 
%% Set Image Parameters                                                     
PxPer_dUnit = 4.7619; %1.3947; %4.7619;% 17.2 on the rescan multiscale ; 1.61 on widefield path; 60x on spinning disk 4.7619px/micron
distUnit = "\mum"; %"\mum"; %"\mum";
FramePer_tUnit = 0.5; %30;%1.1;%11.1
timeUnit =  "sec"; %"min" ; %"min" , "sec"
max_frame = size(IMG,3);

D_factor = [12 12]; % lower values = more arrows, closer 
S_factor = 12; % length/size of arrow. Larger = more noticible

OPTION.colormap_generate = 1;
if OPTION.colormap_generate == 1
    build_cmap = hot;
    cmap_GFP = [build_cmap(:,2) , build_cmap(:,1) , build_cmap(:,3) ];
    
    build_cmap = bone;
    cmap_GFP_norm = [build_cmap(:,1).*0.2 , build_cmap(:,2) , build_cmap(:,3).*0.2];
    
    clear build_cmap

    length_c = 100;
    c_red = [1, 0, 0];
    c_blue = [0, 0, 1];
    c_mid = [0.7, 0.7, 0.7];
    colors_r = [linspace(c_red(1),c_mid(1),length_c)', linspace(c_red(2),c_mid(2),length_c)', linspace(c_red(3),c_mid(3),length_c)'];
    colors_b = [linspace(c_blue(1),c_mid(1),length_c)', linspace(c_blue(2),c_mid(2),length_c)', linspace(c_blue(3),c_mid(3),length_c)'];
    cmap_HiLo = ([colors_r; flipud(colors_b)]);
    clear length_c c_red c_blue c_mid colors_r color_b;
end
%% Enhance/Process Image

OPTION.crop = 0;
OPTION.remove_outliers = 1;
OPTION.rescale  = 1;
OPTION.xy_smooth = 2;
OPTION.bg_subtract = 0;
OPTION.contrast = 1; 
OPTION.temporal_smooth = 0;     temporal_smooth_val = 3;

max_frame = size(IMG,3);
clear IMG; 
    for frame = 1: max_frame
        IMG(:,:,frame) = im2double(IMGog(:,:,frame));
    end
    
tic
if OPTION.crop ==1
    cropim = [];
    mean_im = mean(IMG,3);
    std_im = std(IMG,[],3);
%     crop_preview =(IMG(:,:,1));
    crop_preview = (std_im);
    crop_preview = filloutliers(crop_preview,'nearest','movmedian',8);
    crop_preview = imadjust(crop_preview); 
    figure(3742);
    [temp,rect] = imcrop(crop_preview); 
    for i = 1:max_frame
        cropim(:,:,i) = imcrop(IMG(:,:,i),round(rect));
    end
IMG = cropim;
clear cropim;
clear temp; clear crop_preview; clear rect;
end

if OPTION.remove_outliers == 1
    IMG = filloutliers((IMG),'center','movmedian',8,'ThresholdFactor',3);
%     IMGog_fill = filloutliers((IMGog),'center','movmedian',8,'ThresholdFactor',50);
%     clear im_filled
%     for i = 1:max_frame
%       blank_use = IMG(:,:,i);
%         %  compute local mdedian and standard deviation
%         local_med = medfilt2(blank_use,[5 5]);
%         local_std = stdfilt(blank_use,true(5));
% 
%         % define outliers as those 3 STDs from Median
%         outliers_logical = blank_use > local_med + 3*local_std; 
%         outliers_med_vals = outliers_logical .* local_med; 
% 
%         % store those outliers
%         [outl.x,outl.y] = find (outliers_med_vals);
% 
%         % fill outliers based on the median data
%         blank_fix(:,:,i) = blank_use;
%         for s = 1: size(outl.x,1) 
%             blank_fix(outl.x(s),outl.y(s),i)  = outliers_med_vals(outl.x(s),outl.y(s));
%         end
%     %     imagesc([blank_use>local_med + 4*local_std,blank_fix]); caxis([0 1]);
%     clear outliers_logical; clear outliers_med_vals;  
%     end
%     IMG = blank_fix;
end

if OPTION.rescale ~= 1
    IMG = imresize(IMG,OPTION.rescale);
%     PxPer_dUnit = PxPer_dUnit .* scale_val; 
end

if OPTION.xy_smooth > 0 
    IMG = (imgaussfilt(IMG,OPTION.xy_smooth));
end

if OPTION.bg_subtract ==1
    se = strel('disk',32);
    for i = 1:max_frame
        bg_detect = imopen(IMG(:,:,i),se);
        IMG(:,:,i) = IMG(:,:,i) - bg_detect;
    end
end 

if OPTION.contrast == 1
    vals = prctile(IMG(:),[50 99.9]);
    for i = 1:max_frame
         IMG(:,:,i) = imadjust(IMG(:,:,i), [vals(1) vals(2)],[0 1]);
    end 
end

if OPTION.temporal_smooth == 1
    t_win = gausswin(temporal_smooth_val);
    IMG = filter(t_win,1,IMG);
end

disp(":::Completed Image Preprocessing Enhancement in: " + secs2hms(toc));

figure(2345);
%  sliceViewer(IMG); colormap(cmap_GFP);
i = 5;
frame_prev = [IMG(:,:,i) , im2double(IMGog(:,:,i))];

 imagesc(frame_prev); colormap([0.2 0 0 ; cmap_GFP_norm]); 
%   caxis ([0 1]);
 axis image; colorbar;
%%      TEST: OF on first X frames                                                

OF_obj = opticalFlowFarneback('NumPyramidLevels',3,'FilterSize',50,'NeighborhoodSize',5,'PyramidScale',0.5); %found good results w/ filter siz 30, neighborhoodsiz 3
%  OF_obj = opticalFlowHS('Smoothness',5,'VelocityDifference',10e-10);
% OF_obj = opticalFlowLK('NoiseThreshold',10e-5);
% OF_obj = opticalFlowLKDoG('NoiseThreshold',10e-5,'GradientFilterSigma',10,'ImageFilterSigma',10); %10e-5, 10,10 work on non smoothed

val = 30
if max_frame >= val
    test_max_frame = val;
    disp("Displaying only first " +val+ " frames");
elseif max_frame < val
    test_max_frame = max_frame;
    disp("Displaying all " + max_frame+ " frames");
end


%testing Laplacian of image for smoothing
sig_ma = 0.2;
al_pha = 3;
for i = 1:max_frame
    B(:,:,i) = im2double(locallapfilt(single(IMG(:,:,i)),sig_ma,al_pha,'NumIntensityLevels', 10));
end 


for i = 1:test_max_frame
%  frame = [IMG(:,:,i),im2double(IMGog(:,:,i))];
 frame = [IMG(:,:,i),B(:,:,i)];
 flow_test(i) = estimateFlow(OF_obj,frame);
end

figure(66212);
for i = 1:test_max_frame
%     frame = [IMG(:,:,i),im2double(IMGog(:,:,i))];
    frame = [IMG(:,:,i),B(:,:,i)];
    imagesc(frame); colormap(cmap_GFP_norm); hold on; colorbar; caxis([0 1]);
    B2 = plot(flow_test(i),'DecimationFactor',D_factor,'ScaleFactor',S_factor);
        B2.LineStyleOrder = '-'; %type of line as '-',':','--', or '-.' 
        B2.ColorOrder = [0 0.5 1]; % color of the arrows
        B2.YLim: [0 size(frame,1)]; B2.XLim: [0 size(frame,2)];
        pbaspect([size(frame,2) size(frame,1) 1]);
    hold off;
    pause(0.05);
end

clear flow_test; clear frame;

reset(OF_obj);
%% Calculate Mask w/ and w/o Threshold                                      
tic
Mask_thresh = [];ThrshMask = [];vx_t = [];vy_t = [];flowOri_t = [];flowMag_t = [];

gauss_sig = 3;
pctle = [68 95 99.7 75 50] ; %percentiles represent 1,2,& 3 sigma

se = strel('disk',3,0);

for i = 1:max_frame 
%     frame = (IMG(:,:,i));
%     frame = (B(:,:,i));

    sig_ma = 0.2;
    al_pha = 10;
    frame = im2double(locallapfilt(single(IMG(:,:,i)),sig_ma,al_pha,'NumIntensityLevels', 10));
    
    mean_val_fr = mean(frame(:));
    std_val_fr = std(frame(:));
    ms_thresh_fr = mean_val_fr + (0*std_val_fr);
    pct_thresh_fr = prctile(frame(:),pctle); 
    median_val_fr = 3*median(frame(:));
    
    otsu_thresh_fr = graythresh(frame(:));
    multi_otsu_thresh_fr = multithresh(frame(:),5);
    

    if gauss_sig > OPTION.xy_smooth
        img_use = imgaussfilt(frame, (gauss_sig-OPTION.xy_smooth));
    else
        img_use = frame; 
    end

    ThrshMask = img_use > mean_val_fr;%std_val_fr(1); 
        % use with: mean_val_fr, std_val_fr, or ms_thresh_fr or mean_val, std_val, or pct_thresh  
%     ThrshMask = imquantize(IMG(:,:,i), multi_otsu_thresh_fr(5)); 
        % use with: otsu_thresh, multi_otsu_thresh or otsu_thresh_fr, multi_otsu_thresh_fr

    ThrshMask = im2double(ThrshMask);
    Mask_thresh(:,:,i) = ThrshMask;
end

diff_Mask = diff(IMG,[],3);
diff_Mask(:,:,max_frame) = diff_Mask(:,:,max_frame-1);
% diff_Mask = diff_Mask>0;
% diff_Mask = imdilate(diff_Mask,se);

diff_Mask_protrusive = (diff_Mask);
diff_Mask_protrusive = im2double(diff_Mask_protrusive);
diff_Mask_protrusive = imdilate(diff_Mask_protrusive,se);

diff_Mask_retractive = (diff_Mask)<0;
diff_Mask_retractive = im2double(diff_Mask_retractive);
diff_Mask_retractive = imdilate(diff_Mask_retractive,se);

diff_Mask_all = abs(diff_Mask)>0.025;
diff_Mask_all = im2double(diff_Mask_all);
% diff_Mask_all = imdilate(diff_Mask_all,se);
    %think of not using ABS and instead getting OF vectors for both
    %populations of expansive and contractive areas

% diff_Mask_all = imdilate(diff_Mask_all,se);
% diff_Mask_all = imerode(diff_Mask_all,se);
diff_Mask_all = bwareaopen((diff_Mask_all),50);

Mask_full_FOV = ones(size(IMG,1), size(IMG,2),max_frame);


figure(34523); 
ax(1) = subplot(2,3,1); imagesc(IMGog(:,:,3)); colormap(ax(1), cmap_GFP_norm); title ("Original Image");
ax(1) = subplot(2,3,2); imagesc(IMG(:,:,3)); caxis ([0 1]); colormap(ax(1), cmap_GFP_norm); title ("Processed Image");


ax(2) = subplot(2,3,4); imagesc(Mask_thresh(:,:,3)); caxis([-1 1]); colormap(ax(2), cmap_HiLo); title (["Thresholded ROI";" 'Mask thresh' "]);
ax(2) = subplot(2,3,5); imagesc(diff_Mask_all(:,:,3)); caxis([-1 1]); colormap(ax(2), cmap_HiLo); title(["ROI between frames";" 'diff Mask all' "]);
ax(2) = subplot(2,3,6); imagesc(Mask_full_FOV(:,:,3)); caxis([-1 1]); colormap(ax(2), cmap_HiLo); title(["Full FOV";"'Mask full FOV'"]);

disp(":::Completed Threshold Mask in: " + secs2hms(toc));
%% Optical Flow - Farneback method (30,3,3,68)        

tic
OF_obj = opticalFlowFarneback('NumPyramidLevels',7,'FilterSize',10,'NeighborhoodSize',1,'PyramidScale',0.5); %found good results w/ filter siz 30, neighborhoodsiz 3

Mask_Used = [];
test_1 = 1;
if test_1 == 1
    i = 1;
    while i < max_frame+1
        frame = (IMG(:,:,i));
%         frame = (im2double(IMGog(:,:,i)));
%         frame = (B(:,:,i));
        flow_all(i) = estimateFlow(OF_obj,(frame));

        ThrshMask = diff_Mask_all(:,:,i); %diff_Mask_all %Mask_thresh(:,:,i); %Mask_thresh(:,:,i); %Mask_full_FOV(:,:,i); %diff_Mask_all(:,:,i) .* Mask_thresh(:,:,i);
            %change ;ThrshMask to be where you want to collect OF vectors
        Mask_Used(:,:,i) = ThrshMask;

%         ThrshMask(ThrshMask == 0) = NaN;
%         ThrshMask(find(ThrshMask==0)) = NaN;

        vx_t(:,:,i) = flow_all(i).Vx .* ThrshMask; 
        vy_t(:,:,i) = flow_all(i).Vy .* ThrshMask;
        flow_masked(i) = opticalFlow(vx_t(:,:,i), vy_t(:,:,i));
        flow_masked_ycor(i) = opticalFlow(vx_t(:,:,i), vy_t(:,:,i) .* -1);

        %         flow_masked_ycor(i) = opticalFlow(flow_all(i).Vx .* ThrshMask , -1 .* flow_all(i).Vy .* ThrshMask);
        flowMag_t(:,:,i) = (flow_masked(i).Magnitude); 
            mag0= flow_masked(i).Magnitude~=0;
            mag0 = double(mag0);
            mag0(mag0==0) = NaN;
       
        flowOri_t(:,:,i) = (flow_masked_ycor(i).Orientation .* mag0);
        i = i+1;
    end
end

% Mask_Used(:,:,max_frame) = Mask_Used(:,:,max_frame-1);

test_2 = 0;
if test_2 == 1
    i = 1;
    while i < max_frame
        frame = (IMG(:,:,i));
        flow_all(i) = estimateFlow(OF_obj,(frame));

        ThrshMask = diff_Mask(:,:,i);
        ThrshMask(ThrshMask==0)=NaN;         
        flow_masked(i) = estimateFlow(OF_obj,(frame .* ThrshMask));
        i = i+1;
    end
end

disp(":::Completed Farnback Optical Flow In: " + secs2hms(toc));
reset(OF_obj);

sz = size(IMG,[1 2]);
vx_t(:,:,1) = NaN(sz); 
vy_t(:,:,1) = NaN(sz);
flowOri_t(:,:,1) = NaN(sz);
flowMag_t(:,:,1) = NaN(sz);
%%      TEST: Preview of OF video                                                      
    Prev_Vid_OF2 = 1;
  
if Prev_Vid_OF2 == 1
    for i = 1: max_frame
        figure(2374);
        subplot(2,2,1); image(IMGog(:,:,i)); axis image; title("Original Image");
        subplot(2,2,2); imagesc(IMG(:,:,i)); axis image; title("Processed Image");
        subplot(2,2,3); imagesc(Mask_Used(:,:,i)); axis image; title("Mask");
        subplot(2,2,4); imagesc(IMG(:,:,i)); axis image; colormap([cmap_GFP]);
            hold on
            p1 = plot(flow_masked(i),'DecimationFactor',D_factor,'ScaleFactor',S_factor); axis tight;
            
                p1.LineStyleOrder = '-'; %type of line as '-',':','--', or '-.' 
                p1.ColorOrder = [1 0 1]; % color of the arrows
                p1.YLim: [0 size(IMG,1)]; p1.XLim: [0 size(IMG,2)];
                pbaspect([size(IMG,2) size(IMG,1) 1]);
              
                hold off
            pause(10^-3);
    end
end
%% Process data, on frames: magnitude, orientation                          
OF_direc_frame_bins_counts = [];
n_bins = 36;

hist_edges = [-3.14 : pi/((1+n_bins)/2) : 3.14];


for i = 1:max_frame;   
    holder = flow_masked(i).Orientation; %flowOri_t(:,:,i);
    holder = reshape(holder(~isnan(holder)),[],1);
    [OF_direc_frame_bins_counts(:,i),bins_lab] = histcounts(holder,hist_edges, 'Normalization','Probability');

    holder = flowMag_t(:,:,i); %flow_masked(i).Magnitude; % flowMag_t(:,:,i);
    holder = reshape(holder(~isnan(holder)),[],1);

    OF_median_speeds(i) = median(holder) * FramePer_tUnit * 1/PxPer_dUnit;
    OF_mean_speeds(i) = mean(holder) * FramePer_tUnit * 1/PxPer_dUnit;
    OF_speeds_std(i) = std(holder) * FramePer_tUnit * 1/PxPer_dUnit;
end

for i = 1 : max_frame
    re_vx = reshape( (flow_masked(i).Vx),[],1); 
    re_vx = re_vx(~isnan(re_vx)); 
    av_vx = sum(re_vx) / size(re_vx,1);
    median_vx = median(re_vx); 

    re_vy = reshape(flow_masked(i).Vy,[],1);
    re_vy = re_vy(~isnan(re_vy)); 
    av_vy = sum(re_vy) / size(re_vy,1);
    median_vy = median(re_vy);
    
    av_vmag(i) = sqrt(av_vx.^2 + av_vy.^2);
    median_mag(i) = sqrt(median_vx.^2 + median_vy.^2);
end

holder_ori = reshape(flowOri_t(~isnan(flowOri_t)),[],1); %store orientation of OF vectors from all frames w/o NaNs
holder_mag = reshape(flowMag_t(~isnan(flowMag_t)),[],1); %store magnitude of OF vectors from all frames w/o NaNs
%% Initialize Figure Generation                                             
rnd_Var = randi(1000); 
%% FIGURE: pt1                                                              

figure (rnd_Var);

subdiv_2pi = 8;
subplot(2,3,1); ph = polarhistogram(holder_ori,36,'Normalization','Probability'); 
    title ("Integrated OF Directional Probability"); 
    theta_ticks = [0:rad2deg(2*pi/subdiv_2pi):360];
    thetaticks(theta_ticks);
    [N,D] = rat(wrapTo180(theta_ticks) * 1/(180));
    tickss = N+"\pi/"+D;
    tickss=strrep(tickss,'/1','');
    tickss=strrep(tickss,'0\pi','0');
    tickss=strrep(tickss,'1','');
    thetaticklabels(tickss);

subplot(2,3,2); histogram(holder_mag * FramePer_tUnit * 1/PxPer_dUnit,36,'Normalization','Probability');
    title(["Integrated OF Speed Probability";"Mean Speed: "+ (nanmean(flowMag_t(:))*FramePer_tUnit * 1/PxPer_dUnit) + " "+distUnit+"/"+timeUnit ; "Median Speed: "+ (nanmedian(flowMag_t(:))*FramePer_tUnit * 1/PxPer_dUnit) + " "+distUnit+"/"+timeUnit]);
    axis square; ylabel(["Probability"]); xlabel(["Speed("+distUnit+"/"+timeUnit+")"]); 

subplot(2,3,3); imagesc(sum(Mask_Used,3)); 
    title(["Max Projection of Threshold Mask"]);
    axis image; colormap turbo;  caxis([0 size(IMG,3)]);

    
subplot(2,3,4);
        sz_half = size(OF_direc_frame_bins_counts,1)/2;
        padded_ori = [ OF_direc_frame_bins_counts(sz_half:end,:) ;  OF_direc_frame_bins_counts ; OF_direc_frame_bins_counts(1:sz_half,:)];
        padded_ori = imgaussfilt([padded_ori],1);
     imagesc(padded_ori(sz_half+2:sz_half+2+n_bins,:)); 
        axis square; 
        title(["OF Directional Probability Over Time", n_bins + " bins of size "+360/n_bins+" degree on Y-axis"]); xlabel("Frame"); ylabel("Direction");
        colormap (turbo); %caxis([0 1]);
        cb_ = colorbar;
        ylabel(cb_, 'Probability');
        ytick_inc = [1:(n_bins)./subdiv_2pi*2:n_bins+1];
        yticks(ytick_inc);
        tickss = (-(ytick_inc-1).*360/n_bins)+180;
        [N,D] = rat(wrapTo180(tickss) * 1/(180));
        tickss = N+"\pi/"+D;
        tickss=strrep(tickss,'/1','');
        tickss=strrep(tickss,'0\pi','0');
        tickss=strrep(tickss,'1',''); 
        yticklabels( tickss );

subplot(2,3,5); 
%     plot(OF_mean_speeds);
    shadedErrorBar([2:length(OF_median_speeds)],OF_median_speeds(2:end),OF_speeds_std(2:end),'lineProps','b','transparent',0,'patchSaturation',0.05); 
    axis square;
    title("Median OF Speed Per Frame"); xlabel("Frame"); ylabel("Mean OF Magnitude ("+distUnit+"/"+timeUnit+")");
    axis square; xlim([0 max_frame]); ylim([0 inf]);

subplot(2,3,6); plot(median_mag(2:end)); 
    title("Magnitude of Median OF Velocity"); axis square; xlabel("Frame"); xlim([0 max_frame]); ylim([0 inf]);
%%      Figures Movies display original, OF vectors, magnitude, and orientation  
figure(rnd_Var+2);
step_sz = 10;

i = 1;
while i < max_frame   
    ax(1) = subplot(2,3,1); imagesc((IMGog(:,:,i))); axis image; colormap(ax(1) , cmap_GFP); title("Original Image");
    ax(2) = subplot(2,3,2); imagesc((IMG(:,:,i))); axis image; title("Processed Image"); colormap(ax(2) , cmap_GFP); 
    ax(3) = subplot(2,3,3); imagesc(zeros(size(IMG,1),size(IMG,2))); axis image; hold on;  
         B3 = plot(flow_masked(i),'DecimationFactor',D_factor,'ScaleFactor',S_factor);
            B3.LineStyleOrder = '-'; %type of line as '-',':','--', or '-.' 
            B3.ColorOrder = [0 0.7 1]; % color of the arrows
            B3.YLim: [0 size(IMG,1)]; B3.XLim: [0 size(IMG,2)];
            pbaspect([size(IMG,2) size(IMG,1) 1]);
         axis equal; hold off; colormap(ax(3) ,[0 0 0 ; gray]); title("OF Vectors");

    ax(4) = subplot(2,3,4); imagesc(Mask_Used(:,:,i)); axis image; title("Masked ROI");
    ax(5) = subplot(2,3,5); imagesc(flow_masked(i).Magnitude); caxis([0 max(flowMag_t,[],'all')]); axis image; colormap(ax(5), [0 0 0 ; jet]); title("Masked OF Magnitude"); cb_ = colorbar; ylabel(cb_, "Magnitude"); 
    ax(6) = subplot(2,3,6); imagesc(flow_masked_ycor(i).Orientation); caxis([-pi pi]); axis image; colormap(ax(6) ,[0 0 0 ; hsv]); title("Masked OF Orienation"); cb_ = colorbar; ylabel(cb_, "Directionality (Rad)"); 
    i = i+1;
    pause(0.0001);
end 

%think about adding recipricol charts that are shownn in main figures to
%this one. I.E. OF directional probabity, histogram of speeds

%think about integrating both magnitude and directionality into one single
%chart
%% SAVE or View video OF    

% Initialize video
% myVideo = VideoWriter(pathX+"OF_masked_Overlay_IMGprocessed_of"+fileX,'MPEG-4'); %open video file
% % myVideo = VideoWriter(pathX+"testing_1_2_3",'MPEG-4'); %open video file
% myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
% open(myVideo)


figure(66212);
for i = 1:max_frame
    
    frame = [IMG(:,:,i),im2double(IMGog(:,:,i))];
    ax(1) = subplot(2,2,1); 
    frame = [IMG(:,:,i)];
    imagesc(frame); colormap(cmap_GFP_norm); hold on; colorbar; caxis([0 1]); title("IMG & OF Vectors; Frame: "+i);
    B2 = plot(flow_masked(i),'DecimationFactor',D_factor/3,'ScaleFactor',S_factor/6);
        B2.LineStyleOrder = '-'; %type of line as '-',':','--', or '-.' 
        B2.ColorOrder = [0 0.5 1]; % color of the arrows
        B2.YLim = [0 size(frame,1)]; B2.XLim: [0 size(frame,2)];
        B2.YTick = []; B2.XTick = [];
        pbaspect([size(frame,2) size(frame,1) 1]);
    hold off;
    
    ax(5) = subplot(2,2,2); imagesc(flow_masked(i).Magnitude); caxis([0 max(flowMag_t,[],'all')]); axis image; colormap(ax(5), [0.2 0 0 ; turbo ]); title("Masked OF Magnitude"); cb_ = colorbar; ylabel(cb_, "Magnitude"); 
    ax(6) = subplot(2,2,3); imagesc(flowOri_t(:,:,i)); caxis([-pi pi]); axis image; colormap(ax(6) ,[0.2 0 0 ; hsv]); title("Masked OF Orienation"); cb_ = colorbar; ylabel(cb_, "Directionality (Rad)");
%     ax(6) = subplot(2,2,4); imagesc(flow_masked(i).Orientation); caxis([-pi pi]); axis image; colormap(ax(6) ,[0 0 0 ; hsv]); title("Masked OF Orienation"); cb_ = colorbar; ylabel(cb_, "Directionality (Rad)");

    
%     subplot(2,2,4); polarhistogram(flowOri_t(:,:,i),36,'Normalization','Probability'); 
     subplot(2,2,4); polarhistogram(flow_masked_ycor(i).Orientation,36,'Normalization','Probability');
        title ("OF Directional Probability, Frame "+i); 
        theta_ticks = [0:rad2deg(2*pi/subdiv_2pi):360];
        thetaticks(theta_ticks);
        [N,D] = rat(wrapTo180(theta_ticks) * 1/(180));
        tickss = N+"\pi/"+D;
        tickss=strrep(tickss,'/1','');
        tickss=strrep(tickss,'0\pi','0');
        tickss=strrep(tickss,'1','');
        thetaticklabels(tickss);

    pause(0.0001);
%     frame_vid = getframe(gcf); %get frame for video
%     writeVideo(myVideo, frame_vid); %save frame to video
end

% close(myVideo)


clear flow_test; clear frame;

reset(OF_obj);
%% preview of entire movie
figure(rnd_Var+662); 
for i = 1:max_frame
    
    frame = [IMG(:,:,i),im2double(IMGog(:,:,i))] ;
    ax(1) = subplot(1,1,1); 
    frame = [im2double(IMG(:,:,i))];
    imagesc(frame); colormap([0.2 0 0 ; cmap_GFP_norm]); %colormap(cmap_GFP_norm); 
    hold on; colorbar; caxis([0 1]); title("IMG & OF Vectors; Frame: "+i);
    B2 = plot(flow_masked(i),'DecimationFactor',D_factor/2,'ScaleFactor',S_factor/6);
        B2.LineStyleOrder = '-'; %type of line as '-',':','--', or '-.' 
        B2.ColorOrder = [1 0 0]; % color of the arrows
        B2.YLim = [0 size(frame,1)]; B2.XLim: [0 size(frame,2)];
        B2.YTick = []; B2.XTick = [];
        pbaspect([size(frame,2) size(frame,1) 1]);
    hold off;
  

    pause(0.05);

end

%%  FFT for every orientation
n = size(OF_direc_frame_bins_counts,1);
    for i = 1:n 
    fftransform = fft(OF_direc_frame_bins_counts(i,:)');
    L = size(IMG,3);

    P2 = abs(fftransform/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    P3(:,i) = squeeze(P1);
    
    Fs = 0.5;
    
    f = Fs*(0:(L/2))/L;
    end
    
figure(rnd_Var+653);
plot(f,P3(:,:));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
    
    
    
figure(rnd_Var+652);
imagesc([P3';P3(:,1)']); colormap (turbo);
segm = 10; 
xticks([1 : size(f,2)/segm : size(f,2) , size(f,2) ]);
xticklabels({f(1) : f(end)/segm : f(end)});
xlabel("frequency");
ylabel("Orientation");
yticks(ytick_inc);
yticklabels( tickss );
%%
Y = fft2(OF_direc_frame_bins_counts);
imagesc(abs(fftshift(Y)))

%% Figure: Projection of top % of Magnitude over entire film

top_pct = 15 ;

figure(345+rnd_Var);
for i = 1:max_frame
    holder = reshape(flow_masked(i).Magnitude,[],1);  
    cutoff = prctile(holder,100-top_pct);
    Top_pct_mag(:,:,i) = flow_masked(i).Magnitude > cutoff;
end 

imagesc(sum(Top_pct_mag,3)); colormap([0 0 0 ; turbo]); title("Projection of top "+ top_pct+"% of Magnitude over entire film"); 
axis('image');


%% Function
function ax = phasebar(varargin) 
    % phasebar places a circular donunt-shaped colorbar for phase 
    % from -pi to pi or -180 degrees to 180 degrees. 
    % 
    %% Syntax
    % 
    %  phasebar
    %  phasebar(...,'location',Location) 
    %  phasebar(...,'size',Size) 
    %  phasebar('deg') 
    %  phasebar('rad') 
    %  ax = phasebar(...) 
    % 
    %% Description 
    % 
    % phasebar places a donut-shaped colorbar on the current axes. 
    %
    % phasebar(...,'location',Location) specifies the corner (e.g., 'northeast' or 'ne') 
    % of the current axes in which to place the phasebar. Default location is the upper-right or 'ne' 
    % corner. 
    %
    % phasebar(...,'size',Size) specifies a size fraction of the current axes.  Default is 0.3. 
    %
    % phasebar('deg') plots labels at every 90 degrees. 
    %
    % phasebar('rad') plots labels at every pi/2 radians. 
    %
    % ax = phasebar(...) returns a handle ax of the axes in which the new axes are plotted. 
    % 
    %% Example
    % 
    % Z = 200*peaks(900); 
    % Zw = phasewrap(Z,'degrees'); 
    % imagesc(Zw) 
    % phasemap(12)
    % phasebar('location','se')
    % 
    %% Author Info
    % This function was written by Chad A. Greene of the University of Texas 
    % at Austin's Institute for Geophysics (UTIG), May 2016. 
    % This function includes Kelly Kearney's plotboxpos function as a subfunction. 
    % 
    % If the phasemap function is useful for you, please consider citing our 
    % paper about it: 
    % 
    % Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 
    % 2016. True colors of oceanography: Guidelines for effective and accurate 
    % colormap selection. Oceanography 29(3):9?13. 
    % http://dx.doi.org/10.5670/oceanog.2016.66
    % 
    % See also colorbar and phasemap. 
    %% Set Defaults: 
    usedegrees = false; 
    axsize = 0.3; 
    location = 'northeast'; 
    % Try to automatically determine if current displayed data exist and are in radians or degrees: 
    if max(abs(caxis))>pi
       usedegrees = true;
    else
       usedegrees = false; 
    end
    % If no data are already displayed use radians: 
    if isequal(caxis,[0 1])
       usedegrees = false; 
    end
    %% Parse inputs: 
    tmp = strncmpi(varargin,'location',3); 
    if any(tmp) 
       location = varargin{find(tmp)+1}; 
    end
    tmp = strncmpi(varargin,'size',3); 
    if any(tmp) 
       axsize = varargin{find(tmp)+1}; 
       assert(isscalar(axsize)==1,'Input error: axis size must be a scalar greater than zero and less than one.') 
       assert(axsize>0,'Input error: axis size must be a scalar greater than zero and less than one.') 
       assert(axsize<1,'Input error: axis size must be a scalar greater than zero and less than one.') 
    end
    if any(strncmpi(varargin,'radians',3)); 
       usedegrees = false; 
    end
    if any(strncmpi(varargin,'degrees',3)); 
       usedegrees = true; 
    end
    %% Starting settings: 
    currentAx = gca; 
    cm = colormap; 
    pos = plotboxpos(currentAx); 
    xcol = get(currentAx,'XColor'); 
    % Delete old phasebar if it exists: 
    try
       oldphasebar = findobj(gcf,'tag','phasebar'); 
       delete(oldphasebar); 
    end
    %% Created gridded surface: 
    innerRadius = 10; 
    outerRadius = innerRadius*1.618; 
    [x,y] = meshgrid(linspace(-outerRadius,outerRadius,300));
    [theta,rho] = cart2pol(x,y); 
    % theta = rot90(-theta,3); 
    theta(rho>outerRadius) = nan; 
    theta(rho<innerRadius) = nan; 
    if usedegrees
       theta = theta*180/pi; 
    end
    %% Plot surface: 
    ax = axes; 
    pcolor(x,y,theta)
    shading interp 
    hold on
    % Plot a ring: 
    [xc1,yc1] = pol2cart(linspace(-pi,pi,360),innerRadius); 
    [xc2,yc2] = pol2cart(linspace(-pi,pi,360),outerRadius); 
    plot(xc1,yc1,'-','color',xcol,'linewidth',.2); 
    plot(xc2,yc2,'-','color',xcol,'linewidth',.2); 
    axis image off
    colormap(cm) 
    if usedegrees
       caxis([-180 180]) 
    else
       caxis([-pi pi]) 
    end
    %% Label: 
    [xt,yt] = pol2cart((-1:2)*pi/2+pi/2,innerRadius); 
    if usedegrees
       text(xt(1),yt(1),'0\circ','horiz','right','vert','middle'); 
       text(xt(2),yt(2),'90\circ','horiz','center','vert','top'); 
       text(xt(3),yt(3),'180\circ','horiz','left','vert','middle'); 
       text(xt(4),yt(4),'-90\circ','horiz','center','vert','bottom'); 
    else
       text(xt(1),yt(1),'0','horiz','right','vert','middle'); 
       text(xt(2),yt(2),'\pi/2','horiz','center','vert','top'); 
       text(xt(3),yt(3),'\pi','horiz','left','vert','middle'); 
       text(xt(4),yt(4),'-\pi/2','horiz','center','vert','bottom'); 
    end
    %% Set position of colorwheel: 
    switch lower(location)
       case {'ne','northeast'} 
          set(ax,'position',[pos(1)+(1-axsize)*pos(3) pos(2)+(1-axsize)*pos(4) axsize*pos(3) axsize*pos(4)]); 
       case {'se','southeast'} 
          set(ax,'position',[pos(1)+(1-axsize)*pos(3) pos(2) axsize*pos(3) axsize*pos(4)]); 

       case {'nw','northwest'} 
          set(ax,'position',[pos(1) pos(2)+(1-axsize)*pos(4) axsize*pos(3) axsize*pos(4)]); 

       case {'sw','southwest'} 
          set(ax,'position',[pos(1) pos(2) axsize*pos(3) axsize*pos(4)]); 

       otherwise
          error('Unrecognized axis location.') 
    end

    %% Clean up 
    set(ax,'tag','phasebar')
    % Make starting axes current again: 
    axes(currentAx); 
    uistack(ax,'top'); 
    if nargout==0 
       clear ax
    end
end
%% FUNCTION: Graph With Shaded Error Bars == varargout=shadedErrorBar(x,y,errBar,varargin) 
 function varargout=shadedErrorBar(x,y,errBar,varargin)
    % generate continuous error bar area around a line plot
    %
    % function H=shadedErrorBar(x,y,errBar, ...)
    %
    % Purpose 
    % Makes a 2-d line plot with a pretty shaded error bar made
    % using patch. Error bar color is chosen automatically.
    %
    %
    % Inputs (required)
    % x - vector of x values [optional, can be left empty]
    % y - vector of y values or a matrix of n observations by m cases
    %     where m has length(x);
    % errBar - if a vector we draw symmetric errorbars. If it has a size
    %          of [2,length(x)] then we draw asymmetric error bars with
    %          row 1 being the upper bar and row 2 being the lower bar
    %          (with respect to y -- see demo). ** alternatively ** 
    %          errBar can be a cellArray of two function handles. The 
    %          first defines statistic the line should be and the second 
    %          defines the error bar.
    %
    % Inputs (optional, param/value pairs)
    % 'lineProps' - ['-k' by default] defines the properties of
    %             the data line. e.g.:    
    %             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
    % 'transparent' - [true  by default] if true, the shaded error
    %               bar is made transparent. However, for a transparent
    %               vector image you will need to save as PDF, not EPS,
    %               and set the figure renderer to "painters". An EPS 
    %               will only be transparent if you set the renderer 
    %               to OpenGL, however this makes a raster image.
    % 'patchSaturation'- [0.2 by default] The saturation of the patch color.
    %
    %
    %
    % Outputs
    % H - a structure of handles to the generated plot objects.
    %
    %
    % Examples:
    % y=randn(30,80); 
    % x=1:size(y,2);
    %
    % 1)
    % shadedErrorBar(x,mean(y,1),std(y),'lineProps','g');
    %
    % 2)
    % shadedErrorBar(x,y,{@median,@std},'lineProps',{'r-o','markerfacecolor','r'});
    %
    % 3)
    % shadedErrorBar([],y,{@median,@(x) std(x)*1.96},'lineProps',{'r-o','markerfacecolor','k'});
    %
    % 4)
    % Overlay two transparent lines:
    % clf
    % y=randn(30,80)*10; 
    % x=(1:size(y,2))-40;
    % shadedErrorBar(x,y,{@mean,@std},'lineProps','-r','transparent',1);
    % hold on
    % y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
    % shadedErrorBar(x,y,{@mean,@std},'lineProps','-b','transparent',1);
    % hold off
    %
    %
    % Rob Campbell - November 2009
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parse input arguments
    narginchk(3,inf)
    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('lineProps', '-k', @(x) ischar(x) | iscell(x));
    if (sum( size(ver('MATLAB'))) > 0  )
      params.addParameter('transparent', true, @(x) islogical(x) || x==0 || x==1);
    elseif (sum( size(ver('Octave'))) > 0  )
      params.addParameter('transparent', false, @(x) islogical(x) || x==0 || x==1);
    end
    params.addParameter('patchSaturation', 0.2, @(x) isnumeric(x) && x>=0 && x<=1);
    params.parse(varargin{:});
    %Extract values from the inputParser
    lineProps =  params.Results.lineProps;
    transparent =  params.Results.transparent;
    patchSaturation = params.Results.patchSaturation;
    if ~iscell(lineProps), lineProps={lineProps}; end
    %Process y using function handles if needed to make the error bar dynamically
    if iscell(errBar) 
        fun1=errBar{1};
        fun2=errBar{2};
        errBar=fun2(y);
        y=fun1(y);
    else
        y=y(:).';
    end
    if isempty(x)
        x=1:length(y);
    elseif sum( size(ver('MATLAB'))) > 0 
        x=x(:).';
    end
    %Make upper and lower error bars if only one was specified
    if length(errBar)==length(errBar(:))
        errBar=repmat(errBar(:)',2,1);
    else
        s=size(errBar);
        f=find(s==2);
        if isempty(f), error('errBar has the wrong size'), end
        if f==2, errBar=errBar'; end
    end
    % Check for correct x, errbar formats
    x_size = size(x);
    if (length(x) ~= length(errBar) && sum( size(ver('MATLAB'))) > 0 )
        error('length(x) must equal length(errBar)')
    elseif( ( length(x) ~= length(errBar) && checkOctave_datestr(x) == false ) ...
                && sum( size(ver('Octave'))) > 0  )
        error('length(x) must equal length(errBar) or x must have valid datestr')
    end

    %Log the hold status so we don't change
    initialHoldStatus=ishold;
    if ~initialHoldStatus, hold on,  end
    H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation);
    if ~initialHoldStatus, hold off, end
    if nargout==1
        varargout{1}=H;
    end
    function H = makePlot(x,y,errBar,lineProps,transparent,patchSaturation)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine host application
        if (sum( size(ver('MATLAB'))) > 0  )
          hostName = 'MATLAB';
        elseif (sum(size(ver('Octave'))) > 0)
          hostName = 'Octave';
        end % if

        % Plot to get the parameters of the line
        if hostName == 'MATLAB'
          H.mainLine=plot(x,y,lineProps{:});

        elseif hostName == 'Octave'
          boolxDatestr = checkOctave_datestr(x);
          if boolxDatestr
            x = datenum(x);
            x = x(:).';
            H.mainLine=plot(x,y,lineProps{:});
            datetick(gca);
          else
            H.mainLine=plot(x,y,lineProps{:});
          end
        end
        % Tag the line so we can easily access it
        H.mainLine.Tag = 'shadedErrorBar_mainLine';
        % Work out the color of the shaded region and associated lines.
        % Here we have the option of choosing alpha or a de-saturated
        % solid colour for the patch surface.
        mainLineColor=get(H.mainLine,'color');
        edgeColor=mainLineColor+(1-mainLineColor)*0.55;
        if transparent
            faceAlpha=patchSaturation;
            patchColor=mainLineColor;
        else
            faceAlpha=1;
            patchColor=mainLineColor+(1-mainLineColor)*(1-patchSaturation);
        end
        %Calculate the error bars
        uE=y+errBar(1,:);
        lE=y-errBar(2,:);
        %Make the patch (the shaded error bar)
        yP=[lE,fliplr(uE)];
        xP=[x,fliplr(x)];
        %remove nans otherwise patch won't work
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];

        if isdatetime(x) && strcmp(hostName,'MATLAB')
          H.patch=patch(datenum(xP),yP,1);
        else
          H.patch=patch(xP,yP,1);
        end
        set(H.patch,'facecolor',patchColor, ...
            'edgecolor','none', ...
            'facealpha',faceAlpha, ...
            'HandleVisibility', 'off', ...
            'Tag', 'shadedErrorBar_patch')
        %Make pretty edges around the patch. 
        H.edge(1)=plot(x,lE,'-');
        H.edge(2)=plot(x,uE,'-');
        set([H.edge], 'color',edgeColor, ...
          'HandleVisibility','off', ...
          'Tag', 'shadedErrorBar_edge')
        % Ensure the main line of the plot is above the other plot elements
        if hostName == 'MATLAB'
          if strcmp(get(gca,'YAxisLocation'),'left') %Because re-ordering plot elements with yy plot is a disaster
            uistack(H.mainLine,'top')
          end
        elseif hostName == 'Octave'
          % create the struct from scratch by temp.
          H = struct('mainLine', H.mainLine, ...
          'patch', H.patch, ...
          'edge', H.edge);
        end
    end
%%
    function boolDate = checkOctave_datestr(x)
      %% Simple try/catch for casting datenums, requireing valid datestr
      boolDate = true;
      try
        datenum(x)
      catch
        boolDate = false;
      end
    end
end
%% FUNCTION: Timer == time_string=secs2hms(time_in_secs)
% https://www.mathworks.com/matlabcentral/fileexchange/22817-seconds-to-hours-minutes-seconds
 function time_string=secs2hms(time_in_secs)
    %SECS2HMS - converts a time in seconds to a string giving the time in hours, minutes and second
    %Usage TIMESTRING = SECS2HMS(TIME)]);
    %Example 1: >> secs2hms(7261)
    %>> ans = 2 hours, 1 min, 1.0 sec
    %Example 2: >> tic; pause(61); disp(['program took ' secs2hms(toc)]);
    %>> program took 1 min, 1.0 secs

    time_string='';
    nhours = 0;
    nmins = 0;
    if time_in_secs >= 3600
        nhours = floor(time_in_secs/3600);
        if nhours > 1
            hour_string = ' hours, ';
        else
            hour_string = ' hour, ';
        end
        time_string = [num2str(nhours) hour_string];
    end
    if time_in_secs >= 60
        nmins = floor((time_in_secs - 3600*nhours)/60);
        if nmins > 1
            minute_string = ' mins, ';
        else
            minute_string = ' min, ';
        end
        time_string = [time_string num2str(nmins) minute_string];
    end
    nsecs = time_in_secs - 3600*nhours - 60*nmins;
    time_string = [time_string sprintf('%2.1f', nsecs) ' secs'];
end