%% 
clear
clc 
%Sylvester James Gates III, Losert Lab, University of Maryland College Park
%2021-23
%This code is used to find peaks/spikes from a deltaF/F matrix and then can
%be used to output figures which include bar graphs for the active cells in
%population, percentage of active cells, and percetage of spikes in active
%cells.
%% conditions
tot_time_m = 1; %total time in minutes
fps = 5;            % FPS of video used

%% Use FindPeaks to automate detecting spikes, spike width at half prominence, and height

%create empty vectors which will store the spikes widths and prominence
widths_spikes = NaN(size(DFF,1));
prom_spikes = NaN(size(DFF,1));

%for every cell, use findpeaks with moving mean (10) to find peaks with
%prominence greater than 20. store the spikes widths/ppriminence of each cell on different
%row, with each spikes width/prominence on a different column
zen = DFF*0;

spike_width = 10; % wanted spike width in seconds
    %initally use 30seconds, then 10seconds, then 3seconds
for n = 1:size(DFF,1)
    [pks,locs,w,p] = findpeaks(movmean(DFF(n,:)',1),'MinPeakProminence',0.0075,'MinPeakWidth',fps*0.5,'MaxPeakWidth',fps*spike_width);
    number_of_spikes(n) = size(locs,1);
    widths_spikes(n,1:size(w',2)) = [w'];
    prom_spikes(n,1:size(w',2)) = [p'];
    
    for n_l = 1:size(locs,1)
        zen(n,locs(n_l)) = 1;
    end
end 

%use the stored widths and prominence to get a mean spike width for each
%cell and prominence
widths_spikes_mean_per_cell = nanmean(widths_spikes,2);
prom_spikes_mean_per_cell = nanmean(prom_spikes,2);

%call active cells those with greater than 3 spikes
active_cells_spikes = number_of_spikes(number_of_spikes>3); 
active_cells_mean_spikes = mean(active_cells_spikes);
%% Percentage cells active (ie > 3 spikes)
total_cells = size(DFF,1);
pct = size(active_cells_spikes,2) / total_cells
tot_spike_num = sum(number_of_spikes) 
tot_active_spike_num = sum(active_cells_spikes)
spikes_per_min = tot_spike_num/tot_time_m
%% Histogram of count as spike per min
figure(5347259);
histogram (number_of_spikes/tot_time_m , 'Normalization','probability', 'BinWidth',0.1, 'FaceAlpha',0.15);
hold on; 
%% FIGURES
figure(26543245);
subplot(1,2,1); scatter(prom_spikes_mean_per_cell,1:size(DFF,1),'.'); xlabel("Mean Spike Prominence (Intensity)"); ylabel("Cell Index"); hold on;
subplot(1,2,2); scatter(widths_spikes_mean_per_cell,1:size(DFF,1),'.');  hold on; 
% errorbar(widths_spikes_mean_per_cell,1:size(DFF,1),ans,'horizontal','.');
xlabel("Mean Spike Width (frames)"); ylabel("Cell Index");
%%      NOT NEEDED - TESTING for spike finder
figure(45432); 
nom = 35; 
 
for nom = 1:size(DFF,1)%[1:10:size(DFF,1)]
    % findpeaks(detrend(movmean(DFF(nom,:)',10),3),'MinPeakProminence',10,'Annotate','extents');title("cell "+nom);
    findpeaks(movmean(DFF(nom,:)',fps),'MinPeakProminence',0.0075,'MinPeakWidth',fps*0.5,'MaxPeakWidth',fps*spike_width,'Annotate','extents');
%     findpeaks(movmean(DFF(nom,:)',fps),'MinPeakProminence', mean(reshape(abs(diff(DFF_smooth)),[],1))* 3 ,'MinPeakWidth',fps*0.5,'MaxPeakWidth',fps*spike_width,'Annotate','extents');
%     findpeaks(movmean(DFF(nom,:)',fps),'MinPeakProminence', std(reshape(DFF_smooth(:,:),[],1))*3,'MinPeakWidth',fps*0.5,'MaxPeakWidth',fps*spike_width,'Annotate','extents');
   
ylim([- max(DFF(:)) ,  max(DFF(:))]);
%     hold on
pause(0.2);
end 
%%      NOT NEEDED TESTING is detrending needed? 
nom = 52;
y = detrend(movmean(DFF(nom,:)',10),3);

figure(3453);
subplot(1,3,1); plot(y); ylim([- max(DFF(:)) ,  max(DFF(:))]);

subplot(1,3,2); plot(movmean(DFF(nom,:)',10));ylim([- max(DFF(:)) ,  max(DFF(:))]);

subplot(1,3,3); plot(smoothdata(DFF(nom,:)','sgolay',10));ylim([- max(DFF(:)) ,  max(DFF(:))]);
%% Use spike widths and prominence and scatter plot

% take prominence/widths and linearize the vector
lin_prom = prom_spikes';
lin_prom = lin_prom(:);

lin_width = widths_spikes';
lin_width = lin_width(:);

%remove NaNs from the vector
lin_prom = lin_prom(~isnan(lin_prom));
lin_width = lin_width(~isnan(lin_width));

%create cell labels which follow the linearized versions, to be used as
%colormap in plot
for n = 1:size(DFF,1)
    labss(n,:) = prom_spikes(n,:) * 0 + n;
end 
labss_d = labss';
labss_d = labss_d(:);
labss_d = labss_d(~isnan(labss_d));
%% FIGURES
%plot spikeprominence vs width
figure(4567);
scatter (lin_prom,lin_width,[],labss_d,'filled'); 
colormap(turbo); 
a  = colorbar;
a.Label.String = 'Cell Index';

xlabel("Prominence (Intensity)");
ylabel("Width (half prominence)");
%%      NOT NEEDED total number of cells/spikes vs prominence or 1/duration
% 
% subplot(1,2,1); plot(total_number_cells_p); xlabel("Prominence"); ylabel("Total Number of Cells ")
% subplot(1,2,2); plot(total_number_cells_w); xlabel("Width"); ylabel("Total Number of Cells")
%% total number of cells/spikes vs prominence or 1/duration

for val = 1:max(lin_prom(:))
    hol = lin_prom>=val;
    total_number_spikes_p(1,val) = sum(hol(:) == 1);
end 
    nrm_tns_p = total_number_spikes_p/max(total_number_spikes_p);

n = 99;
tot_n = n + 1;
refs = 0: max(1./lin_width(:))/99 :  max(1./lin_width(:));

for val = 1:tot_n
    hol = (1 ./ lin_width) >= refs(val);
    total_number_spikes_w(1,val) = sum(hol(:) == 1);
end 
    nrm_tns_w = total_number_spikes_w/max(total_number_spikes_w);
%% FIGURES
figure(2546);
subplot(2,2,1); plot(nrm_tns_p); xlabel("Prominence"); ylabel("% Total Number of Spikes > Prominence "); hold on; 
subplot(2,2,2); plot(refs,nrm_tns_w); xlabel("1/Width(frame duration)"); ylabel("% Total Number of Spikes > 1/Width(frame duration)"); hold on; 
    
subplot(2,2,3); plot(total_number_spikes_p); xlabel("Prominence"); ylabel("Total Number of Spikes > Prominence"); hold on; 
subplot(2,2,4); plot(refs,total_number_spikes_w); xlabel("Width"); ylabel("Total Number of Spikes"); hold on; 

%%
%save variables as needed to

