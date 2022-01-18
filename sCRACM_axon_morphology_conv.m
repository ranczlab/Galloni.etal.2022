clear
clc

%% Morphology
cd data     
all_dendrites = csvread('Morphology comparison - Glt V2 wedge analysis.csv',3,2);
oblique_dendrites = csvread('Morphology comparison - Glt V2 oblique wedge analysis.csv',3,2);
cd ..

morphology_scale = -300:10:700;
morphology_scale_expanded = -2000:10:2000;
morphology_scale_interpolated = -1000:1000/48:1000;

distance = 0:10:1000;
smoothing_factor = 5;

n_morphologies = size(all_dendrites,1);
n_distance = size(all_dendrites,2);

%Align tufts by moving 0s from the end of each row to the beginning
flipped_data = fliplr(all_dendrites);
flipped_obliques = fliplr(oblique_dendrites);
flipped_scale = fliplr(morphology_scale);
for cell = 1:n_morphologies
    tuft_end = find(flipped_data(cell,:)~=0,1);
    soma_depth(cell) = flipped_scale(tuft_end);
    shift = tuft_end-2;
    all_dendrites(cell,:) = circshift(all_dendrites(cell,:),shift); %Align morphology plots by pia
    oblique_dendrites(cell,:) = circshift(oblique_dendrites(cell,:),shift); %Align oblique plots by pia
end

depth = mean(soma_depth)
depth_sem = std(soma_depth)/sqrt(n_morphologies)

dendrites_mean = mean(all_dendrites);
dendrites_error = std(all_dendrites)/sqrt(n_morphologies);
obliques_mean = mean(oblique_dendrites);
obliques_error = std(oblique_dendrites)/sqrt(n_morphologies);
basal_tuft = all_dendrites - oblique_dendrites;
basal_tuft_mean = mean(basal_tuft);
basal_tuft_error = std(basal_tuft)/sqrt(n_morphologies);
tuft_mean = mean(basal_tuft);
tuft_error = std(basal_tuft)/sqrt(n_morphologies);
basal_mean = mean(basal_tuft);
basal_error = std(basal_tuft)/sqrt(n_morphologies);
basal_mean(74:end) = 0;
basal_error(74:end) = 0;
tuft_mean(1:73) = 0;
tuft_error(1:73) = 0;

dendrites_mean = smoothdata(dendrites_mean,'gaussian',smoothing_factor);
dendrites_error = smoothdata(dendrites_error,'gaussian',smoothing_factor);
dendrites_mean = fliplr(dendrites_mean);
dendrites_error = fliplr(dendrites_error);

% shadedErrorBar(distance,dendrites_mean,dendrites_error,'lineprops','k',...
%    'transparent',1,'patchSaturation',0.3,'LineWidth',3);
% 
% dendrites_mean = smoothdata(dendrites_mean,'gaussian',smoothing_factor);
% dendrites_error = smoothdata(dendrites_error,'gaussian',smoothing_factor);
% obliques_mean = smoothdata(obliques_mean,'gaussian',smoothing_factor);
% obliques_error = smoothdata(obliques_error,'gaussian',smoothing_factor);
% basal_mean = smoothdata(basal_mean,'gaussian',smoothing_factor);
% basal_error = smoothdata(basal_error,'gaussian',smoothing_factor);
% tuft_mean = smoothdata(tuft_mean,'gaussian',smoothing_factor);
% tuft_error = smoothdata(tuft_error,'gaussian',smoothing_factor);               
% shadedErrorBar(distance,dendrites_mean,dendrites_error,'lineprops','k',...
%     'transparent',1,'patchSaturation',0.3,'LineWidth',3);
% shadedErrorBar(distance,obliques_mean,obliques_error,'lineprops','r',...
%     'transparent',1,'patchSaturation',0.3,'LineWidth',3);
% shadedErrorBar(distance,basal_mean,basal_error,'lineprops','b',...
%     'transparent',1,'patchSaturation',0.3,'LineWidth',3);
% shadedErrorBar(distance,tuft_mean,tuft_error,'lineprops','g',...
%     'transparent',1,'patchSaturation',0.3,'LineWidth',3);

set(gca,'TickDir','out');
box off
xlabel('Distance (µm)')
ylabel('Dendrite length')
title('Glt 3D wedge analysis')
set(gca,'FontSize',13)


%% Axons
cd data
dataORB = csvread('sCRACM_axons_ORB.csv',1,3);
dataATN = csvread('sCRACM_axons_ATN.csv',1,3);
dataRSP = csvread('sCRACM_axons_RSP.csv',1,3);
dataVISp = csvread('sCRACM_axons_VISp.csv',1,3);
dataACA = csvread('sCRACM_axons_ACA.csv',1,3);
dataLP = csvread('sCRACM_axons_LP.csv',1,3);
cd ..

all_data = {dataVISp,dataRSP,dataACA,dataORB,dataATN,dataLP};
plot_number = 1;

for plot_number = 1:6
    colours = {'r','k','m','c','g','y'};
    data = all_data{plot_number};
    smoothing = 149;

    n_rows = size(data,1);
    n_cols = size(data,2);
    pixel_values_all = zeros(n_rows/2,101);

    for i = 1:2:n_rows %for each slice
        dist = data(i,:);
        pixel_values = data(i+1,:);

        [m,idx] = max(dist);
        dist = dist(1:idx); %remove 0s at the end
        pixel_values = pixel_values(1:idx); %remove 0s at the end

        % Interpolate to new distance values
        new_dist = 0:10:m;
        pixel_values = smoothdata(pixel_values,'gaussian',smoothing);
        pixel_values = interp1(dist,pixel_values,new_dist);

        % Baseline and normalize pixel values
        baseline = min(pixel_values);
        pixel_values = pixel_values - baseline;
        pixel_values = pixel_values / max(pixel_values);

        % Store interpolated data
        if length(pixel_values) < 101
            pixel_values_all((i+1)/2,1:length(pixel_values)) = pixel_values;
        else
            pixel_values_all((i+1)/2,:) = pixel_values(1:101);
        end
    end

    mean_values = mean(pixel_values_all);
    sem_axons = std(pixel_values_all)/sqrt(n_rows/2);
    dist = 0:10:1000;

    % Normalize curves
    dendrites_error = dendrites_error/max(dendrites_mean);
    dendrites_mean = dendrites_mean/max(dendrites_mean);
    sem_axons = sem_axons/max(mean_values);
    mean_values = mean_values/max(mean_values);

    % Convolve (multiply) axon profile and morpologies
    convolved_axons = mean_values .* dendrites_mean;
    sem = sem_axons .* dendrites_mean;
    
    sem = sem/max(convolved_axons);
    convolved_axons = convolved_axons/max(convolved_axons);

    a = convolved_axons';
    c = colours{plot_number};
    shadedErrorBar(dist,convolved_axons,sem,'lineprops',colours{plot_number})
%     shadedErrorBar(dist,mean_values,sem_axons,'lineprops',colours{plot_number})
%     shadedErrorBar(distance,dendrites_mean,dendrites_error,'lineprops','k')
end

xlim([0,1000])
%ylim([0,1])
ylabel('fluorescence intensity')
xlabel('distance from pia (µm)')
set(gca,'Fontsize',15)
set(gca,'TickDir','out');
box off    


