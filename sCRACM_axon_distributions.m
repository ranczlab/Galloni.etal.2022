clear
clc

% Load data
cd data
dataORB = csvread('sCRACM_axons_ORB.csv',1,3);
dataATN = csvread('sCRACM_axons_ATN.csv',1,3);
dataRSP = csvread('sCRACM_axons_RSP.csv',1,3);
dataVISp = csvread('sCRACM_axons_VISp.csv',1,3);
dataACA = csvread('sCRACM_axons_ACA.csv',1,3);
dataLP = csvread('sCRACM_axons_LP.csv',1,3);

%Load layer boundaries
layers_ORB = csvread('sCRACM axon boundaries - ORB axons.csv',2,2);
layers_ATN = csvread('sCRACM axon boundaries - ATN axons.csv',2,2);
layers_RSP = csvread('sCRACM axon boundaries - RSP axons.csv',2,2);
layers_VISp = csvread('sCRACM axon boundaries - VISp axons.csv',2,2);
layers_ACA = csvread('sCRACM axon boundaries - ACA axons.csv',2,2);
layers_LP = csvread('sCRACM axon boundaries - LP axons.csv',2,2);
cd ..


all_data = {dataVISp,dataRSP,dataACA,dataORB,dataATN,dataLP};
all_layers = {layers_VISp,layers_RSP,layers_ACA,layers_ORB,layers_ATN,layers_LP};

for plot_number = 1:6
    colours = {'r','k','m','c','g','y'};
    data = all_data{plot_number};
    layers = all_layers{plot_number};
    smoothing = 149;

    n_rows = size(data,1);
    n_cols = size(data,2);
    pixel_values_all = zeros(n_rows/2,n_cols);

    counter = 1;
    for i = 1:2:n_rows %for each slice
        dist = data(i,:);
        pixel_values = data(i+1,:);

        [m,idx] = max(dist);
        dist = dist(1:idx); %remove 0s at the end
        pixel_values = pixel_values(1:idx); %remove 0s at the end

        % Interpolate to new distance values
        new_dist = 0:0.25:m;
        pixel_values = smoothdata(pixel_values,'gaussian',smoothing);
        pixel_values = interp1(dist,pixel_values,new_dist);

        % Baseline and normalize pixel values
        baseline = layers(counter,5);
        pixel_values = max(pixel_values - baseline,0); %set values smaller than baseline to 0
        pixel_values = pixel_values / max(pixel_values);

        % Store interpolated data
        pixel_values_all((i+1)/2,1:size(pixel_values,2)) = pixel_values;
        
        counter = counter+1;
    end

    mean_values = mean(pixel_values_all);
    sem = std(pixel_values_all)/sqrt(n_rows/2);
    sem = sem/max(mean_values);
    mean_values = mean_values/max(mean_values);
    
    dist = 0:0.25:n_cols;
    dist = dist(1:size(mean_values,2));
    
    shadedErrorBar(dist,mean_values,sem,'lineprops',colours{plot_number})
end

xlim([0,1000])
%ylim([0,1])

ylabel('fluorescence intensity (norm.)')
xlabel('distance from pia (µm)')
set(gca,'Fontsize',15)
set(gca,'TickDir','out');
box off
    
    
    
    
    