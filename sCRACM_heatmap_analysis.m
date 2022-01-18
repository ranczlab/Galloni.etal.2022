clear
clc

%% Plotting parameters

plot_type = 'heatmap';     % Options: 'heatmap','profile','total_input','morphology'
input = 'ACA';             % Options: 'all','RSP','VISp','V2','ORB','LP','ATN','ACA'
analysed_cell = 'all';     % Options: 'all', or enter cell number to see a particular cell
alignment = 'soma';         % Options: 'soma','pia','all'
cutoff = 7;                % filter signals smaller than cutoff * SD

tuft_only = 0;          
oblique_only = 0;
basal_only = 0;

% Requires plot_type = 'heatmap':
correlation_depth = 0;     % Plot soma depth vs peak location
correlation_total = 0;     % Plot total input vs peak location

stacked_heatmaps = 1;      % Show heatmaps for all cells side-by-side
projection = 'h';          % Options: 'v' = vertical,'h' = horizontal
sort_by = 'peakloc';       % Options: 'depth','total','peakloc','location','mediolateral','allen_slice'

% Requires plot_type = 'profile':
if strcmp(alignment,'pia')
    convolve_morphology = 0;
elseif strcmp(alignment,'soma')
    % Scale input profiles by morphology profile
    convolve_morphology = 1;  % Options: 1=regular profile, 2=split by dendritic domain  
end

% Requires plot_type = 'morphology':
plot_proportions = 0;      % Plot average morphology as proportion of dendrites

%% Load data
cd data;

filename1 = 'sCRACM area - VISp inputs.csv';
filename2 = 'sCRACM area - V2 inputs.csv';
filename3 = 'sCRACM area - RSP inputs.csv';
filename4 = 'sCRACM area - ORB inputs.csv';
filename5 = 'sCRACM area - LP inputs.csv';
filename6 = 'sCRACM area - ATN inputs.csv';
filename7 = 'sCRACM area - ACA inputs.csv';

all_dendrites = csvread('Morphology comparison - Glt V2 wedge analysis.csv',3,2);
oblique_dendrites = csvread('Morphology comparison - Glt V2 oblique wedge analysis.csv',3,2);

filenames = {filename1,filename2,filename3,filename4,filename5,filename6,filename7};
allData = {};
areas={};

for i = 1:length(filenames)
    allData{i} = csvread(filenames{i},3,6); %Load each dataset into a separate cell, ignore first 3 rows/7 cols
    areas{i} = readtable(filenames{i},'VariableNamingRule','preserve');
end

if strcmp(input,'all')
    all_inputs = {'RSP','VISp','V2','ORB','ATN','ACA','LP'};
else
    all_inputs = {input};
end

profile_horizontal = {};
profile_vertical = {};

cd ..;

%% Calculate average morphology (to convolve with heatmaps)
if convolve_morphology ~= 0
    % cd'~/Github_repos/Rancz-Lab/Galloni.etal.2021';
    all_dendrites(find(all_dendrites==0)) = 0.1; %Replace 0s with 0.1 to avoid NaNs when dividing
    n_morphologies = size(all_dendrites,1);
    n_distance = size(all_dendrites,2);

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
    basal_mean(50:end) = 0; %everything above a certain distance belongs to the tuft
    basal_error(50:end) = 0;
    tuft_mean(1:49) = 0; %everything below a certain distance belongs to the soma
    tuft_error(1:49) = 0;
    
    %Calculate soma depth of morphologies
    for cell = 1:n_morphologies
        morphology_soma_depth = length(find(all_dendrites(cell,31:end)>0)) * 10;
        morphology_lengths(cell) = morphology_soma_depth;
    end
    morphology_depth = mean(morphology_lengths);

    %Calculate proportion of each class of dendrites for each distance
    proportion_oblique = obliques_mean ./ dendrites_mean;
    proportion_basal = basal_mean ./ dendrites_mean;
    proportion_tuft = tuft_mean ./ dendrites_mean;
    proportion_total = proportion_oblique + proportion_basal + proportion_tuft;    
    
    %Expand scale of morphology map to fit expanded grid
    proportion_oblique = [zeros(1,170),proportion_oblique,zeros(1,130)]; 
    proportion_basal = [ones(1,170),proportion_basal,zeros(1,130)]; 
    proportion_tuft = [zeros(1,170),proportion_tuft,ones(1,130)]; 
    morphology_scale = -300:10:700;
    morphology_scale_expanded = -2000:10:2000;
    morphology_scale_interpolated = -1000:1000/48:1000;
    
end


%% Generate plots  
for input = all_inputs
    input_number = find(contains(all_inputs,input));
    input = input{1};
    
    input_names = {'VISp','V2','RSP','ORB','LP','ATN','ACA'};
    index = find(strcmp(input_names,{input}));
    data_matrix = allData{index};

    %% Extract relevant data values

    area_names = areas{index}(4:3:end,4).(1);
    allen_slice_number = areas{index}(4:3:end,5).(1);
    mediolateral_distance = areas{index}(4:3:end,6).(1);       
    
    grid_order = data_matrix(1,13:end)';                                       %Grid location of each stim in the order they were presented        
    quadrant_vec = data_matrix(2:3:end,5);                                     %Quadrant of the sCRACM spot the soma is in (1=top-right, 4=bottom-right)
    soma_locations = data_matrix(2:3:end,2:3);                                 %row(1-24) & column(1-12) grid coordinates of soma for each cell 
    soma_row_vec = soma_locations(:,1);
    soma_col_vec = soma_locations(:,2);
    pia_row_vec = data_matrix(2:3:end,6);                                      %row(1-24) grid coordinates of pia at soma column
    cortex_edges_left = data_matrix(2:3:end,7:8);                              %row(1-24) & column(1-12) coordinates of pia at left edge of grid
    cortex_edges_right = data_matrix(2:3:end,10:11);                            %row(1-24) & column(1-12) coordinates of pia at right edge of grid
    pia_delta_row = cortex_edges_right(:,1) - cortex_edges_left(:,1);
    pia_delta_col = cortex_edges_right(:,2) - cortex_edges_left(:,2);
    flip_vector = data_matrix(2:3:end,1);
    
    EPSC_currents = data_matrix(2:3:end,13:300);                                       %Unfiltered EPSC values (min in 20ms window after stim)
    EPSC_areas = data_matrix(3:3:end,13:300);                                          %EPSC areas (30ms window after stim)
    SDs = data_matrix(4:3:end,13:300);                                         %Standard deviation in 40ms window before stim
    n_cells = size(EPSC_currents,1);                                                   %Number of cells
    fullfield_current = data_matrix(2:3:end,302);
    fullfield_area = data_matrix(3:3:end,302);
    
    %Remove EPSC_currents smaller than cutoff*SD
    for cell = 1:n_cells
        for epsc = 1:size(EPSC_currents,2)
            if EPSC_currents(cell,epsc) > (-cutoff * SDs(cell,epsc))
                EPSC_currents(cell,epsc) = 0;
                EPSC_areas(cell,epsc) = 0;
            end

            if EPSC_areas(cell,epsc) > 0   %remove any "negative" areas (i.e. where there is no EPSC)
                EPSC_currents(cell,epsc) = 0;
                EPSC_areas(cell,epsc) = 0;
            end
        end
    end


    %% Generate and average sCRACM heatmaps for all cells
    gridrows = 24;
    gridcols = 12;
    grid_expansion = 2; %Increase the resolution of the grid, to allow for increased precision in soma location

    %Initialise data matrices:
    average_matrix = zeros(4*gridrows+1,4*gridcols+1); %+1 makes grid symmetrical around soma, x4 is for expanded grid (x2) and padding (x2)
    average_matrix_basal = zeros(4*gridrows+1,4*gridcols+1);
    average_matrix_apical = zeros(4*gridrows+1,4*gridcols+1);
    average_morphology_matrix = zeros(4*gridrows+1,4*gridcols+1);
    variance_matrix = zeros(4*gridrows+1,4*gridcols+1); 
    variance_matrix_basal = zeros(4*gridrows+1,4*gridcols+1); 
    variance_matrix_apical = zeros(4*gridrows+1,4*gridcols+1); 
    all_heatmaps = zeros(n_cells,4*gridrows+1,4*gridcols+1);                   
    all_heatmaps_basal = zeros(n_cells,4*gridrows+1,4*gridcols+1);
    all_heatmaps_apical = zeros(n_cells,4*gridrows+1,4*gridcols+1);
    compressed_heatmaps_cols = zeros(4*gridrows+1,n_cells);
    compressed_heatmaps_rows = zeros(n_cells,4*gridcols+1);
    all_morphologies = zeros(n_cells,4*gridrows+1,4*gridcols+1); 

    total_input = zeros(1,n_cells);

    expanded_grid_center = [grid_expansion*gridrows+1,grid_expansion*gridcols+1]; %+1 makes row & col numbers odd to have a symmetric center point
    expanded_pia_loc = [grid_expansion*3*gridcols,grid_expansion*gridcols+1];  %*3 positions pia row at 3/4 of the full grid

    %Define coordinates of soma in expanded grid:
    fx_row = @(x) round(1./(1+1*exp(-10*(x-2.5)))-1);                          %Logistic function to convert quadrant number to row shift, i.e. shift row by 0 or -1 depending on quadrant
    fx_col = @(x) 0.5*(x-2.5).^2 - 9/8;                                        %Parabolic function to convert quadrant number to column shifts, i.e. shift col by 0 or -1 depending on quadrant
    soma_row = grid_expansion*soma_row_vec + fx_row(quadrant_vec);             %without quadrant correction, the soma is by default in quadrant 4                                                                           
    soma_col = grid_expansion*soma_col_vec + fx_col(quadrant_vec);             %the correction applies -1 or 0 to row/col depending on the correct quadrant
            
    %Define pia coordinates in expanded grid 
    pia_row = grid_expansion*pia_row_vec;
    
    %Scale average morphology to soma depth for each cell
    if convolve_morphology == 0
        morphology_scaling_factor = 1;
    else
        morphology_scaling_factor = (pia_row-soma_row)*(1000/48) / morphology_depth;
    end
    
    %Pad data matrix with zeros to center the soma within a larger grid
    if strcmp(alignment,'soma') 
        pad_top = expanded_grid_center(1) - soma_row;
    elseif strcmp(alignment,'pia') 
        pad_top = floor(expanded_pia_loc(1) - pia_row);
        soma_location = pad_top + soma_row;
    end
    pad_left = expanded_grid_center(2) - soma_col;
    pad_right = size(average_matrix,2) - (pad_left+grid_expansion*gridcols);
    pad_bottom = size(average_matrix,1) - (pad_top+grid_expansion*gridrows); 

    if strcmp(analysed_cell,'all')
        for cell_number = 1:n_cells
            data = EPSC_areas(cell_number,:); %data vector contains EPSC area values for sCRACM stimulation of one cell
            data = -data; %turn inward currents into positive numbers

            %Sort data according to stim position
            data_ordered = zeros(1,size(data,2));
            for i=1:length(data)
                position = find(grid_order==i); %index in stim_pos at which number i appears
                data_ordered(i) = data(position);
            end

            %Reshape data from vector to matrix
            heatmap = vec2mat(data_ordered,gridcols);
            if flip_vector(cell_number) == 1
                heatmap = fliplr(heatmap);
            end

            %Normalise data to max current for each cell
            heatmap = heatmap / max(max(heatmap));

            %Align heatmaps:
            %Expand data matrix (by grid_expansion=2) to increase sCRACM resolution
            heatmap_expanded = kron(heatmap,ones(grid_expansion,grid_expansion));

            %Pad data matrix with zeros to center the soma within a larger grid
            heatmap_padded = padarray(heatmap_expanded,[pad_top(cell_number),pad_left(cell_number)],'pre');
            heatmap_padded = padarray(heatmap_padded,[pad_bottom(cell_number),pad_right(cell_number)],'post');

            %% Cutoffs for apical and basal input    
            if strcmp(alignment,'pia')
                basal_only_cutoff = 60;  % 250 µm cutoff
                apical_only_cutoff = 59; % 250 µm cutoff

                % Make sure apical and basal vertical projections are normalized to the same number
                heatmap_padded_basal = heatmap_padded;
                heatmap_padded_basal(basal_only_cutoff:end,:) = 0;
                normalization_basal = max(sum(heatmap_padded_basal));
                heatmap_padded_oblique = heatmap_padded; %added to avoid errors when calling for oblique heatmap further down
                heatmap_padded_oblique(basal_only_cutoff:end,:) = 0; %added to avoid errors when calling for oblique heatmap further down
                normalization_oblique = max(sum(heatmap_padded_oblique)); %added to avoid errors when calling for oblique heatmap further down
                heatmap_padded_tuft = heatmap_padded;
                heatmap_padded_tuft(1:apical_only_cutoff,:) = 0;
                normalization_tuft = max(sum(heatmap_padded_tuft));

                normalization_rows = max([normalization_basal,normalization_tuft]); %normalize to the larger of basal or apical inputs
                   
            elseif strcmp(alignment,'soma')  
                %Expand scale of morphology map to fit soma depth for each cell (in the expanded grid)
                morphology_scale = (-300:10:700)*morphology_scaling_factor(cell_number);
                morphology_scale_expanded = (-2000:10:2000)*morphology_scaling_factor(cell_number);
                morphology_scale_interpolated = -1000:1000/48:1000;
                
                %Sample morphology at same distances as sCRACM maps
                proportion_oblique_scaled = interp1(morphology_scale_expanded,proportion_oblique,morphology_scale_interpolated);
                proportion_basal_scaled = interp1(morphology_scale_expanded,proportion_basal,morphology_scale_interpolated);
                proportion_tuft_scaled = interp1(morphology_scale_expanded,proportion_tuft,morphology_scale_interpolated);

                %Extend morphology profiles to 2D
                proportion_oblique_2D = ones(97,49).* proportion_oblique_scaled';
                proportion_basal_2D = ones(97,49).* proportion_basal_scaled';
                proportion_tuft_2D = ones(97,49).* proportion_tuft_scaled';

                %Convolve heatmaps with morphology and normalize to the same number
                heatmap_padded_basal = heatmap_padded .* proportion_basal_2D;
                normalization_basal = max(sum(heatmap_padded_basal));
                heatmap_padded_tuft = heatmap_padded .* proportion_tuft_2D;
                normalization_tuft = max(sum(heatmap_padded_tuft));
                heatmap_padded_oblique = heatmap_padded .* proportion_oblique_2D;
                normalization_oblique = max(sum(heatmap_padded_oblique));

                normalization_rows = max([normalization_basal,normalization_tuft,normalization_oblique]); %normalize to the larger of basal,oblique, or tuft inputs
            end
            
            if basal_only == 1
                heatmap_padded = heatmap_padded_basal;
            elseif tuft_only == 1
                heatmap_padded = heatmap_padded_tuft;
            elseif oblique_only == 1 
                heatmap_padded = heatmap_padded_oblique;
            elseif (basal_only + tuft_only + oblique_only) == 0
                normalization_rows = max(sum(heatmap_padded));
            end
            
            %Save heatmap and generate projected/compressed heatmaps
            all_heatmaps(cell_number,:,:) = heatmap_padded;
            all_heatmaps_basal(cell_number,:,:) = heatmap_padded_basal;
            all_heatmaps_tuft(cell_number,:,:) = heatmap_padded_tuft;
            all_heatmaps_oblique(cell_number,:,:) = heatmap_padded_oblique;
            
            total_input(cell_number) = -fullfield_area(cell_number);
            
            %Compress heatmaps along X or Y dimensions and stack them, to show individual cells
            compressed_heatmap_row = sum(heatmap_padded) / normalization_rows; %compress and max normalise
            compressed_heatmap_column = sum(heatmap_padded') / max(sum(heatmap_padded'));
            compressed_heatmaps_cols(:,cell_number) = compressed_heatmap_column;     
            compressed_heatmaps_rows(cell_number,:) = compressed_heatmap_row;  

            compressed_heatmap_row_basal = sum(heatmap_padded_basal) / normalization_rows; %compress and max normalise
            compressed_heatmap_column_basal = sum(heatmap_padded_basal') / max(sum(heatmap_padded_basal'));
            compressed_heatmaps_cols_basal(:,cell_number) = compressed_heatmap_column_basal;     
            compressed_heatmaps_rows_basal(cell_number,:) = compressed_heatmap_row_basal;  
            
            compressed_heatmap_row_oblique = sum(heatmap_padded_oblique) / normalization_rows; %compress and max normalise
            compressed_heatmap_column_oblique = sum(heatmap_padded_oblique') / max(sum(heatmap_padded_oblique'));
            compressed_heatmaps_cols_oblique(:,cell_number) = compressed_heatmap_column_oblique;     
            compressed_heatmaps_rows_oblique(cell_number,:) = compressed_heatmap_row_oblique;  
            
            compressed_heatmap_row_tuft = sum(heatmap_padded_tuft) / normalization_rows; %compress and max normalise
            compressed_heatmap_column_tuft = sum(heatmap_padded_tuft') / max(sum(heatmap_padded_tuft'));
            compressed_heatmaps_cols_tuft(:,cell_number) = compressed_heatmap_column_tuft;     
            compressed_heatmaps_rows_tuft(cell_number,:) = compressed_heatmap_row_tuft;     
        end

        % Generate average heatmaps
        averages = mean(all_heatmaps);
        averages_basal = mean(all_heatmaps_basal);
        averages_tuft = mean(all_heatmaps_tuft);
        averages_oblique = mean(all_heatmaps_oblique);

        variances = var(all_heatmaps);
        variances_basal = var(all_heatmaps_basal);
        variances_tuft = var(all_heatmaps_tuft);
        variances_oblique = var(all_heatmaps_oblique);

        variance_matrix(:,:) = variances(1,:,:);                               %map 3D array to 2D matrix
        variance_matrix_basal(:,:) = variances_basal(1,:,:);
        variance_matrix_tuft(:,:) = variances_tuft(1,:,:);
        variance_matrix_oblique(:,:) = variances_oblique(1,:,:);

        average_matrix(:,:) = averages(1,:,:);                                 %map 3D array to 2D matrix
        average_matrix_basal(:,:) = averages_basal(1,:,:);
        average_matrix_tuft(:,:) = averages_tuft(1,:,:);
        average_matrix_oblique(:,:) = averages_oblique(1,:,:);

        %% Display average heatmap
        
        display_matrix = average_matrix;  
        
        %Crop expanded matrix to actual grid dimensions
        crop_bottom = 31;
        if strcmp(alignment,'pia')
            crop_bottom = 24;
        end
        window_height = 48;
        window_width = 24;
        crop_top = size(heatmap_padded,1) - window_height - crop_bottom;
        crop_left = 13;
        crop_right = size(heatmap_padded,2) - window_width - crop_left; 

        display_matrix = display_matrix(crop_bottom:97-crop_top,crop_left:49-crop_right);   
        avg_matrix_basal = average_matrix_basal(crop_bottom:97-crop_top,crop_left:49-crop_right);
        avg_matrix_tuft = average_matrix_tuft(crop_bottom:97-crop_top,crop_left:49-crop_right);
        avg_matrix_oblique = average_matrix_oblique(crop_bottom:97-crop_top,crop_left:49-crop_right);

        var_matrix = variance_matrix(crop_bottom:97-crop_top,crop_left:49-crop_right);
        var_matrix_basal = variance_matrix_basal(crop_bottom:97-crop_top,crop_left:49-crop_right);
        var_matrix_tuft = variance_matrix_tuft(crop_bottom:97-crop_top,crop_left:49-crop_right);
        var_matrix_oblique = variance_matrix_oblique(crop_bottom:97-crop_top,crop_left:49-crop_right);
        
        stdev_matrix = sqrt(var_matrix);
        stdev_matrix_basal = sqrt(var_matrix_basal);
        stdev_matrix_tuft = sqrt(var_matrix_tuft);
        stdev_matrix_oblique = sqrt(var_matrix_oblique);
        
        sem_matrix = stdev_matrix/sqrt(n_cells);
        sem_matrix_basal = stdev_matrix_basal/sqrt(n_cells);
        sem_matrix_tuft = stdev_matrix_tuft/sqrt(n_cells);
        sem_matrix_oblique = stdev_matrix_oblique/sqrt(n_cells);

        unrotated_soma_depth = soma_row+pad_top;
        unrotated_pia_depth = pia_row+pad_top;
        soma_depth = flip(size(average_matrix,1) - (unrotated_soma_depth)) - crop_top;
        pia_depth = flip(size(average_matrix,1) - (unrotated_pia_depth)) - crop_top;

        %Rotate average heatmap to get pia on top in final image (sCRACM data is acquired "upside down")
        rotated_heatmap = imrotate(display_matrix,180); 
        rotated_heatmap_basal = imrotate(avg_matrix_basal,180);
        rotated_heatmap_tuft = imrotate(avg_matrix_tuft,180);
        rotated_heatmap_oblique = imrotate(avg_matrix_oblique,180);
        
        rotated_var = imrotate(var_matrix,180);
        rotated_var_basal = imrotate(var_matrix_basal,180);
        rotated_var_tuft = imrotate(var_matrix_tuft,180);
        rotated_var_oblique = imrotate(var_matrix_oblique,180);
        
        rotated_stdev = imrotate(stdev_matrix,180);
        rotated_stdev_basal = imrotate(stdev_matrix_basal,180);
        rotated_stdev_tuft = imrotate(stdev_matrix_tuft,180);
        rotated_stdev_oblique = imrotate(stdev_matrix_oblique,180);
        
        rotated_sem = imrotate(sem_matrix,180);
        rotated_sem_basal = imrotate(sem_matrix_basal,180);
        rotated_sem_tuft = imrotate(sem_matrix_tuft,180);
        rotated_sem_oblique = imrotate(sem_matrix_oblique,180);

%%
        if strcmp(plot_type,'heatmap') 
            if stacked_heatmaps == 1
                if strcmp(projection,'h') %h = horizontal
                    display_matrix = compressed_heatmaps_cols(crop_bottom:97-crop_top,:);
                    display_matrix_basal = compressed_heatmaps_cols_basal(crop_bottom:97-crop_top,:);
                    display_matrix_oblique = compressed_heatmaps_cols_oblique(crop_bottom:97-crop_top,:);
                    display_matrix_tuft = compressed_heatmaps_cols_tuft(crop_bottom:97-crop_top,:);
                elseif strcmp(projection,'v') %v  = vertical
                    display_matrix = imrotate(compressed_heatmaps_rows(:,crop_left:49-crop_right),90);
                    display_matrix_basal = imrotate(compressed_heatmaps_rows_basal(:,crop_left:49-crop_right),90);
                    display_matrix_oblique = imrotate(compressed_heatmaps_rows_oblique(:,crop_left:49-crop_right),90);
                    display_matrix_tuft = imrotate(compressed_heatmaps_rows_tuft(:,crop_left:49-crop_right),90);
                end

                rotated_heatmap = imrotate(display_matrix,180);
                rotated_heatmap_basal = imrotate(display_matrix_basal,180);
                rotated_heatmap_oblique = imrotate(display_matrix_oblique,180);
                rotated_heatmap_tuft = imrotate(display_matrix_tuft,180);

                if correlation_depth == 1
                    if strcmp(alignment,'soma')
                        display("!! Change alignment to 'pia' !!")
                        break
                    end
                    %Calculate peak location and soma depth (in µm)
                    peakloc = zeros(n_cells,1);
                    peakloc_basal = zeros(n_cells,1);
                    peakloc_oblique = zeros(n_cells,1);
                    peakloc_tuft = zeros(n_cells,1);
                    for cell = 1:n_cells
                        cell_heatmap = rotated_heatmap(:,cell);
                        [peak, peakloc(cell)] = max(cell_heatmap);
                        cell_heatmap_basal = rotated_heatmap_basal(:,cell);
                        [peak, peakloc_basal(cell)] = max(cell_heatmap_basal);
                        cell_heatmap_oblique = rotated_heatmap_oblique(:,cell);
                        [peak, peakloc_oblique(cell)] = max(cell_heatmap_oblique);
                        cell_heatmap_tuft = rotated_heatmap_tuft(:,cell);
                        [peak, peakloc_tuft(cell)] = max(cell_heatmap_tuft);
                    end
                    
                    %Remove cells with no basal/tuft/oblique input (default peakloc = 1)
                    peakloc_basal(find(peakloc_basal==1)) = [];
                    peakloc_oblique(find(peakloc_oblique==1)) = [];
                    peakloc_tuft(find(peakloc_tuft==1)) = [];
                    soma_depth_basal = soma_depth(find(peakloc_basal ~= 1));
                    soma_depth_oblique = soma_depth(find(peakloc_oblique ~= 1));
                    soma_depth_tuft = soma_depth(find(peakloc_tuft ~= 1));
                    
                    %Convert distance from pixel number to µm
                    peakloc = peakloc*1000/48 - 1000/96; %1000/96 places the point in the middle of the first pixel
                    peakloc_basal = peakloc_basal*1000/48 - 1000/96;
                    peakloc_oblique = peakloc_oblique*1000/48 - 1000/96;
                    peakloc_tuft = peakloc_tuft*1000/48 - 1000/96;
                    soma_depth = soma_depth*1000/48; %the soma is already in the middle of its pixel by default
                    soma_depth_basal = soma_depth_basal*1000/48;
                    soma_depth_oblique = soma_depth_oblique*1000/48;
                    soma_depth_tuft = soma_depth_tuft*1000/48;
                                        
                    mean_depth = mean(soma_depth)
                    sem_depth = std(soma_depth)/sqrt(n_cells)

                    %Scatter soma depth vs peak location
                    all_inputs = {'RSP','VISp','V2','ORB','ATN','ACA','LP'};
                    input_number = find(contains(all_inputs,input));
                    colours = {'k','r','b','c','g','m','y'};

                    hold on
                    if correlation_total == 1
                        peakloc = (soma_depth - peakloc);
                        total_input = total_input/1000;
                        scatter(total_input,peakloc,150,'filled',colours{input_number})
                        mdl = fitlm(total_input,peakloc)
                        x_max = ceil(max(total_input)/10)*10;
                        x_regression = [0,x_max];
                        y_regression = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_regression;
                        plot(x_regression,y_regression,'Color','k','Linestyle','--','LineWidth',2)
                        slope = num2str(mdl.Coefficients.Estimate(2));
                        intercept = num2str(mdl.Coefficients.Estimate(1));
                        pval = num2str(mdl.Coefficients.pValue(2));
                        r2 = num2str(mdl.Rsquared.Ordinary);
                        txt1 = strcat('p =',pval);
                        txt2 = strcat('y =',slope,'x ',intercept);
                        txt3 = strcat('r2 =',r2);
                        text(1,200,txt1,'FontSize',30)
                        text(1,250,txt2,'FontSize',30)
                        text(1,300,txt3,'FontSize',30)
                        ylim([-100,500])
                        %xlim([0,2])
                        min(total_input)
                        max(total_input)
                        ylabel('Input distance from soma (µm)','FontSize',30)
                        xlabel('full-field response (pC)','FontSize',30)
                        
                    else
                        scatter(soma_depth_basal,peakloc_basal,150,'filled','b')
                        mdl = fitlm(soma_depth_basal,peakloc_basal)
                        x_regression = [300,700];
                        y_regression = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_regression;
                        plot(x_regression,y_regression,'Color','b','Linestyle','--','LineWidth',2)
                        slope_basal = num2str(mdl.Coefficients.Estimate(2));
                        intercept_basal = num2str(mdl.Coefficients.Estimate(1));
                        pval_basal = num2str(mdl.Coefficients.pValue(2))
                        r2_basal = num2str(mdl.Rsquared.Ordinary)
                        txt1 = strcat('p =',pval_basal);
                        txt2 = strcat('y =',slope_basal,'x ',intercept_basal);
                        txt3 = strcat('r2 =',r2_basal);
                        text(310,680,txt1,'FontSize',30,'Color','b')
                        text(310,730,txt2,'FontSize',30,'Color','b')
                        text(310,780,txt3,'FontSize',30,'Color','b')    
                        
                        scatter(soma_depth_tuft,peakloc_tuft,150,'filled','g')
                        mdl = fitlm(soma_depth_tuft,peakloc_tuft)
                        x_regression = [300,700];
                        y_regression = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_regression;
                        plot(x_regression,y_regression,'Color','g','Linestyle','--','LineWidth',2)
                        slope_tuft = num2str(mdl.Coefficients.Estimate(2));
                        intercept_tuft = num2str(mdl.Coefficients.Estimate(1));
                        pval_tuft = num2str(mdl.Coefficients.pValue(2))
                        r2_tuft = num2str(mdl.Rsquared.Ordinary)
                        txt1 = strcat('p =',pval_tuft);
                        txt2 = strcat('y =',slope_tuft,'x ',intercept_tuft);
                        txt3 = strcat('r2 =',r2_tuft);
                        text(310,10,txt1,'FontSize',30,'Color','g')
                        text(310,60,txt2,'FontSize',30,'Color','g')
                        text(310,110,txt3,'FontSize',30,'Color','g')  
                        ylim([0,800])
                        set(gca,'Ydir','reverse')
                        ylabel('Input distance from pia (µm)','FontSize',30)
                        xlabel('Soma depth (µm)','FontSize',30)
                    end
                    
                    box off
                    set(gca,'tickdir','out')
                    set(gcf, 'Position', [700, 1000, 550, 500])
                    set(gca,'FontSize',30)
                    break
                    
                elseif strcmp(sort_by,'total')
                    %Sort columns by total input (fullfield stim response):                     
                    [sorted_total,total_index] = sort(total_input);
                    soma_depth = soma_depth(total_index);
                    pia_depth = pia_depth(total_index);

                    rotated_heatmap = rotated_heatmap(:,total_index);

                elseif strcmp(sort_by,'depth')    %Sort columns by soma depth:
                    if strcmp(alignment,'pia') 
                        [sorted_depth,depth_index] = sort(soma_depth);
                        soma_depth = sorted_depth;
                        pia_depth = pia_depth(depth_index);
                    elseif strcmp(alignment,'soma') 
                        [sorted_depth,depth_index] = sort(pia_depth);
                        pia_depth = sorted_depth;
                        soma_depth = soma_depth(depth_index);
                    end
                    rotated_heatmap = rotated_heatmap(:,depth_index);

                elseif strcmp(sort_by,'peakloc')    %Sort columns by distance of the peak from the pia or soma 
                    peak_location = zeros(1,n_cells);
                    for cell = 1:n_cells
                        cell_heatmap = rotated_heatmap(:,cell);
                        [peak, peak_location(cell)] = max(cell_heatmap);
                    end

                    peak_location(find(peak_location==1))=NaN;
                    display('horizontal projection distance:')
                    nanmean(soma_depth' - peak_location)*1000/48

                    display('vertical projection distance:')
                    nanmean(peak_location-13)*1000/48                        

                    [sorted_location,sorted_index] = sort(peak_location);
                    
                    %convert to um:
                    if strcmp(projection,'v')
                        peakloc = (sorted_location-13) * 1000/48
                    elseif strcmp(projection,'h')
                        peakloc = sorted_location * 1000/48
                    end
                    
                elseif strcmp(sort_by,'location')
                    VISam_index = find(strcmp(area_names,'VISam'));
                    VISpm_index = find(strcmp(area_names,'VISpm')); 
                    RSPagl_index = find(strcmp(area_names,'RSPagl'));
  
                    total_index = [VISam_index;VISpm_index;RSPagl_index];
                    other = [1:n_cells]';
                    other(total_index) = [];
                    total_index = [total_index;other];
                    
                    soma_depth = soma_depth(total_index);
                    pia_depth = pia_depth(total_index);
                    rotated_heatmap = rotated_heatmap(:,total_index); 
                    
                    %Add blank columns to separate groups of cells by location
                    blank = 0.8*ones(size(rotated_heatmap,1),1);
                    l1 = length(VISam_index);
                    l2 = length(VISpm_index);
                    l3 = length(RSPagl_index);
                    rotated_heatmap = [rotated_heatmap(:,1:l1),blank,...
                        rotated_heatmap(:,1+l1:l1+l2),blank,...
                        rotated_heatmap(:,1+l1+l2:l1+l2+l3),blank,...
                        rotated_heatmap(:,1+l1+l2+l3:end)];
                    blank_depth = 24;
                    soma_depth = [soma_depth(1:l1);blank_depth;...
                        soma_depth(1+l1:l1+l2);blank_depth;...
                        soma_depth(1+l1+l2:l1+l2+l3);blank_depth;...
                        soma_depth(1+l1+l2+l3:end)];                    
                    pia_depth = [pia_depth;0;0;0];
                    display('VISam, VISpm, RSPagl')
                    
                elseif strcmp(sort_by,'mediolateral') %Sort by measured mediolateral distance 
                    [sorted_total,total_index] = sort(mediolateral_distance);
                    soma_depth = soma_depth(total_index);
                    pia_depth = pia_depth(total_index);
                    rotated_heatmap = rotated_heatmap(:,total_index);
                    
                elseif strcmp(sort_by,'allen_slice') %Sort by approximate Allen slice
                    [sorted_total,total_index] = sort(allen_slice_number);
                    soma_depth = soma_depth(total_index);
                    pia_depth = pia_depth(total_index);
                    rotated_heatmap = rotated_heatmap(:,total_index);  
                end
            end

            if stacked_heatmaps == 1 && strcmp(projection,'v') %v  = vertical
                rotated_heatmap = imrotate(rotated_heatmap,-90);
            end
            imagesc(rotated_heatmap);

            % cd'~/Github_repos/Rancz-Lab/Galloni.etal.2021';
            colormap(inferno()); 

            set(gcf, 'Position', [10, 1000, 1400, 1000])
            daspect([1 1 1])
            colorbar
            set(colorbar,'YTick',[])  

            %Add soma locations:
            if stacked_heatmaps == 1
                hold on
                plot(soma_depth,'LineWidth',6,'Color','b')
                plot(pia_depth+1,'LineWidth',6,'Color','w')
                daspect([1 2 1])

                if strcmp(projection,'v')
                    daspect([2 1 1])
                end
            elseif stacked_heatmaps == 0
                hold on
                x = round((expanded_grid_center(2)*ones(n_cells,1))/2);
                if strcmp(alignment,'pia') 
                    soma_depth = soma_depth + randn(n_cells,1)/6; %jitter by small amount to show overlapping cells
                end
                scatter(x,soma_depth,700,'^','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerFaceAlpha',0.9)
                scatter(x,pia_depth,100,'o','MarkerFaceColor','w','MarkerEdgeColor','k')
            end

            set(gca,'XTick',[], 'YTick', [])
            box off
            axis off


        elseif strcmp(plot_type,'profile')
            % cd'~/Github_repos/Rancz-Lab/Galloni.etal.2021';
            colours = {'k','r','b','c','g','m','y'};
            smoothing_factor = 3;

            if convolve_morphology == 2                     
                %Define profiles
                profile_mean = mean(rotated_heatmap,2);
                profile_mean = fliplr(profile_mean');
                profile_basal = mean(rotated_heatmap_basal,2);
                profile_basal = fliplr(profile_basal');
                profile_tuft = mean(rotated_heatmap_tuft,2);
                profile_tuft = fliplr(profile_tuft');
                profile_oblique = mean(rotated_heatmap_oblique,2);
                profile_oblique = fliplr(profile_oblique');

                all_profiles = profile_basal + profile_tuft + profile_oblique;

                morphology_scale = -300:10:700;
                sCRACM_depth_min = (crop_bottom-50)*1000/48;
                sCRACM_depth_max = (crop_bottom-2)*1000/48; 
                step_size = 1000/48;
                sCRACM_scale = sCRACM_depth_min:step_size:sCRACM_depth_max; %soma-aligned input depth (based on crop window, chosen to be near -300:700µm)

%                 %********
%                 %Comment this out to avoid printing
%                 area_all = trapz(sCRACM_scale,profile_mean)
%                 area_tuft = trapz(sCRACM_scale,profile_tuft)
%                 area_basal = trapz(sCRACM_scale,profile_basal)
%                 area_oblique = trapz(sCRACM_scale,profile_oblique)
%                 total_area = area_tuft + area_basal + area_oblique
%                 area_proportion_tuft = area_tuft/total_area
%                 area_proportion_oblique = area_oblique/total_area
%                 area_proportion_basal = area_basal/total_area
%                 %********        
                
                %Smooth profile
                profile_mean = smoothdata(profile_mean,'gaussian',smoothing_factor);
                profile_tuft = smoothdata(profile_tuft,'gaussian',smoothing_factor);
                profile_basal = smoothdata(profile_basal,'gaussian',smoothing_factor);
                profile_oblique = smoothdata(profile_oblique,'gaussian',smoothing_factor);
                        
                hold on

                plot(sCRACM_scale,profile_mean,'k','linewidth',2)
                plot(sCRACM_scale,profile_tuft,'g','linewidth',2)
                plot(sCRACM_scale,profile_basal,'b','linewidth',2)
                plot(sCRACM_scale,profile_oblique,'r','linewidth',2)
                xlabel('Distance from soma (µm)','FontSize',30);
                xlim([-400,600])
                ax = gca;
                %ax.YAxis.Visible = 'off';
                ax.TickDir = 'out';
                ax.FontSize = 20;
                box off       
                
            else
                subplot(2,2,[1,3])
                     if convolve_morphology == 1
                        %Define profile
                        morphology_scale = -300:10:700;
                        sCRACM_depth_min = (crop_bottom-50)*1000/48;
                        sCRACM_depth_max = (crop_bottom-2)*1000/48; 
                        step_size = 1000/48;
                        sCRACM_scale = sCRACM_depth_min:step_size:sCRACM_depth_max; %soma-aligned input depth (based on crop window, chosen to be near -300:700µm)
                        sCRACM_scale = fliplr(sCRACM_scale);                    

                        profile_mean = mean(rotated_heatmap,2);
                        error_horizontal = mean(rotated_sem,2);
                        profile_tuft = mean(rotated_heatmap_tuft,2);
                        error_tuft = mean(rotated_sem_tuft,2);
                        profile_basal = mean(rotated_heatmap_basal,2);
                        error_basal = mean(rotated_sem_basal,2);
                        profile_oblique = mean(rotated_heatmap_oblique,2);
                        error_oblique = mean(rotated_sem_oblique,2);

                        [apical_peak, apical_index] = max(profile_mean(1:13));
                        apical_location = apical_index*1000/48 - 1000/96
                        [perisomatic_peak, perisomatic_index] = max(profile_mean(14:end));
                        perisomatic_location = (13+perisomatic_index)*1000/48 - 1000/96

                        area_all = trapz(sCRACM_scale,profile_mean);
                        area_tuft = trapz(sCRACM_scale,profile_tuft);
                        area_basal = trapz(sCRACM_scale,profile_basal);
                        area_oblique = trapz(sCRACM_scale,profile_oblique);
                        total_area = area_tuft + area_basal + area_oblique;
                        area_proportion_tuft = area_tuft/total_area;
                        area_proportion_oblique = area_oblique/total_area;
                        area_proportion_basal = area_basal/total_area;

                        [tuft_peak, tuft_index] = max(profile_tuft);
                        tuft_location = sCRACM_scale(tuft_index)
                        [oblique_peak, oblique_index] = max(profile_oblique);
                        oblique_location = sCRACM_scale(oblique_index)
                        [basal_peak, basal_index] = max(profile_basal);
                        basal_location = sCRACM_scale(basal_index)

                        %Smooth profile
                        profile_mean = smoothdata(profile_mean,'gaussian',smoothing_factor);
                        error_horizontal = smoothdata(error_horizontal,'gaussian',smoothing_factor);

                        %Normalize profile
                        error_horizontal = error_horizontal/max(profile_mean);
                        profile_horizontal{input_number} = profile_mean/max(profile_mean);
                        
                        if strcmp(alignment,'pia')
                            sCRACM_scale = 1000/96:step_size:97*1000/96;
                        end

                        hold on
                        shadedErrorBar(sCRACM_scale,profile_horizontal{input_number},error_horizontal,'lineprops',colours{input_number},...
                            'transparent',1,'patchSaturation',0.3,'LineWidth',3);

                        ax = gca;
                        ax.TickDir = 'out';
                        ax.XAxis.FontSize = 15;
                        ax.YAxis.Visible = 'off';
                        if strcmp(alignment,'pia')
                            xlim([0,1000])
                            xlabel('Distance from pia (µm)','FontSize',30);
                            set(gca, 'XDir','reverse')
                        elseif strcmp(alignment,'soma')
                            xlim([-400,600])
                            set(gca,'XTick',[-400:200:600])
                            xticklabels(-400:200:600);
                            xlabel('Distance from soma (µm)','FontSize',30);
                        end
                        view([90 -90]);
                        box off                    

                    else
                        %Define profile
                        depth = [0:1000/48:1000]; %depth from pia
                        profile_mean = mean(rotated_heatmap,2);
                        error_horizontal = mean(rotated_sem,2);

                        [apical_peak, apical_index] = max(profile_mean(1:13));
                        apical_location = (apical_index-1)*1000/48 %-1 because the first point is at distance 0
                        [perisomatic_peak, perisomatic_index] = max(profile_mean(14:end));
                        perisomatic_location = (13+perisomatic_index-1)*1000/48

                        %Smooth profile
                        profile_mean = smoothdata(profile_mean,'gaussian',smoothing_factor);
                        error_horizontal = smoothdata(error_horizontal,'gaussian',smoothing_factor);

                        %Normalize profile
                        error_horizontal = error_horizontal/max(profile_mean);
                        profile_horizontal{input_number} = profile_mean/max(profile_mean);

                        a = profile_horizontal{input_number}
                        shadedErrorBar(depth,profile_horizontal{input_number},error_horizontal,'lineprops',colours{input_number},...
                            'transparent',1,'patchSaturation',0.3,'LineWidth',3);

                        ax = gca;
                        ax.TickDir = 'out';
                        ax.XAxis.FontSize = 20;
                        ax.YAxis.Visible = 'off';
                        set(gca,'XTick',[0:250:1000])
                        if strcmp(alignment,'pia')
                            xticklabels(0:250:1000);
                            xlabel('Distance from pia (µm)','FontSize',30);
                        elseif strcmp(alignment,'soma')
                            %xticklabels(500:-250:-500);
                            xticklabels(700:-250:-300);
                            xlabel('Distance from soma (µm)','FontSize',30);
                        end
                        view([90 90]);
                        box off
                    end


                subplot(2,2,2)
                    %Define profile
                    profile_mean = mean(rotated_heatmap);
                    profile_basal = mean(rotated_heatmap_basal);
                    profile_tuft = mean(rotated_heatmap_tuft);
                    profile_oblique = mean(rotated_heatmap_oblique);
                    error_vertical = mean(rotated_sem);

                    %Smooth profile
                    error_vertical = smoothdata(error_vertical,'gaussian',smoothing_factor);
                    profile_mean = smoothdata(profile_mean,'gaussian',smoothing_factor);
                    profile_tuft = smoothdata(profile_tuft,'gaussian',smoothing_factor);
                    profile_basal = smoothdata(profile_basal,'gaussian',smoothing_factor);
                    profile_oblique = smoothdata(profile_oblique,'gaussian',smoothing_factor);
                    
                    %Measure peak locations
                    [lateral_peak, lateral_index] = max(profile_mean);
                    lateral_location = (lateral_index-13)*1000/48   
                    [lateral_peak_tuft, lateral_index] = max(profile_tuft);
                    lateral_location_tuft = (lateral_index-13)*1000/48   
                    [lateral_peak_basal, lateral_index] = max(profile_basal);
                    lateral_location_basal = (lateral_index-13)*1000/48  
                    [lateral_peak_oblique, lateral_index] = max(profile_oblique);
                    lateral_location_oblique = (lateral_index-13)*1000/48

                    %Normalize profile
                    normalization_basal = max(profile_basal);
                    normalization_tuft = max(profile_tuft);
                    normalization_oblique = max(profile_oblique);
                    normalization_number = max([normalization_basal,normalization_tuft,normalization_oblique]);

                    if (basal_only + tuft_only + oblique_only) == 0
                        normalization_number = max(profile_mean);
                    end

                    error_vertical = error_vertical/normalization_number;
                    profile_vertical{input_number} = profile_mean/normalization_number;

                    width = [-250:500/24:250];
                    shadedErrorBar(width,profile_vertical{input_number},error_vertical,'lineprops',colours{input_number},...
                        'transparent',1,'patchSaturation',0.3,'LineWidth',3);
                    ax = gca;
                    ax.TickDir = 'out';
                    ax.XAxis.FontSize = 20;
                    ax.YAxis.Visible = 'off';
                    %set(gca,'XTick',[0:100:500])
                    %xticklabels(-250:100:250);
                    xlabel('Distance from soma (µm)','FontSize',30);
                    box off

                set(gcf, 'Position', [10, 1000, 1000, 1000])
            end

        elseif strcmp(plot_type,'total_input')
            % cd'~/Github_repos/Rancz-Lab/Galloni.etal.2021';
            total_input = total_input/1000;
            
            colours = {'k','r','b','c','g','m','y'};
            all_inputs = {'RSP','VISp','V2','ORB','ATN','ACA','LP'};
            limits = {[0,12],[0,2],[0,25],[0,25],[0,6],[0,25],[0,2]};
            input_number = find(contains(all_inputs,input));

            hold on
            boxplot(total_input,'Colors','k','Whisker',10)
            scatter(0.5,total_input,100,'MarkerFaceColor',colours{input_number},'MarkerEdgeColor','none')
            ylim(limits{input_number})
              
            mean_current = mean(-fullfield_current)
            mean_area = mean(total_input)
            SEM = std(total_input)/sqrt(length(total_input))

            set(gca,'TickDir','out')
            box off
            title('Total input')
            
        elseif strcmp(plot_type,'morphologies_sholl')
            % cd'~/Github_repos/Rancz-Lab/Galloni.etal.2021'     

            sholl_radius = 0:10:680;
            sholl_data_Glt = sholl_data_Glt(:,26:end);
            n_cells = size(sholl_data_Glt,1);
            n_sholl = size(sholl_data_Glt,2);

            if strcmp(alignment,'pia') %Move 0s from the end of each row to the beginning
                for cell = 1:n_cells
                    shift = n_sholl - find(sholl_data_Glt(cell,:) == 0,1)+1;
                    sholl_data_Glt(cell,:) = circshift(sholl_data_Glt(cell,:),shift);
                end
                sholl_data_Glt = [sholl_data_Glt,zeros(n_cells,1)]; %Append zeros to make plot symmetrical
            elseif strcmp(alignment,'soma')
                sholl_data_Glt = [zeros(n_cells,1),sholl_data_Glt]; %Append zeros to make plot symmetrical
            end

            sholl_mean = mean(sholl_data_Glt);
            sholl_error = std(sholl_data_Glt)/sqrt(n_cells);

            shadedErrorBar(sholl_radius,sholl_mean,sholl_error,'lineprops','k',...
                'transparent',1,'patchSaturation',0.3,'LineWidth',3);

            xlim([0,700])
            ylim([0,35])
            set(gca,'TickDir','out');
            box off
            xlabel('Distance (µm)')
            ylabel('Intersections')
            title('Glt Sholl plot')
            set(gca,'FontSize',13)

        elseif strcmp(plot_type,'morphology')
            % cd'~/Github_repos/Rancz-Lab/Galloni.etal.2021'     
            plot_data = all_dendrites;
            plot_obliques = oblique_dendrites;
            plot_data(find(plot_data==0)) = 0.1; %Replace 0s with 0.1 to avoid NaNs when dividing

            distance = -300:10:700;

            smoothing_factor = 5;

            n_cells = size(plot_data,1);
            n_distance = size(plot_data,2);
            if strcmp(alignment,'pia') %Align tufts by moving 0s from the end of each row to the beginning
                for cell = 1:n_cells
                    flipped_data = fliplr(plot_data);
                    flipped_obliques = fliplr(plot_obliques);
                    shift = find(flipped_data(cell,:)~=0,1)-2;
                    plot_data(cell,:) = circshift(plot_data(cell,:),shift); %Align morphology plots by pia
                    plot_obliques(cell,:) = circshift(plot_obliques(cell,:),shift); %Align oblique plots the same way
                end
            end

            dendrites_mean = mean(plot_data);
            dendrites_error = std(plot_data)/sqrt(n_cells);
            obliques_mean = mean(plot_obliques);
            obliques_error = std(plot_obliques)/sqrt(n_cells);
            basal_tuft = plot_data - plot_obliques;
            basal_tuft_mean = mean(basal_tuft);
            basal_tuft_error = std(basal_tuft)/sqrt(n_cells);
            tuft_mean = mean(basal_tuft);
            tuft_error = std(basal_tuft)/sqrt(n_cells);
            basal_mean = mean(basal_tuft);
            basal_error = std(basal_tuft)/sqrt(n_cells);
            basal_mean(50:end) = 0;
            basal_error(50:end) = 0;
            tuft_mean(1:49) = 0;
            tuft_error(1:49) = 0;

            proportion_oblique = obliques_mean ./ dendrites_mean;
            proportion_basal = basal_mean ./ dendrites_mean;
            proportion_tuft = tuft_mean ./ dendrites_mean;
            proportion_total = proportion_oblique + proportion_basal + proportion_tuft;
            
            area_all = trapz(distance,dendrites_mean);
            area_tuft = trapz(distance,tuft_mean)/area_all
            area_oblique = trapz(distance,obliques_mean)/area_all
            area_basal = trapz(distance,basal_mean)/area_all
            

            if plot_proportions == 1
                hold on
                plot(distance,proportion_oblique,'r','linewidth',2)
                plot(distance,proportion_basal,'b','linewidth',2)
                plot(distance,proportion_tuft,'g','linewidth',2)
                xlabel('Distance (µm)')
                ylabel('Proportion of dendrites')
                break
            end

            dendrites_mean = smoothdata(dendrites_mean,'gaussian',smoothing_factor);
            dendrites_error = smoothdata(dendrites_error,'gaussian',smoothing_factor);
            obliques_mean = smoothdata(obliques_mean,'gaussian',smoothing_factor);
            obliques_error = smoothdata(obliques_error,'gaussian',smoothing_factor);
            basal_mean = smoothdata(basal_mean,'gaussian',smoothing_factor);
            basal_error = smoothdata(basal_error,'gaussian',smoothing_factor);
            tuft_mean = smoothdata(tuft_mean,'gaussian',smoothing_factor);
            tuft_error = smoothdata(tuft_error,'gaussian',smoothing_factor);               

            %Convert from dendritic length to probability density function (area=1)
%                 area = trapz(distance,dendrites_mean);
%                 dendrites_mean = dendrites_mean/area;
%                 dendrites_error = dendrites_error/area;
%                 obliques_mean = obliques_mean/area;
%                 obliques_error = obliques_error/area;
%                 basal_mean = basal_mean/area;
%                 basal_error = basal_error/area;
%                 tuft_mean = tuft_mean/area;
%                 tuft_error = tuft_error/area;

            shadedErrorBar(distance,dendrites_mean,dendrites_error,'lineprops','k',...
                'transparent',1,'patchSaturation',0.3,'LineWidth',3);
            shadedErrorBar(distance,obliques_mean,obliques_error,'lineprops','r',...
                'transparent',1,'patchSaturation',0.3,'LineWidth',3);
            shadedErrorBar(distance,basal_mean,basal_error,'lineprops','b',...
                'transparent',1,'patchSaturation',0.3,'LineWidth',3);
            shadedErrorBar(distance,tuft_mean,tuft_error,'lineprops','g',...
                'transparent',1,'patchSaturation',0.3,'LineWidth',3);

            %xlim([-200,700])
            %ylim([0,320])
            set(gca,'TickDir','out');
            box off
            xlabel('Distance (µm)')
            ylabel('Dendrite length')
            title('Glt 3D wedge analysis')
            set(gca,'FontSize',13)
        end


    elseif strcmp(analysed_cell,'all') == 0
        if analysed_cell == -1
            data = EPSC_areas(end,:);                                           %data vector contains EPSC amplitude values for sCRACM stimulation of one cell
        else
            data = EPSC_areas(analysed_cell,:);                                 
        end

        %Turn inward currents into positive numbers
        data = -data;

        data_ordered = zeros(1,size(data,2));

        %Sort data according to stim position
        for i=1:length(data)
            position = find(grid_order==i);                                %index in stim_pos at which number i appears
            data_ordered(i) = data(position);
        end

        heatmap = vec2mat(data_ordered,gridcols);

        %% Display heatmap
        % cd'~/Github_repos/Rancz-Lab/Galloni.etal.2021';
        rotated_heatmap = imrotate(heatmap,180);
        imagesc(rotated_heatmap)
        %colormap(hot); 
        colormap(inferno());
        set(gcf, 'Position', [10, 1000, 500, 1000])
        axis off
        axis tight %make axes square
        colorbar
    end
end

% cd'~/Github_repos/Rancz-Lab/Galloni.etal.2021';

