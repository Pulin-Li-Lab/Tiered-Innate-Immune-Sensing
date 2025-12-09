%% HCR analysis in lung tissue sections 
clear all; close all;

datapath = pwd;
% get all files in maxZ folder
allfiles = dir('maxZ/airway/');
plaque_values1 = {};
plaque_values2 = {};
plaque_keys = [];
for i = 3:length(allfiles)
    imfilename = strcat(datapath, '/maxZ/airway/', allfiles(i).name);
    NP = imread(imfilename, 3);
    IFN = imread(imfilename, 2);
    ISG = imread(imfilename, 4);

    % segmentation from cellpose
    cellpose_file = strcat(datapath, '/maxZ/cellpose/', 'pos', string(i-3), '_cell_masks.txt');
    mask = dlmread(cellpose_file);
    cell_mask = imdilate(mask, strel('disk', 5));
    cellpose_boundary = strcat(datapath, '/maxZ/cellpose/', 'pos', string(i-3), '_cell_outlines.txt');
    cell_outline = dlmread(cellpose_boundary);

     %% Dot counting (global)
    NP_adj = mat2gray(NP);
    IFN_adj = mat2gray(IFN);
    ISG_adj = mat2gray(ISG);
    % subtract background 
    NP_background = imopen(NP_adj,strel('disk',10));
    NP_update = imsubtract(NP_adj, NP_background);
    IFN_background = imopen(IFN_adj,strel('disk',10));
    IFN_update = imsubtract(IFN_adj, IFN_background);
    ISG_background = imopen(ISG_adj,strel('disk',10));
    ISG_update = imsubtract(ISG_adj, ISG_background);
    % count dots both nuclei and cytoplasm 
    [NPpos, NPval, NPnum] = count_dots(NP_update, 7,1);
    [IFNpos, IFNval, IFNnum] = count_dots(IFN_update, 7,1);
    [ISGpos, ISGval, ISGnum] = count_dots(ISG_update, 7,1);

    %% Dot quantification (single-cell) 
    num_region = unique(cell_mask);
    region_cell = zeros(length(num_region)-1, 4);
    
    for n = 1 : length(region_cell)
        region_ind = cell_mask == num_region(n+1,:);
        region_isg = ISGpos.*region_ind;
        count_isg = length(unique(region_isg))-1; 

        region_ifn = IFNpos.*region_ind;
        count_ifn = length(unique(region_ifn))-1;

        region_np = NPpos.*region_ind;
        count_np = length(unique(region_np))-1;

        region_cell(n,1) = num_region(n+1,:);
        region_cell(n, 2) = count_np;
        region_cell(n, 3) = count_ifn;
        region_cell(n, 4) = count_isg;
    end

    np_quant = regionprops(cell_mask, NP, 'Area', 'MeanIntensity');
    np_tot = [np_quant.MeanIntensity]';

    ifn_quant = regionprops(cell_mask, IFN, 'Area', 'MeanIntensity');
    ifn_tot = [ifn_quant.MeanIntensity]';

    isg_quant = regionprops(cell_mask, ISG, 'Area', 'MeanIntensity');
    isg_tot = [isg_quant.MeanIntensity]';

    allcentroid = regionprops(cell_mask, 'Centroid');
    xy_coord = zeros(length(allcentroid), 2);

    for j = 1:length(allcentroid)
        xy_coord(j,1) = allcentroid(j).Centroid(1);
        xy_coord(j,2) = allcentroid(j).Centroid(2);
    end

    cell_state = cat(2, np_tot, ifn_tot, isg_tot, xy_coord);

    plaque_values1{i-2} = cell_state;
    plaque_values2{i-2} = region_cell;

end
