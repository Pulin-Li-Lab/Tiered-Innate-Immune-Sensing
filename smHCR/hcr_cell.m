%% Cell culture HCR
close all; clear all;

%%

allpos = 1:19;

huvec_wsn = zeros(1,4);
frac_ifn_huvec = zeros(20,1);

for i = 1:length(allpos)
    NP = zeros(2048, 2048, 10);
    IFN = zeros(2048, 2048, 10);
    RIG = zeros(2048, 2048, 10);
    pos = allpos(i);
    datapath = strcat(pwd, '/HUVEC_WSN/');
    namepre='230406';
    for j = 3:13
        % get NP
        NPfilename = strcat(namepre,'xy',num2str(pos,'%02d'),'z',num2str(j,'%02d'), 'c2','.tif');
        NP(:,:,j) = imread(fullfile(datapath,NPfilename));
        % get IFN
        IFNfilename = strcat(namepre,'xy',num2str(pos,'%02d'),'z',num2str(j,'%02d'), 'c4','.tif');
        IFN(:,:,j) = imread(fullfile(datapath,IFNfilename));
        % get RIG-I
        RIGfilename = strcat(namepre,'xy',num2str(pos,'%02d'),'z',num2str(j,'%02d'), 'c3','.tif');
        RIG(:,:,j) = imread(fullfile(datapath,RIGfilename));

    end

    zNP = max(NP, [], 3);
    zIFN = max(IFN, [], 3);
    zRIG = max(RIG, [], 3);


    %% segment cell boundary -- use output from Cellpose
    maskfile = strcat(pwd, '/HUVEC_WSN/','cellpose/', ...
        namepre, 'xy', num2str(pos, '%02d'), 'nuclei_masks.txt');
    cyto_mask = dlmread(maskfile);

    outlinefile = strcat(pwd, '/HUVEC_WSN/','cellpose/', ...
        namepre, 'xy', num2str(pos, '%02d'), 'nuclei_outlines.txt');
    cyto_outline = dlmread(outlinefile);

    cyto = imdilate(cyto_mask, strel('disk',30));
%     cyto_outline = imdilate(cyto_outline, strel('disk',10));
    %% Dot counting (global)
%     zIFN (zIFN >= 10000) = 0;
    NP_adj = mat2gray(zNP);
    IFN_adj = mat2gray(zIFN, [100 5000]);
    RIG_adj = mat2gray(zRIG);
%     subtract background 
    NP_background = imopen(NP_adj,strel('disk',10));
    NP_update = imsubtract(NP_adj, NP_background);
    IFN_background = imopen(IFN_adj,strel('disk',10));
    IFN_update = imsubtract(IFN_adj, IFN_background);
    RIG_background = imopen(RIG_adj,strel('disk',10));
    RIG_update = imsubtract(RIG_adj, RIG_background);
    % count dots both nuclei and cytoplasm 
    [NPpos, NPval, NPnum] = count_dots(NP_update, 7,1);
    [IFNpos, IFNval, IFNnum] = count_dots(IFN_update, 7,1);
    [RIGpos, RIGval, RIGnum] = count_dots(RIG_update, 7,1);
    overlay1 = imoverlay(NP_adj, cyto_outline);
    NPoverlay = imoverlay(overlay1, NPpos, [0 1 1]);
    overlay2 = imoverlay(IFN_adj, cyto_outline);
    IFNoverlay = imoverlay(overlay2, IFNpos, [.3 1 .3]);
    overlay3 = imoverlay(RIG_adj, cyto_outline);
    RIGoverlay = imoverlay(overlay3, RIGpos, [1 .3 .3]);
    fig1 = figure(); 
    subplot(1,3,1); imshow(NPoverlay);
    subplot(1,3,2); imshow(IFNoverlay);
    subplot(1,3,3); imshow(RIGoverlay);
%     savefig(fig1, strcat(datapath,'xy', num2str(pos,'%02d'), '_dots.fig'))
%     close(fig1)
    
    %% Dot quantification (single-cell) 
    num_region = unique(cyto);
    region_cell = zeros(length(num_region)-1, 4);
    
    for n = 1 : length(region_cell)
        region_ind = cyto == num_region(n+1,:);
        region_rig = RIGpos.*region_ind;
        count_rig = length(unique(region_rig))-1; 

        region_ifn = IFNpos.*region_ind;
        count_ifn = length(unique(region_ifn))-1;

        region_np = NPpos.*region_ind;
        count_np = length(unique(region_np))-1;

        region_cell(n,1) = num_region(n+1,:);
        region_cell(n, 2) = count_np;
        region_cell(n, 3) = count_ifn;
        region_cell(n, 4) = count_rig;
    end

end
% 
% tosave1 = strcat(datapath, '6hpi_A2_quant.mat');
% save(tosave1, 'a549_control');


