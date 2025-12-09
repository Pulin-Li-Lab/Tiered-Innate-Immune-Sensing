%% 6hpi HCR
% c1: DAPI - 405
% c2: NP - 488
% c3: IFN - 546
% c4: IRF3 - 594
% c5: NFkB - 647

close all; clear all;

%%

allpos = 1:196;

allcells = zeros(1,9);

for i = 1:length(allpos)
    NP = zeros(2048, 2048);
    IFN = zeros(2048, 2048);
    IRF3 = zeros(2048, 2048);
    NFkB = zeros(2048, 2048);
    pos = allpos(i);
    datapath = strcat('/Volumes/T9/240114_3_6hpi_IF_HCR/6hpi_NP_IFN_IRF3_NFkB-p3/');
    namepre='240114';
    % get NP
    NPfilename = strcat(namepre,'xy',num2str(pos,'%03d'), 'c2', '_spots.tif');
    NP = imread(fullfile(datapath,NPfilename));
    % get IFN
    IFNfilename = strcat(namepre,'xy',num2str(pos,'%03d'), 'c3','_spots.tif');
    IFN = imread(fullfile(datapath,IFNfilename));
    % get IRF3
    IRFfilename = strcat(namepre,'xy',num2str(pos,'%03d'), 'c4','.tif');
    IRF3 = imread(fullfile(datapath,IRFfilename));
    % get NFkB
    NFkBfilename = strcat(namepre,'xy',num2str(pos,'%03d'),'c5','.tif');
    NFkB = imread(fullfile(datapath,NFkBfilename));

    NF_adj = NFkB - 150; %background value 300-RIG / 350-IRF3
    IRF_adj = IRF3 - 500;

    %% segment cell boundary -- use output from Cellpose
    maskfile = strcat(pwd, '/6hpi_NP_IFN_IRF3_NFkB_p3/','cellpose/cyto/', ...
        namepre, 'xy', num2str(pos, '%02d'), 'cyto_masks.txt');
    cyto_mask = dlmread(maskfile);

    nucleimaskfile = strcat(pwd, '/6hpi_NP_IFN_IRF3_NFkB_p3/','cellpose/nuclei/', ...
        namepre, 'xy', num2str(pos, '%02d'), 'nuclei_masks.txt');
    nuclei_mask = dlmread(nucleimaskfile);

    %% Align nuclei and cytoplasm label 
    neg_nuclei = imcomplement(imbinarize(nuclei_mask));
    cyto = neg_nuclei .* cyto_mask;
    nuclei = cyto_mask - cyto;

    %% Dot counting (global)
    NP_adj = mat2gray(NP);
    IFN_adj = mat2gray(IFN);
%     subtract background 
    NP_background = imopen(NP_adj,strel('disk',10));
    NP_update = imsubtract(NP_adj, NP_background);
    IFN_background = imopen(IFN_adj,strel('disk',10));
    IFN_update = imsubtract(IFN_adj, IFN_background);
    % count dots both nuclei and cytoplasm 
    [NPpos, NPval, NPnum] = count_dots(NP_update, 7,1);
    [IFNpos, IFNval, IFNnum] = count_dots(IFN_update, 7,1);
    [RIGpos, RIGval, RIGnum] = count_dots(RIG_update, 7,1);
    overlay1 = imoverlay(NP_adj, cyto_outline);
    NPoverlay = imoverlay(region, NPpos, [1 .3 .3]);
    overlay2 = imoverlay(IFN_adj, cyto_outline);
    IFNoverlay = imoverlay(region, IFNpos, [1 .3 .3]);
    fig1 = figure(); 
    subplot(1,2,1); imshow(NPoverlay);
    subplot(1,2,2); imshow(IFNoverlay);
    savefig(fig1, strcat(datapath,'xy', num2str(pos,'%02d'), '_dots.fig'))
    close(fig1)
    
    %% Dot quantification (single-cell) 
    num_region = unique(nuclei);
    region_cell = zeros(length(num_region)-1, 9);
    
    for n = 1 : length(region_cell)
        region_ind = cyto_mask == num_region(n+1,:);
        region_ifn = IFNpos.*region_ind;
        count_ifn = length(unique(region_ifn))-1;

        region_np = NPpos.*region_ind;
        count_np = length(unique(region_np))-1;

        region_cell(n,1) = num_region(n+1,:);
        region_cell(n, 2) = count_np;
        region_cell(n, 3) = count_ifn;

    end
    nf_val_total = regionprops(cyto, NF_adj, 'MeanIntensity', 'Area');
    nf_avg = [nf_val_total.MeanIntensity];
    cyto_val = [nf_val_total.Area];
    nf_total = nf_avg .* cyto_val;
    nf_val_nuclei = regionprops(nuclei, NF_adj, 'MeanIntensity', 'Area');
    nf_avg_nuclei = [nf_val_nuclei.MeanIntensity];
    nuclei_val = [nf_val_nuclei.Area];
    nf_nuclei = nf_avg_nuclei .* nuclei_val;
    to_keep = find(~isnan(nf_avg_nuclei));
    nf_avg = nf_avg(to_keep);
    region_cell(:,4) = nf_avg_nuclei(~isnan(nf_avg_nuclei));
    region_cell(:,5) = nf_nuclei(~isnan(nf_nuclei))';
    region_cell(:,6) = region_cell(:,4) ./ nf_avg';

    irf_val_total = regionprops(cyto, IRF_adj, 'MeanIntensity', 'Area');
    irf_avg = [irf_val_total.MeanIntensity];
    irf_total = irf_avg .* cyto_val;
    irf_val_nuclei = regionprops(nuclei, IRF_adj, 'MeanIntensity', 'Area');
    irf_avg_nuclei = [irf_val_nuclei.MeanIntensity];
    irf_nuclei = irf_avg_nuclei .* nuclei_val;
    irf_avg = irf_avg(to_keep);

    region_cell(:,7) = irf_avg_nuclei(~isnan(irf_avg_nuclei));
    region_cell(:,8) = irf_nuclei(~isnan(irf_nuclei));
    region_cell(:,9) = region_cell(:,7)./irf_avg';

    allcells = cat(1, allcells, region_cell);

end
 
tosave1 = strcat(datapath, '6hpi_NP_IFN_IRF3_NFkB_p3_airlocalize_thresh100_240125.mat');
save(tosave1, 'allcells');
