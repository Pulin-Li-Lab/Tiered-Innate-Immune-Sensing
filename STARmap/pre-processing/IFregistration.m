
% Register protein images (registration round) to round 1

% parameters 
%withLR = 'yes';
ref_round = 1; %round with DAPI channel
%input_dim = [2048 2048 30 4 1]; %[x, y, z, num_channels, num_rounds], only need to load round 1 since all others are registered to it
run_id = 'IF';
disk = 'cluster'; % {'cluster', 'local'}
diskpath = '/stanley/WangLab';


if strcmpi(disk,'local')
     diskpath = '/Volumes/stanley_WangLab';
end

input_path = fullfile(diskpath, 'connie/08.lung/01.data');
output_path = fullfile(diskpath, 'connie/08.lung/02.processed_data');

addpath(fullfile(diskpath, 'morgan/code/starmap-matlab/src'));
addpath(fullfile(diskpath, 'morgan/code/starmap-matlab/archive'));
% 
useGPU = false;
% 
curr_data_dir = tile;
curr_out_path = fullfile(output_path, tile, run_id)
% 
if ~exist(curr_out_path, 'dir')
     mkdir(curr_out_path)
end

sdata = new_STARMapDataset_zf(input_path,output_path, 'useGPU', useGPU);
sdata.log = fopen(fullfile(curr_out_path, 'registration_log.txt'), 'w');
%sdata = sdata.LoadRawImages('sub_dir', curr_data_dir, 'input_dim', input_dim);
%%sdata = sdata.LoadRawImages('input_dim', input_dim);
%sdata = sdata.SwapChannels; % nuclei registration doesn't use sequencing channels, so no need to swap
%sdata = sdata.MinMaxNormalize;
%sdata = sdata.HistEqualize('Method', "inter_round");
%sdata = sdata.HistEqualize('Method', "intra_round");
%sdata = sdata.MorphoRecon('Method', "2d", 'radius', 3);

%%% registration
sdata = sdata.NucleiRegistrationProtein(run_id, curr_data_dir, 1);
%sdata = sdata.test_GlobalRegistration('useGPU', useGPU, 'ref_round', ref_round); % easier to incorporate GPU
%if strcmpi(withLR, 'yes')
%     sdata = sdata.xxx_LocalRegistration('Iterations', 50, 'AccumulatedFieldSmoothing', 1, 'ref_round',ref_round);
%end

%%% spot finding
%sdata = sdata.LoadCodebook;
%sdata = sdata.SpotFinding('Method', "max3d", 'ref_index', ref_round,  'qualityThreshold', 0.01, ...
%        'volumeThreshold', 30, 'barcodeMethod', "iteration", 'showPlots', false);
%sdata.allSpots = load('centroids.mat').centroids + 1;
%sdata = sdata.ReadsExtraction('showPlots', false, 'voxelSize', [2 2 1]); 
%sdata = sdata.ReadsFiltration('showPlots', false);

% %%% visualization
% bknd_img = max(sdata.registeredImages(:,:,:,:,1), [], 4);
% bknd_img = max(bknd_img, [], 3);
% plot_centroids(sdata.goodSpots, bknd_img, 2, 'blue');

        
% Save round 1 image
%output_dir = fullfile(output_path, 'output', run_id, tile, 'round1_merged');
%if ~exist(output_dir, 'dir')
%   mkdir(output_dir);
%end
%r1_img = max(sdata.registeredImages(:,:,:,:,1), [], 4);
%r1_img_name = fullfile(output_dir, sprintf("%s.tif", curr_data_dir));
%SaveSingleTiff(r1_img, r1_img_name)

% Save cell image
output_dir = fullfile(output_path, run_id);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
sub_dirs = ["dapi", "np", "cd45"];
SaveCellImg(output_dir, sdata.proteinImages, curr_data_dir, sub_dirs);

%%% output results in format x, y, z, gene
%goodSpots = sdata.goodSpots;
%goodReads = sdata.goodReads;
%goodReads = cellfun(@(x) sdata.seqToGene(x), goodReads, 'UniformOutput', false);

%save(fullfile(curr_out_path, strcat('goodPoints_max3d.mat')), 'goodReads', 'goodSpots');
disp(strcat(input_path,'\nRegistration Finished!'))
