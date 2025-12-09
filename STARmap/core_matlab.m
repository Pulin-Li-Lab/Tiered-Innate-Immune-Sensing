function out=core_matlab(tile,sqrt_pieces,mode,t,input_dim,run_id)
    %tile='100_Position350';
    %sqrt_pieces = 2;
    %mode = 'stitching'; % {'global_registration','local_registration','stitching'}
    %t = 1;
    % parameters 
    % input_dim = [100 100 35 4 9]; % [coords_y, coords_x, coords_z, channel, round] % [row, col, z, channel, round]
    ref_round = 1;
    % run_id = 'Brain_batch_v2';
    disk = 'cluster'; % {'cluster', 'local'}
    diskpath = '/stanley/WangLab';
    %arr = split(run_id, '_');
    tissue = run_id;

    if strcmpi(disk,'local')
         diskpath = '/Volumes/stanley_WangLab';
    end
    
    input_path = fullfile(diskpath, 'connie', tissue, '01.data');
    output_path = fullfile(diskpath, 'connie', tissue, '02.processed_data');
    
    addpath(fullfile(diskpath, 'tangzefang/WangLab/starmap/starmap-matlab/src'));
    addpath(fullfile(diskpath, 'tangzefang/WangLab/starmap/starmap-matlab/archive'));
    % 
    useGPU = false;
    % 
    curr_data_dir = tile;
    curr_out_path = fullfile(output_path, tile);
    % 
    if ~exist(curr_out_path, 'dir')
         mkdir(curr_out_path);
         fileattrib(curr_out_path, '+w', 'g'); % allows write permissions for group of users
    end
    
    if strcmp(mode,'global_registration')
    
        sdata = new_STARMapDataset_zf(input_path,output_path, 'useGPU', useGPU);
        sdata.log = fopen(fullfile(curr_out_path, 'log.txt'), 'w');
        sdata = sdata.LoadRawImages('sub_dir', curr_data_dir, 'input_dim', input_dim);
        % sdata = sdata.LoadRawImages('input_dim', input_dim);
        sdata = sdata.SwapChannels; % !!
        sdata = sdata.MinMaxNormalize;
        sdata = sdata.HistEqualize('Method', "inter_round");
        sdata = sdata.HistEqualize('Method', "intra_round");
        % sdata = sdata.MorphoRecon('Method', "2d", 'radius', 5);
        sdata = sdata.MorphoRecon('Method', "2d", 'radius', 6);
    
        %%% registration
        sdata = sdata.test_GlobalRegistration('useGPU', useGPU, 'ref_round', ref_round); % easier to incorporate GPU
    
    
        % Save round 1 image
        output_dir = fullfile(output_path, tile, 'round1_merged');
        if ~exist(output_dir, 'dir')
           mkdir(output_dir);
        end
        r1_img = max(sdata.registeredImages(:,:,:,:,ref_round), [], 4);
        r1_img_name = fullfile(output_dir, sprintf("%s.tif", curr_data_dir));
        SaveSingleTiff(r1_img, r1_img_name);
    
        fclose(sdata.log);
    
    
        coords_mat = table([],[],[],[],[],[],[],[],[],[],[],'VariableNames',{'t','ind_x','ind_y','scoords_x','scoords_y','ecoords_x','ecoords_y','upperleft_x','upperleft_y','inputdim_x','inputdim_y'});
        sub_order = [];
        for i = 0:(sqrt_pieces-1)
            for j = 0:(sqrt_pieces-1)
                sub_order = [sub_order;[i,j]];
            end
        end
        tile_size = floor(input_dim(1) / sqrt_pieces);
        overlap_half = floor(tile_size * 0.1);
    
        upper_left = [0,0];
        for t=1:size(sub_order,1)
            tile_idx = sub_order(t,:);
            start_coords_x = tile_idx(1) * tile_size - overlap_half + 1;
            end_coords_x = (tile_idx(1)+1) * tile_size + overlap_half;
            start_coords_y = tile_idx(2) * tile_size - overlap_half + 1;
            end_coords_y = (tile_idx(2)+1) * tile_size + overlap_half;
            %% compensate in edge
            if tile_idx(1) == 0
                start_coords_x = start_coords_x + overlap_half;
            end
            if tile_idx(2) == 0
                start_coords_y = start_coords_y + overlap_half;
            end
            %% compensate in edge
            if tile_idx(1) == sqrt_pieces - 1
                end_coords_x = input_dim(1);
            end
            if tile_idx(2) == sqrt_pieces - 1
                end_coords_y = input_dim(2);
            end
            upper_left(1) = tile_idx(1) * tile_size;
            upper_left(2) = tile_idx(2) * tile_size;    
    
    
            input_dim_t = input_dim;
            input_dim_t(1:2) = [end_coords_x - start_coords_x + 1,end_coords_y - start_coords_y + 1];
            disp([tile_idx,start_coords_x,end_coords_x,start_coords_y,end_coords_y,upper_left(1:2),input_dim_t(1:2)]);
            coords_mat_t = table(t,tile_idx(1),tile_idx(2),start_coords_x,start_coords_y,end_coords_x,end_coords_y,upper_left(1),upper_left(2),input_dim_t(1),input_dim_t(2),'VariableNames',{'t','ind_x','ind_y','scoords_x','scoords_y','ecoords_x','ecoords_y','upperleft_x','upperleft_y','inputdim_x','inputdim_y'});
            coords_mat = [coords_mat;coords_mat_t];
            t_output = sdata.registeredImages(start_coords_y:end_coords_y,start_coords_x:end_coords_x,:,:,:); %% row - y , col - x [row, col, z, :,:]
    
        %     t_output = sdata.registeredImages;
            save(fullfile(curr_out_path, strcat('registeredImages_','t',num2str(t),'.mat')), 't_output');
        end
        writetable(coords_mat,fullfile(curr_out_path, strcat('coords_mat.csv')),'Delimiter',',','QuoteStrings',false);
    end
    
    
    
    if strcmp(mode,'local_registration')
        coords_mat =readtable(fullfile(curr_out_path,'coords_mat.csv'),'ReadVariableNames',true,'TextType','string');
        goodSpots = table([],[],[],[],'VariableNames',{'x','y','z','Gene'});
    
        input_dim_t = input_dim;
        tile_idx = table2array(coords_mat(t,2:3));
        start_coords_x = table2array(coords_mat(t,4));
        start_coords_y = table2array(coords_mat(t,5));
        upper_left = table2array(coords_mat(t,8:9));
        input_dim_t(1:2) = table2array(coords_mat(t,10:11));
        % [tile_idx(1),tile_idx(2),start_coords_x,start_coords_y,end_coords_x,end_coords_y,upper_left(1),upper_left(2),input_dim_t(1),input_dim_t(2)] = table2array(coords_mat(t,2:11));
    
        sdata_t = new_STARMapDataset_zf(input_path,output_path, 'useGPU', useGPU);
        sdata_t.log = fopen(fullfile(curr_out_path, strcat('log_t',num2str(t),'.txt')), 'w');
        fprintf(sdata_t.log, strcat('log_t',num2str(t),':\n'));
        load(fullfile(curr_out_path, strcat('registeredImages_','t',num2str(t),'.mat')));
        sdata_t.registeredImages = t_output;
        t_output = [];
%         sdata_t.registeredImages = sdata.registeredImages(start_coords_y:end_coords_y,start_coords_x:end_coords_x,:,:,:); %% row - y , col - x [row, col, z, :,:]
        sdata_t.dims = input_dim_t;
        sdata_t.dimX = input_dim_t(1);
        sdata_t.dimY = input_dim_t(2);
        sdata_t.dimZ = input_dim_t(3);
        sdata_t.Nchannel = input_dim_t(4);
        sdata_t.Nround = input_dim_t(5);
    
        sdata_t = sdata_t.xxx_LocalRegistration('Iterations', 50, 'AccumulatedFieldSmoothing', 1, 'ref_round',ref_round);
    
        %sdata_t = sdata_t.LoadCodebook('remove_index', 5);
        sdata_t = sdata_t.LoadCodebook()
        sdata_t = sdata_t.SpotFinding('Method', "max3d", 'ref_index', ref_round, 'showPlots', false);
        sdata_t = sdata_t.ReadsExtraction('voxelSize', [1 1 1]);
        %sdata_t = sdata_t.ReadsFiltration('mode', "duo", 'endBases', ['C', 'T'], 'split_loc', 5, 'showPlots', false);
        sdata_t = sdata_t.ReadsFiltration('showPlots', false);
    
        if size(sdata_t.goodSpots,1) > 0
            sdata_t.goodSpots(:,1) = sdata_t.goodSpots(:,1) + start_coords_x - 1;
            sdata_t.goodSpots(:,2) = sdata_t.goodSpots(:,2) + start_coords_y - 1;
            goodSpots_t = [table(sdata_t.goodSpots(:,1),sdata_t.goodSpots(:,2),sdata_t.goodSpots(:,3),'VariableNames',{'x','y','z'}),cell2table(cellfun(@(x) sdata_t.seqToGene(x), sdata_t.goodReads, 'UniformOutput', false),'VariableNames',{'Gene'})];
        %     kept_idx_t = (sdata_t.goodSpots(:,1) > upper_left(1)) | (sdata_t.goodSpots(:,2) > upper_left(2));
        %     goodSpots_t = goodSpots_t(kept_idx_t,:);
        else
            goodSpots_t = table([],[],[],[],'VariableNames',{'x','y','z','Gene'});
        end
        writetable(goodSpots_t,fullfile(curr_out_path, strcat('goodPoints_max3d_t',num2str(t),'.csv')),'Delimiter',',','QuoteStrings',false);
    
        fclose(sdata_t.log);
    end
    
    
    
    if strcmp(mode,'stitching')
        coords_mat =readtable(fullfile(curr_out_path,'coords_mat.csv'),'ReadVariableNames',true,'TextType','string');
        goodSpots = table([],[],[],[],'VariableNames',{'x','y','z','Gene'});
        for t=1:size(coords_mat,1)
            start_coords_x = table2array(coords_mat(t,4));
            start_coords_y = table2array(coords_mat(t,5));
            upper_left = table2array(coords_mat(t,8:9));
            goodSpots_t = readtable(fullfile(curr_out_path,strcat('goodPoints_max3d_t',num2str(t),'.csv')),'ReadVariableNames',true,'TextType','string');
            %% filter + stitching
            % filter
            if size(goodSpots,1) > 0
                goodSpots = goodSpots((table2array(goodSpots(:,1)) <= upper_left(1)) | (table2array(goodSpots(:,2)) <= upper_left(2)),:);
            end
            if size(goodSpots_t,1) > 0
                kept_idx_t = (table2array(goodSpots_t(:,1)) > upper_left(1)) | (table2array(goodSpots_t(:,2)) > upper_left(2));
                goodSpots_t = goodSpots_t(kept_idx_t,:);
            else
                goodSpots_t = table([],[],[],[],'VariableNames',{'x','y','z','Gene'});
            end
            goodSpots = [goodSpots;goodSpots_t];
        end
    
        writetable(goodSpots,fullfile(curr_out_path, strcat('goodPoints_max3d.csv')),'Delimiter',',','QuoteStrings',false);
        disp(strcat(input_path,"  Work Finished!!!!!"))
    end
