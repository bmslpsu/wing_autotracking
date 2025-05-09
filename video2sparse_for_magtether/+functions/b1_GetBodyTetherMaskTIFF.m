function [return_body_mask_path, return_tether_mask_path] = b1_GetBodyTetherMaskTIFF(pathVector, nameVector, bgPathVector, bgNameVector, start_frame, end_frame, tether_points, bflag,fly_work_dir, fly_folder)
    % get binary body mask
    
    % Specify the path to your AVI file
    % fly_work_dir = 'G:\BMS Lab\wael chapter3';
    % fly_folder = 'fly_1';
    aviRootPath = fullfile(fly_work_dir, fly_folder);
    num_vid = length(pathVector);
    % bflag = 3; % flag for bottom view, bot view cam's number, 0 = no bot view cam.
    maskBuffer = 15; % the buffer for creating masks
    % start_frame = 3000;
    % end_frame = 7000;
    % has_tether = 1; % 1-has tether;0-no tether
    
    %get directories of avis and bgs
    bgSize = cell(num_vid,1);
    background = cell(num_vid,1);

    se = strel('squar',4);
    for ii = 1:num_vid
        %get size of the bgs
        I = imread(fullfile(bgPathVector{ii},bgNameVector{ii}));
        I = imcomplement(I);
        I = imdilate(I, se);
        bg_size = size(I);
        bgSize{ii} = bg_size;
        background{ii} = I;
    end
    
    
    % get body masks (tif)
    return_body_mask_path = cell(length(nameVector),1);
    for ii = 1:num_vid
        fprintf('start body #: %d\n',ii);
        cur_vid_path = fullfile(pathVector{ii},nameVector{ii});
        Vreader = VideoReader(cur_vid_path);
        if end_frame > Vreader.NumFrames
            error('End_frame exceeds the number of frames. Please check')
        end
        sI = start_frame+maskBuffer;
        eI = end_frame-maskBuffer;
    
        % check if the vid is bot view
        bottomFlag = 0;
        if bflag == ii
            bottomFlag = 1; 
        end
    
        % get body masks saving path
        sep_path = strsplit(pathVector{ii},filesep);
        saveBodyMaskPath = fullfile(aviRootPath, 'body_mask', sep_path{end-1});
        if ~isfolder(saveBodyMaskPath)
            mkdir(saveBodyMaskPath)
        end
        return_body_mask_path{ii} = saveBodyMaskPath;
        utils.make_body_mask_image(Vreader, maskBuffer, sI, eI, bottomFlag, tether_points{ii}, ...
            background{ii}, saveBodyMaskPath); %function for body mask creation
    
        fprintf('%d. body masks got %s.\n', ii, saveBodyMaskPath);
    end
    
    % get tether masks (tif)
    return_tether_mask_path = cell(length(nameVector),1);
    showplot = false;
    for ii = 1:num_vid
        fprintf('start tether #: %d\n',ii);
        cur_vid_path = fullfile(pathVector{ii},nameVector{ii});
        Vreader = VideoReader(cur_vid_path);
        if end_frame > Vreader.NumFrames
            error('End_frame exceeds the number of frames. Please check')
        end
        sI = start_frame+maskBuffer;
        eI = end_frame-maskBuffer;
    
        % check if the vid is bot view
        bottomFlag = 0;
        if bflag == ii
            bottomFlag = 1; 
        end
    
        % get tether masks saving path
        sep_path = strsplit(pathVector{ii},filesep);
        saveTetherMaskPath = fullfile(aviRootPath, 'tether_mask', sep_path{end-1});
        if ~isfolder(saveTetherMaskPath)
            mkdir(saveTetherMaskPath)
        end
        return_tether_mask_path{ii} = saveTetherMaskPath;
        utils.make_tether_mask_image(Vreader, maskBuffer, sI, eI, bottomFlag, tether_points{ii}, ...
            background{ii}, saveTetherMaskPath); %function for tether mask creation
    
        fprintf('%d. tether masks got %s.\n', ii, saveTetherMaskPath);
    end
    disp('Complete! get body and tether masks')

end
