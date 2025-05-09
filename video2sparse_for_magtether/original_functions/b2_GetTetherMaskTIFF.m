clc;
clear;
% get binary tether mask

% Specify the path to your AVI file
fly_work_dir = 'G:\BMS Lab\wael chapter3';
fly_folder = 'fly_1';
aviRootPath = fullfile(fly_work_dir, fly_folder);
num_vid = 3;
bflag = 3; % flag for bottom view, bot view cam's number, 0 = no bot view cam.
maskBuffer = 15; % the buffer for creating masks
start_frame = 3000;
end_frame = 7000;


%get directories of avis and bgs
nameVector = cell(num_vid,1);
pathVector = cell(num_vid,1);
bgPathVector = cell(num_vid,1);
bgNameVector = cell(num_vid,1);
bgSize = cell(num_vid,1);
background = cell(num_vid,1);
for ii = 1:num_vid
    [filename, filepath] = utils.select_avi(aviRootPath, ii);
    nameVector{ii} = string(filename);
    pathVector{ii} = string(filepath);
end
se = strel('squar',4);
for ii = 1:num_vid
    fprintf('Start for video %d.\n', ii);
    [filename, filepath] = utils.select_bg(aviRootPath, ii);
    bgPathVector{ii} = string(filepath);
    bgNameVector{ii} = string(filename);
    %get size of the bgs
    I = imread(fullfile(filepath,filename));
    I = imcomplement(I);
    I = imdilate(I, se);
    bg_size = size(I);
    bgSize{ii} = bg_size;
    background{ii} = I;
end

% Find initial tether points
tether_points = cell(num_vid,1);
for ii = 1:num_vid
    cur_vid_path = fullfile(pathVector{ii},nameVector{ii});
    tether_points{ii} = utils.set_tether_points(cur_vid_path);
end

% get tether masks (tif)
showplot = false;
for ii = 1:num_vid
    fprintf('start #: %d\n',ii);
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
    utils.make_tether_mask_image(Vreader, maskBuffer, sI, eI, bottomFlag, tether_points{ii}, ...
        background{ii}, saveTetherMaskPath); %function for tether mask creation

    fprintf('%d. tether masks got %s.\n', ii, saveTetherMaskPath);
end
disp('Complete!')





