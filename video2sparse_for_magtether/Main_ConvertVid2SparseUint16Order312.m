clc;
clear;
close all;

%% init
fly_work_dir = 'C:\Users\yqz5970\OneDrive - The Pennsylvania State University\Desktop';
fly_folder = 'videos';
aviRootPath = fullfile(fly_work_dir, fly_folder);
num_vid = 1;
start_frame = 100;
end_frame = 130;
actionStartFrame = 1; % at which frame the action starts (the start_frame is the first frame here)
sparse_order312 = 1; % if to order the sparse in the order of 312

%% get directories of avis and bgs
nameVector = cell(num_vid,1); % avi names
pathVector = cell(num_vid,1); % avi dirs
bgPathVector = cell(num_vid,1); % bg dirs
bgNameVector = cell(num_vid,1); % bg names
for ii = 1:num_vid
    [filename, filepath] = utils.select_avi(aviRootPath, ii);
    nameVector{ii} = string(filename);
    pathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    [filename, filepath] = utils.select_bg(aviRootPath, ii);
    bgPathVector{ii} = string(filepath);
    bgNameVector{ii} = string(filename);
end

%% Find initial tether points
tether_points = cell(num_vid,1);
for ii = 1:num_vid
    cur_vid_path = fullfile(pathVector{ii},nameVector{ii});
    tether_points{ii} = utils.set_tether_points(cur_vid_path);
end

%% convert video to tiff directly in Uint8
comp_tiff_dirs = functions.a1_vid2tiff_Uint8(pathVector, nameVector, bgPathVector, bgNameVector, start_frame, end_frame, actionStartFrame,fly_work_dir, fly_folder);

%% get body and tether mask from comp tiff
bflag = 3; % flag for bottom view, bot view cam's number, 0 = no bot view cam.
[body_mask_paths, tether_mask_paths] = functions.b1_GetBodyTetherMaskTIFF(pathVector, nameVector, bgPathVector, bgNameVector, start_frame, end_frame, tether_points, bflag,fly_work_dir, fly_folder);

%% get wing and whole fly mask
[wing_tiff_dirs,whole_fly_tiff_dirs] = functions.c1_GetWingandWholeMaskTIFF_Uint8(comp_tiff_dirs, tether_mask_paths, body_mask_paths,fly_work_dir, fly_folder);

%% Fill wing occlusion and return the whole fly, black background and original color of the wings and bodys (bot white wings and white bodys.)
[filled_whole_fly_paths] = functions.d1_FillWingOcclusion_Uint8(wing_tiff_dirs, tether_mask_paths, whole_fly_tiff_dirs,fly_work_dir, fly_folder);

%% convert tiff to mat file
[mat_paths, mat_file_names] = functions.e_tiff2mat_uint8(filled_whole_fly_paths,fly_work_dir, fly_folder);

%% convert mat to sparse
if sparse_order312
    % in order of 312, where 3 is the z-axis cam (means to put the z-axis cam at the first position, then the side cameras)
    functions.f_mat2sparse_uint8touint16_order312(mat_paths, mat_file_names,fly_work_dir, fly_folder)
else
    % in the original order of 123
    functions.f_mat2sparse_uint8touint16(matDirVector, matNameVector,fly_work_dir, fly_folder)
end
