clc;
clear;

% Specify the path to AVI file
fly_work_dir = 'G:\BMS Lab\wael chapter3';
fly_folder = 'fly_1';
tifRootPath = fullfile(fly_work_dir, fly_folder);
num_vid = 1;

% get directories of tifs
wingPathVector = cell(num_vid,1);
wholeFlyPathVector = cell(num_vid,1);
for ii = 1:num_vid
    filepath = utils.select_tif_folder(fullfile(tifRootPath,'Wing_Occ_Filled_Uint8'), ii, 'Occ Filled Wing Mask Tifs');
    wingPathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(fullfile(tifRootPath,'Wing_Occ_Filled_Uint8'), ii, 'Occ Filled Wing Mask Tifs');
    wingPathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(fullfile(tifRootPath,'body_mask'), ii, 'Body Mask Tifs');
    wholeFlyPathVector{ii} = string(filepath);
end

% check the number of tif files and size of the first tiff
wing_file_names = cell(num_vid,1);
body_file_names = cell(num_vid,1);
for ii = 1:num_vid
    wing_path = wingPathVector{ii};
    wing_file_names{ii} = dir(fullfile(wing_path, '*.tif'));
    body_path = wholeFlyPathVector{ii};
    body_file_names{ii} = dir(fullfile(body_path, '*.tif'));

    utils.compare_numbers_size(wing_path, body_path, wing_file_names{ii}, body_file_names{ii})
end

% combine wing and body
for ii = 1:num_vid

    % get tif saving path
    sep_path = strsplit(wingPathVector{ii},filesep);
    saveTiffPath = fullfile(tifRootPath, 'Fly_Occ_Filled_Uint8', sep_path{end});
    if ~isfolder(saveTiffPath)
        mkdir(saveTiffPath)
    end
    
    count = 0;
    t_time_elapsed = 0;
    s_time_elapsed = 0;
    for jj = 1:length(body_file_names{ii})
        tic

        wing_tif_dir = fullfile(wingPathVector{ii}, wing_file_names{ii}(jj).name);
        body_tif_dir = fullfile(wholeFlyPathVector{ii}, body_file_names{ii}(jj).name);
        wing_tif = imread(wing_tif_dir);
        body_tif = imread(body_tif_dir);

        % modify body to 40/255 and keep edge
        % define gray scale
        innerGrayValue = 40;
        binaryImage = body_tif; % binarize image
        
        % % % edgeImage = bwperim(binaryImage, 8); % use bwperim to detect edge
        % % % % imshow(edgeImage);
        % % % body_40_image = uint8(body_tif).*40;% set original body to 40/255 gray
        % % % % body_40_image(edgeImage) = 150; % set edge gray scale for body edge
        % % % % imshow(body_40_image);  
        % % % 
        % % % % add body and wing
        % % % whole_fly = imadd(body_40_image, wing_tif);


        % save
        imwrite(whole_fly, fullfile(saveTiffPath, body_file_names{ii}(jj).name));

        % estimate how long need to run
        t_time_elapsed = t_time_elapsed+s_time_elapsed;
        avg_time = t_time_elapsed/jj;
        est_time = (length(body_file_names{ii})-jj)*(avg_time);

        % print progress
        if rem(jj, 100) == 0 || jj == length(body_file_names{ii})
            fprintf(repmat('\b',1,count));
            count = fprintf('folder #: %d, current frame: %s, percentage: %.2f%%, eta: %.2f (sec)\n', ii, body_file_names{ii}(jj).name, (jj)/(length(body_file_names{ii}))*100, est_time);
        end 
        s_time_elapsed = toc;
    end
    % copy vid data from Occ tif
    m_file = dir(fullfile(wingPathVector{ii}, '*.mat'));
    copyfile(fullfile(wingPathVector{ii}, m_file(1).name), fullfile(saveTiffPath,m_file(1).name));
    fprintf('number: %d. wing mask get, and saved to %s.\n', ii, saveTiffPath);
end
disp('All done!')