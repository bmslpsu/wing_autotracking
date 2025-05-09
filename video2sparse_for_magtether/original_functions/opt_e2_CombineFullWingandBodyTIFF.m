clc;
clear;

% Specify the path to AVI file
fly_work_dir = 'G:\BMS Lab\wael chapter3';
fly_folder = 'fly_1';
tifRootPath = fullfile(fly_work_dir, fly_folder);
num_vid = 3;

% get directories of tifs
wingPathVector = cell(num_vid,1);
bodyPathVector = cell(num_vid,1);
for ii = 1:num_vid
    filepath = utils.select_tif_folder(tifRootPath, ii, 'Occ Filled Wing Mask Tifs');
    wingPathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(tifRootPath, ii, 'Body Mask Tifs');
    bodyPathVector{ii} = string(filepath);
end

% check the number of tif files and size of the first tiff
wing_file_names = cell(num_vid,1);
body_file_names = cell(num_vid,1);
for ii = 1:num_vid
    wing_path = wingPathVector{ii};
    wing_file_names{ii} = dir(fullfile(wing_path, '*.tif'));
    body_path = bodyPathVector{ii};
    body_file_names{ii} = dir(fullfile(body_path, '*.tif'));

    utils.compare_numbers_size(wing_path, body_path, wing_file_names{ii}, body_file_names{ii})
end

% combine wing and body
for ii = 1:num_vid

    % get tif saving path
    sep_path = strsplit(wingPathVector{ii},filesep);
    saveTiffPath = fullfile(tifRootPath, 'Whole_Fly_Uint8', sep_path{end});
    if ~isfolder(saveTiffPath)
        mkdir(saveTiffPath)
    end
    
    count = 0;
    t_time_elapsed = 0;
    s_time_elapsed = 0;
    for jj = 1:length(body_file_names{ii})
        tic

        wing_tif_dir = fullfile(wingPathVector{ii}, wing_file_names{ii}(jj).name);
        body_tif_dir = fullfile(bodyPathVector{ii}, body_file_names{ii}(jj).name);
        wing_tif = imread(wing_tif_dir);
        body_tif = imread(body_tif_dir);

        % combine body and wing output uint8 for gray images
        wing_uint8 = uint8(wing_tif) .* round(255/2);
        body_uint8 = uint8(body_tif) .* 255;
        whole_fly = imadd(wing_uint8, body_uint8);
        % imshow(whole_fly)


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