clc;
clear;

% Specify the path to your AVI file
fly_work_dir = 'G:\BMS Lab\wael chapter3';
fly_folder = 'fly_1';
tifRootPath = fullfile(fly_work_dir, fly_folder);
num_vid = 3;

% get directories of tifs
pathVector = cell(num_vid,1);
tetherPathVector = cell(num_vid,1);
bodyPathVector = cell(num_vid,1);
for ii = 1:num_vid
    filepath = utils.select_tif_folder(tifRootPath, ii, 'Comp Mask Tifs');
    pathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(tifRootPath, ii, 'Tether Mask Tifs');
    tetherPathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(tifRootPath, ii, 'Body Mask Tifs');
    bodyPathVector{ii} = string(filepath);
end

% check the number of tif files and size of the first tiff
comp_file_names = cell(num_vid,1);
tether_file_names = cell(num_vid,1);
body_file_names = cell(num_vid,1);
for ii = 1:num_vid
    comp_path = pathVector{ii};
    comp_file_names{ii} = dir(fullfile(comp_path, '*.tif'));
    tether_path = tetherPathVector{ii};
    tether_file_names{ii} = dir(fullfile(tether_path, '*.tif'));
    body_path = tetherPathVector{ii};
    body_file_names{ii} = dir(fullfile(body_path, '*.tif'));

    utils.compare_numbers_size(comp_path, tether_path, comp_file_names{ii}, tether_file_names{ii})
    utils.compare_numbers_size(comp_path, body_path, comp_file_names{ii}, body_file_names{ii})
end


% substract tether and body from comp tif
for ii = 1:num_vid

    % get tif saving path
    sep_path = strsplit(pathVector{ii},filesep);
    saveTiffPath = fullfile(tifRootPath, 'wing_mask', sep_path{end});
    if ~isfolder(saveTiffPath)
        mkdir(saveTiffPath)
    end
    
    count = 0;
    t_time_elapsed = 0;
    s_time_elapsed = 0;
    for jj = 1:length(comp_file_names{ii})
        tic
        comp_tif_dir = fullfile(pathVector{ii}, comp_file_names{ii}(jj).name);
        tether_tif_dir = fullfile(tetherPathVector{ii}, tether_file_names{ii}(jj).name);
        body_tif_dir = fullfile(bodyPathVector{ii}, body_file_names{ii}(jj).name);
        comp_tif = imread(comp_tif_dir);
        tether_tif = imread(tether_tif_dir);
        body_tif = imread(body_tif_dir);

        % get removed body and tether tif
        wo_tether_tif = imsubtract(comp_tif, tether_tif);
        wo_tether_tif = imbinarize(wo_tether_tif);
        wo_tether_tif = imfill(wo_tether_tif, 'holes');
        wo_tether_body_tif = imsubtract(wo_tether_tif, body_tif);
        wo_tether_body_tif = imbinarize(wo_tether_body_tif);
        wo_tether_body_tif = imfill(wo_tether_body_tif, 'holes');

        % save
        imwrite(wo_tether_body_tif, fullfile(saveTiffPath, comp_file_names{ii}(jj).name));

        % estimate how long need to run
        t_time_elapsed = t_time_elapsed+s_time_elapsed;
        avg_time = t_time_elapsed/jj;
        est_time = (length(comp_file_names{ii})-jj)*(avg_time);

        % print progress
        if rem(jj, 100) == 0 || jj == length(comp_file_names{ii})
            fprintf(repmat('\b',1,count));
            count = fprintf('folder #: %d, current frame: %s, percentage: %.2f%%, eta: %.2f (sec)\n', ii, comp_file_names{ii}(jj).name, (jj)/(length(comp_file_names{ii}))*100, est_time);
        end 
        s_time_elapsed = toc;
    end
    % copy vid data from Occ tif
    m_file = dir(fullfile(pathVector{ii}, '*.mat'));
    copyfile(fullfile(pathVector{ii}, m_file(1).name), fullfile(saveTiffPath,m_file(1).name));
    fprintf('number: %d. wing mask get, and saved to %s.\n', ii, saveTiffPath);
end
disp('All done!')