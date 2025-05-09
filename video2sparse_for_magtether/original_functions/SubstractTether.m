clc;
clear;

% init para
fly_work_dir = 'G:\BMS Lab\wael chapter3';
fly_folder = 'fly_1';
root = fullfile(fly_work_dir, fly_folder);
num_vid = 3;

%get directories of complete tiffs and tether tiffs
compPathVector = cell(num_vid,1);
tetherPathVector = cell(num_vid,1);
for ii = 1:num_vid
    filepath = utils.select_tif_folder(root, ii, 'comp');
    compPathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(root, ii, 'tether');
    tetherPathVector{ii} = string(filepath);
end

% check the number of tif files and size of the first tiff
comp_file_names = cell(num_vid,1);
tether_file_names = cell(num_vid,1);
for ii = 1:num_vid
    comp_path = compPathVector{ii};
    comp_file_names{ii} = dir(fullfile(comp_path, '*.tif'));
    tether_path = tetherPathVector{ii};
    tether_file_names{ii} = dir(fullfile(tether_path, '*.tif'));

    % check number of files
    if length(comp_file_names) ~= length(tether_file_names)
        error('the number of frames are not the same, comp#: %d, tether#: %d', length(comp_file_names), length(tether_file_names))
    end
    %check size of first tiff
    fcomp_size = size(imread(fullfile(comp_path, comp_file_names{ii}(1).name)));
    ftether_size = size(imread(fullfile(tether_path, tether_file_names{ii}(1).name)));
    if (fcomp_size) ~= (ftether_size)
        error('size of the frames are not the same, comp#: [%d, %d], tether#: [%d, %d]', fcomp_size(1), fcomp_size(2), ftether_size(1), ftether_size(2))
    end
end

% substract tether from comp tif
for ii = 1:num_vid

    % get tif saving path
    sep_path = strsplit(compPathVector{ii},filesep);
    saveTiffPath = fullfile(root, 'without_tether', sep_path{end});
    if ~isfolder(saveTiffPath)
        mkdir(saveTiffPath)
    end
    
    count = 0;
    for jj = 1:length(comp_file_names{ii})
        comp_tif_dir = fullfile(compPathVector{ii}, comp_file_names{ii}(jj).name);
        tether_tif_dir = fullfile(tetherPathVector{ii}, tether_file_names{ii}(jj).name);
        comp_tif = imread(comp_tif_dir);
        tether_tif = imread(tether_tif_dir);

        % get removed tether tif
        without_tether_tif = imsubtract(comp_tif, tether_tif);
        without_tether_tif = imbinarize(without_tether_tif);
        without_tether_tif = imfill(without_tether_tif, 'holes');

        % save
        imwrite(without_tether_tif, fullfile(saveTiffPath, comp_file_names{ii}(jj).name));

        % print progress
        if rem(jj, 100) == 0 || jj == length(comp_file_names{ii})
            fprintf(repmat('\b',1,count));
            count = fprintf('folder #: %d, current frame: %s, percentage: %.2f%%\n', ii, comp_file_names{ii}(jj).name, (jj)/(length(comp_file_names{ii}))*100);
        end 
    end
    % copy vid data from Occ tif
    m_file = dir(fullfile(compPathVector{ii}, '*.mat'));
    copyfile(fullfile(compPathVector{ii}, m_file(1).name), fullfile(saveTiffPath,m_file(1).name));
    fprintf('number: %d. substract tether completed, and saved to %s.\n', ii, saveTiffPath);
end
disp('All done!')
