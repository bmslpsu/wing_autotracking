clc;
clear;

% Specify the path to your AVI file
fly_work_dir = 'G:\BMS Lab\wael chapter3';
fly_folder = 'fly_1';
tifRootPath = fullfile(fly_work_dir, fly_folder);
num_vid = 3;
se = strel('square',2);
se_fill = strel('disk',150);

% get directories of tifs
pathVector = cell(num_vid,1);
tetherPathVector = cell(num_vid,1);
for ii = 1:num_vid
    filepath = utils.select_tif_folder(tifRootPath, ii, 'Wing Occlusion Tifs');
    pathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(tifRootPath, ii, 'Tether Mask Tifs');
    tetherPathVector{ii} = string(filepath);
end

% check the number of tif files and size of the first tiff
occ_file_names = cell(num_vid,1);
tether_file_names = cell(num_vid,1);
for ii = 1:num_vid
    occ_path = pathVector{ii};
    occ_file_names{ii} = dir(fullfile(occ_path, '*.tif'));
    tether_path = tetherPathVector{ii};
    tether_file_names{ii} = dir(fullfile(tether_path, '*.tif'));
    
    % % check number of files
    % if length(occ_file_names) ~= length(tether_file_names)
    %     error('the number of frames are not the same, comp#: %d, tether#: %d', length(occ_file_names), length(tether_file_names))
    % end
    % %check size of first tiff
    % fcomp_size = size(imread(fullfile(occ_path, occ_file_names{ii}(1).name)));
    % ftether_size = size(imread(fullfile(tether_path, tether_file_names{ii}(1).name)));
    % if (fcomp_size) ~= (ftether_size)
    %     error('size of the frames are not the same, comp#: [%d, %d], tether#: [%d, %d]', fcomp_size(1), fcomp_size(2), ftether_size(1), ftether_size(2))
    % end
    utils.compare_numbers_size(occ_path, tether_path, occ_file_names{ii}, tether_file_names{ii})
end

% start filling occlusion
for ii = 1:num_vid
    fprintf('start folder %d\n', ii);

    % get tif saving path
    sep_path = strsplit(pathVector{ii},filesep);
    saveTiffPath = fullfile(tifRootPath, 'Occ_Filled', sep_path{end});
    if ~isfolder(saveTiffPath)
        mkdir(saveTiffPath)
    end
    
    count = 0;
    t_time_elapsed = 0;
    s_time_elapsed = 0;
    for jj = 1:length(occ_file_names{ii})
        tic
        occ_tif_dir = fullfile(pathVector{ii}, occ_file_names{ii}(jj).name);
        tether_tif_dir = fullfile(tetherPathVector{ii}, tether_file_names{ii}(jj).name);
        occ_tif = imread(occ_tif_dir);
        tether_tif = imread(tether_tif_dir);
        tether_tif_morph = bwmorph(tether_tif, 'thicken');

        % get occlusion range
        mask_intersect = and(tether_tif_morph, occ_tif);
        mask_intersect = imdilate(mask_intersect,se); % expand tether mask
        % imshow(mask_intersect)
        CC = bwconncomp(mask_intersect);

        if CC.NumObjects >= 1
            if_fix = 1; % automated fix
        else
            if_fix = 0; % no need for fix
        end

        % Fill the blank
        if if_fix
            labeledImage = bwlabel(mask_intersect);

            % only keep large areas
            stats = regionprops(labeledImage, 'Area'); % get proporty
            threshold = 10;% only keep area greater than threshold
            largeAreaIndices = [stats.Area] > threshold;
            % init an image to save areas want to keep (white pixel = 1)
            resultImage = zeros(size(labeledImage));
            for i = 1:length(stats)
                if largeAreaIndices(i)
                    resultImage(labeledImage == i) = 1;
                end
            end
            resultImage = imbinarize(resultImage);
            % imshow(resultImage)

            % performs a closing to reduce tether occlusion
            mask_intersect = imclose(resultImage, se_fill);  
            % imshow(mask_intersect)
            filled_wing = or(occ_tif, mask_intersect);
            % imshow(filled_wing)
        else
            filled_wing = occ_tif;
        end


        % save 
        imwrite(filled_wing, fullfile(saveTiffPath,occ_file_names{ii}(jj).name));

        % estimate how long need to run
        t_time_elapsed = t_time_elapsed+s_time_elapsed;
        avg_time = t_time_elapsed/jj;
        est_time = (length(occ_file_names{ii})-jj)*(avg_time);

        if rem(jj, 50) == 0 || jj == (length(occ_file_names{ii}))
            fprintf(repmat('\b',1,count));
            count = fprintf('completed vid #: %d, percentage: %.2f %%, eta: %.2f (sec)\n', ii, (jj)/(length(occ_file_names{ii}))*100, est_time);
        end 
        s_time_elapsed = toc;
    end

    % copy vid data from Occ tif
    m_file = dir(fullfile(occ_path, '*.mat'));
    copyfile(fullfile(occ_path, m_file(1).name), fullfile(saveTiffPath,m_file(1).name));
    fprintf('number: %d. Fill wing completed, and saved to %s.\n', ii, saveTiffPath);
end
disp('All done')