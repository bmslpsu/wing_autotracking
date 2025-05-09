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
wholeFlyPathVector = cell(num_vid,1);
for ii = 1:num_vid
    filepath = utils.select_tif_folder(fullfile(tifRootPath, 'wing_mask_Uint8'), ii, 'Wing Occlusion Uint8 Tifs');
    pathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(fullfile(tifRootPath, 'tether_mask'), ii, 'Tether Mask Tifs');
    tetherPathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    filepath = utils.select_tif_folder(fullfile(tifRootPath, 'whole_fly_mask_Uint8'), ii, 'Whole Fly Occlusion Uint8 Tifs');
    wholeFlyPathVector{ii} = string(filepath);
end

% check the number of tif files and size of the first tiff
occ_file_names = cell(num_vid,1);
tether_file_names = cell(num_vid,1);
for ii = 1:num_vid
    occ_path = pathVector{ii};
    occ_file_names{ii} = dir(fullfile(occ_path, '*.tif'));
    tether_path = tetherPathVector{ii};
    tether_file_names{ii} = dir(fullfile(tether_path, '*.tif'));
    utils.compare_numbers_size(occ_path, tether_path, occ_file_names{ii}, tether_file_names{ii})
end
fly_file_names = cell(num_vid,1);
for ii = 1:num_vid
    occ_path = pathVector{ii};
    occ_file_names{ii} = dir(fullfile(occ_path, '*.tif'));
    fly_path = wholeFlyPathVector{ii};
    fly_file_names{ii} = dir(fullfile(fly_path, '*.tif'));
    utils.compare_numbers_size(occ_path, fly_path, occ_file_names{ii}, fly_file_names{ii})
end

% start filling occlusion
parfor ii = 1:num_vid
    print_cont = fprintf('start folder %d\n', ii);

    % get tif saving path
    sep_path = strsplit(pathVector{ii},filesep);
    saveTiffPath = fullfile(tifRootPath, 'Fly_Occ_Filled', sep_path{end});
    if ~isfolder(saveTiffPath)
        mkdir(saveTiffPath)
    end
    
    count = 0;
    t_time_elapsed = 0;
    s_time_elapsed = 0;
    for jj = 1:length(occ_file_names{ii})
        tic
        occ_tif_dir = fullfile(pathVector{ii}, occ_file_names{ii}(jj).name);
        fly_tif_dir = fullfile(wholeFlyPathVector{ii}, fly_file_names{ii}(jj).name);
        tether_tif_dir = fullfile(tetherPathVector{ii}, tether_file_names{ii}(jj).name);
        occ_tif = imread(occ_tif_dir);
        fly_tif = imread(fly_tif_dir);
        tether_tif = imread(tether_tif_dir);
        tether_tif_morph = bwmorph(tether_tif, 'thicken');

        % get occlusion range
%         occ_bin = imbinarize(occ_tif); % intersection of tether and wing
%         mask_intersect = and(tether_tif_morph, occ_bin);
        fly_bin = imbinarize(fly_tif); % intersection of tether and fly
        mask_intersect = and(tether_tif_morph, fly_bin);
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
            % imshow(resultImage)

            % performs a closing to reduce tether occlusion
            mask_intersect = imclose(resultImage, se_fill);

            % get mean gray value for wings
            blackPixels = (occ_tif == 0);
            modifiedImage = double(occ_tif);
            modifiedImage(blackPixels) = NaN;
            averageGray = mean(modifiedImage(:), 'omitnan');
            
            mask_intersect = uint8(mask_intersect) * 255;
            % change non black pixels to averageGray
            nonBlackPixels = (mask_intersect ~= 0);
            mask_intersect(nonBlackPixels) = averageGray;
            % imshow(mask_intersect)

            % imshow(mask_intersect)
            img_mask = (fly_bin > 0)&(mask_intersect > 0);
            mask_intersect(img_mask) = 0;
            % filled_wing = imadd(occ_tif, mask_intersect);
            % imshow(filled_wing)

            % combine different(filled) area of wing only and whole fly
            % img_mask = (filled_wing == fly_tif);
            % fly_tif = fly_tif .* uint8(~img_mask) + filled_wing .* uint8(img_mask);
            fly_tif = imadd(fly_tif, mask_intersect);
        end

        % save 
        imwrite(fly_tif, fullfile(saveTiffPath,occ_file_names{ii}(jj).name));

        % estimate how long need to run
        t_time_elapsed = t_time_elapsed+s_time_elapsed;
        avg_time = t_time_elapsed/jj;
        est_time = (length(occ_file_names{ii})-jj)*(avg_time);

        if rem(jj, 10) == 0 || jj == (length(occ_file_names{ii}))
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