function  [return_wing_tiff_dir,return_whole_fly_tiff_dir]= c1_GetWingandWholeMaskTIFF_Uint8(compPathVector, tetherPathVector, bodyPathVector,fly_work_dir, fly_folder)
    % Specify the path to your AVI file
    % fly_work_dir = 'G:\BMS Lab\wael chapter3';
    % fly_folder = 'fly_1';
    tifRootPath = fullfile(fly_work_dir, fly_folder);
    num_vid = length(compPathVector);
    NeedWholeFly = 1;
    NeedWingOnly = 1;
    
    % check the number of tif files and size of the first tiff
    comp_file_names = cell(num_vid,1);
    tether_file_names = cell(num_vid,1);
    body_file_names = cell(num_vid,1);
    for ii = 1:num_vid
        comp_path = compPathVector{ii};
        comp_file_names{ii} = dir(fullfile(comp_path, '*.tif'));
        tether_path = tetherPathVector{ii};
        tether_file_names{ii} = dir(fullfile(tether_path, '*.tif'));
        body_path = tetherPathVector{ii};
        body_file_names{ii} = dir(fullfile(body_path, '*.tif'));
    
        utils.compare_numbers_size(comp_path, tether_path, comp_file_names{ii}, tether_file_names{ii})
        utils.compare_numbers_size(comp_path, body_path, comp_file_names{ii}, body_file_names{ii})
    end
    
    return_wing_tiff_dir = cell(num_vid,1);
    return_whole_fly_tiff_dir = cell(num_vid,1);
    % substract tether and/or body from comp tif
    for ii = 1:num_vid
    
        % get wing tif saving path
        sep_path = strsplit(compPathVector{ii},filesep);
        saveWingTiffPath = fullfile(tifRootPath, 'wing_mask_Uint8', sep_path{end});
        if ~isfolder(saveWingTiffPath)
            mkdir(saveWingTiffPath)
        end
        return_wing_tiff_dir{ii} = saveWingTiffPath;
    
        % get whole fly tif saving path
        sep_path = strsplit(compPathVector{ii},filesep);
        saveWholeTiffPath = fullfile(tifRootPath, 'whole_fly_mask_Uint8', sep_path{end});
        if ~isfolder(saveWholeTiffPath)
            mkdir(saveWholeTiffPath)
        end
        return_whole_fly_tiff_dir{ii} = saveWholeTiffPath;
        
        count = 0;
        t_time_elapsed = 0;
        s_time_elapsed = 0;
        for jj = 1:length(comp_file_names{ii})
            tic
            comp_tif_dir = fullfile(compPathVector{ii}, comp_file_names{ii}(jj).name);
            tether_tif_dir = fullfile(tetherPathVector{ii}, tether_file_names{ii}(jj).name);
            body_tif_dir = fullfile(bodyPathVector{ii}, body_file_names{ii}(jj).name);
            comp_tif = imread(comp_tif_dir);
            tether_tif = imread(tether_tif_dir);
            tether_uint8 = uint8(tether_tif) * 255;
            body_tif = imread(body_tif_dir);
            body_uint8 = uint8(body_tif) * 255;
    
            % get removed body and tether tif
            wo_tether_tif = imsubtract(comp_tif, tether_uint8);
            % wo_tether_tif = imbinarize(wo_tether_tif);
            % wo_tether_tif = imfill(wo_tether_tif, 'holes');
    
            % save
            if NeedWholeFly
                imwrite(wo_tether_tif, fullfile(saveWholeTiffPath, comp_file_names{ii}(jj).name));
            end
            if NeedWingOnly
                wo_tether_body_tif = imsubtract(wo_tether_tif, body_uint8);
                % imshow(wo_tether_body_tif)
                % wo_tether_body_tif = imbinarize(wo_tether_body_tif);
                % wo_tether_body_tif = imfill(wo_tether_body_tif, 'holes');
                imwrite(wo_tether_body_tif, fullfile(saveWingTiffPath, comp_file_names{ii}(jj).name));
            end
    
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
        m_file = dir(fullfile(compPathVector{ii}, '*.mat'));
        copyfile(fullfile(compPathVector{ii}, m_file(1).name), fullfile(saveWingTiffPath,m_file(1).name));
        fprintf('number: %d. wing mask get, and saved to %s.\n', ii, saveWingTiffPath);
        fprintf('number: %d. whole wing mask get, and saved to %s.\n', ii, saveWholeTiffPath);
    end
    disp('All done! Get wing and whole mask tiff')
end