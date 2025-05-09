function [return_mat_path,return_mat_file_name] = e_tiff2mat_uint8(wholeFlyOccFilledPathVector,fly_work_dir, fly_folder)

    % Specify the directory containing TIFF images
%     fly_work_dir = 'G:\BMS Lab\wael chapter3';
%     fly_folder = 'fly_1';
    tifDir_o = fullfile(fly_work_dir, fly_folder);
    num_vid = length(wholeFlyOccFilledPathVector);
    metaFileName = 'metaData.mat';

    % convert tiff 2 mat files and save mat files
    return_mat_path = cell(num_vid,1);
    return_mat_file_name = cell(num_vid,1);
    for ii = 1:length(wholeFlyOccFilledPathVector)
        tifDir = wholeFlyOccFilledPathVector{ii};
        % Get a list of all TIFF files in the directory
        tifFiles = dir(fullfile(tifDir, '*.tif'));
        % sort by numbers in the tifFiles
        numericValues = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), {tifFiles.name});
        [sortedValues, sortOrder] = sort(numericValues);
        tifFiles = tifFiles(sortOrder);

        % Initialize a cell array to store image data
        frames = cell(1, numel(tifFiles));

        % Loop through each TIFF file
        count = 0;
        for i = 1:numel(tifFiles)
            % Read the TIFF image
            tifFileName = fullfile(tifDir, tifFiles(i).name);
            img = imread(tifFileName);

            % Store the image data in the cell array
            frames{i} = img;

            if rem(i, 100) == 0 || i == numel(tifFiles)
                fprintf(repmat('\b',1,count));
                count = fprintf('folder #: %d, current frame: %d, percentage: %.2f %%\n', ii, i, (i)/(numel(tifFiles))*100);
            end 
        end
        % convert cell to struct
        frames = cell2struct(frames, 'indIm');

        % Save the cell array to a MAT file
        % [saveMatDir, ~, ~] = fileparts(tifDir_o);
        saveMatDir = fullfile(tifDir_o, 'Mat_Frames_Uint8');
        if ~isfolder(saveMatDir)
                mkdir(saveMatDir)
        end
        sep_path = strsplit(tifDir,filesep);
        matFileName = [sep_path{end}, '.mat'];
        metaData = load(fullfile(tifDir, metaFileName));
        save(fullfile(saveMatDir, matFileName), 'frames', 'metaData');
        return_mat_path{ii} = saveMatDir;
        return_mat_file_name{ii} = matFileName;

        fprintf('Folder: %d completed and saved to %s.\n', ii, saveMatDir);  
    end
    disp(['All TIFF images have been converted and saved to mat']);
end 