clc;
clear;

% Specify the directory containing TIFF images
fly_work_dir = 'L:\BMS_Lab\wael_chapter3';
fly_folder = 'fly_5';
tifDir_o = fullfile(fly_work_dir, fly_folder);
num_vid = 3;
metaFileName = 'metaData.mat';

% get directories of tiffs
folderVector = cell(num_vid,1);
for ii = 1:num_vid
    selectedFolder = utils.select_tif_folder(fullfile(tifDir_o, 'Whole_Fly_Occ_Filled_Uint8'), ii, 'Ready2Convert2Mat');
    folderVector{ii} = string(selectedFolder);
end

% convert tiff 2 mat files and save mat files
for ii = 1:length(folderVector)
    tifDir = folderVector{ii};
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
    
    fprintf('Folder: %d completed and saved to %s.\n', ii, saveMatDir);  
end
disp(['All TIFF images have been converted and saved to mat']);