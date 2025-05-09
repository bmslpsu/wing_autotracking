clc;
clear;

% initialize
fly_work_dir = 'G:\BMS Lab\wael chapter3';
fly_folder = 'fly_1';
matData_root = fullfile(fly_work_dir, fly_folder);
num_vid = 3;
bit_num = 16;

% load mat files and directories in to cells
matNameVector = cell(num_vid,1);
matDirVector = cell(num_vid,1);
for ii = 1:num_vid
    [filename, filepath] = utils.select_mat(fullfile(matData_root, 'Mat_Frames_Uint8'), ii);
    matNameVector{ii} = string(filename);
    matDirVector{ii} = string(filepath);
end

% comvert mat files to sparse files and save sparse
for ii = 1:length(matNameVector)
    frames = cell(length(matNameVector),1);
    matData = load(fullfile(matDirVector{ii},matNameVector{ii}));

    count = 0;
    for jj = 1:length(matData.frames)
        if bit_num == 16
            matData.frames(jj).indIm = utils.convert_uint8touint16(matData.frames(jj).indIm);
        end 
        sparseData = sparse(double(matData.frames(jj).indIm));
        [row, col, values] = find(sparseData);
        frames{jj} = [row col values];

        if rem(jj, 100) == 0 || jj == length(matData.frames)
            fprintf(repmat('\b',1,count));
            count = fprintf('mat #: %d, current frame: %d, percentage: %.2f %%\n', ii, jj, (jj)/(length(matData.frames))*100);
        end         
    end
    % convert cell to struct
    frames = cell2struct(frames, 'indIm',size(frames{1},2));
    % do not forget to convert bg
    if bit_num == 16
        matData.metaData.bg = utils.convert_uint8touint16(matData.metaData.bg);
    end 
    % [saveSparseDir, ~, ~] = fileparts(matData_root);
    saveSparseDir = fullfile(matData_root, 'Sparse_Frames_Uint16');
    if ~isfolder(saveSparseDir)
            mkdir(saveSparseDir)
    end
    fore_name = sprintf('sparse_uint%d_', bit_num);
    sparseDataName = fore_name+matNameVector{ii};
    metaData = matData.metaData;
    save(fullfile(saveSparseDir,sparseDataName), "frames", 'metaData')

    fprintf('Folder: %d completed and saved to %s.\n', ii, saveSparseDir); 
end
disp(['All mat data have been converted and saved to sparse']);