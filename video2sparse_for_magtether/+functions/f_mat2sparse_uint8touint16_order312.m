function f_mat2sparse_uint8touint16_order312(matDirVector, matNameVector,fly_work_dir, fly_folder)

    % initialize
%     fly_work_dir = 'G:\BMS Lab\wael chapter3';
%     fly_folder = 'fly_1';
    matData_root = fullfile(fly_work_dir, fly_folder);
    num_vid = length(matDirVector);
    bit_num = 16;

    % order matNameVector
    charVector = cell(num_vid,1);
    for ii = 1:length(matNameVector)
        charVector{ii} = char(matNameVector{ii});
    end
    charVector = sort(charVector);
    for ii = 1:length(matNameVector)
        matNameVector{ii} = string(charVector{ii});
    end

    % comvert mat files to sparse files and save sparse
    for ii = 1:length(matNameVector)
        frames = cell(length(matNameVector),1);
        matData = load(fullfile(matDirVector{ii},matNameVector{ii}));

        % get new order for the videos.
        switch ii
            case 1
                new_order = 2;
            case 2
                new_order = 3;
            case 3
                new_order = 1;
        end

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
        saveSparseDir = fullfile(matData_root, 'Sparse_Frames_Uint16_Order312');
        if ~isfolder(saveSparseDir)
                mkdir(saveSparseDir)
        end
        fore_name = sprintf('%d_sparse_order312_uint%d_', new_order,bit_num);
        sparseDataName = fore_name+matNameVector{ii};
        metaData = matData.metaData;
        save(fullfile(saveSparseDir,sparseDataName), "frames", 'metaData')

        fprintf('Folder: %d completed and saved to %s.\n', ii, saveSparseDir); 
    end
    disp(['All mat data have been converted and saved to sparse']);
end