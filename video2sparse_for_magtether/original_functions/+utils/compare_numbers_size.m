function []= compare_numbers_size(path_1, path_2, file_names_1, file_names_2)

    % check number of files
    if length(file_names_1) ~= length(file_names_2)
        error('the number of frames are not the same, comp#: %d, tether#: %d', length(file_names_1), length(file_names_2))
    end
    %check size of first tiff
    size_1 = size(imread(fullfile(path_1, file_names_1(1).name)));
    size_2 = size(imread(fullfile(path_2, file_names_2(1).name)));
    if (size_1) ~= (size_2)
        error('size of the frames are not the same, comp#: [%d, %d], tether#: [%d, %d]', size_1(1), size_1(2), size_2(1), size_2(2))
    end
end