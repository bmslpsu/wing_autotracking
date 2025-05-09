function [filename, filepath] = select_mat(root, cam_num)
    % select an avi file
    [filename, filepath] = uigetfile(fullfile(root,'*.mat'), ['select a mat file from camera #' num2str(cam_num)]);
    
    % check if selected a file
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('cancelled');
        error('not selected a mat file');
    else
        disp(['mat file selected: ', fullfile(filepath, filename)]);
    end
end