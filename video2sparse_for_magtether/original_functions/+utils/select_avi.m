function [filename, filepath] = select_avi(root, cam_num)
    % select an avi file
    [filename, filepath] = uigetfile(fullfile(root,'*.avi'), ['select an avi file from camera #' num2str(cam_num)]);
    
    % check if selected a file
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('cancelled');
        error('not selected a file');
    else
        disp(['avi file selected: ', fullfile(filepath, filename)]);
    end
end