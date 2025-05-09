function [filename, filepath] = select_bg(root, cam_num)
    % select an avi file
    [filename, filepath] = uigetfile(fullfile(root,'*.tif'), ['select a bg(tif) file from camera #' num2str(cam_num)]);
    
    % check if selected a file
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('cancelled');
        error('not selected a background');
    else
        disp(['bg file selected: ', fullfile(filepath, filename)]);
    end
end