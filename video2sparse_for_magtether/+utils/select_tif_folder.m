function selectedFolder = select_tif_folder(root, cam_num, file_name)
    % select an tif file
    selectedFolder = uigetdir(root, ['select a ' file_name ' tif directory from camera #' num2str(cam_num)]);
    
    % check if selected a file
    if selectedFolder == 0
        disp('cancelled');
        error('not selected a tif folder');
    else
        disp(['folder selected: ', selectedFolder]);
    end
end