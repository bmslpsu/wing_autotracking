function [return_tiff_path] = a1_vid2tiff_Uint8(pathVector, nameVector, bgPathVector, bgNameVector, start_frame, end_frame, actionStartFrame,fly_work_dir, fly_folder)
    % save video to tiff series. Where the tiff is substracted by background.
    % The output is uint8
    
    % Specify the path to your AVI file
    cameraSpeed = 8000;
    % fly_work_dir = 'G:\BMS Lab\wael chapter3';
    % fly_folder = 'fly_1';
    aviRootPath = fullfile(fly_work_dir, fly_folder);
    % num_vid = 3;
    bufferFramesForMask = 15;
    % start_frame = 3000;
    % end_frame = 7000;
    
    % save to metaData
    isFlipped = logical(0);
    % actionStartFrame = 1; % at which frame the action starts (the start_frame is the first frame here)
    
    return_tiff_path = cell(length(nameVector),1);
    % convert vid 2 tiff and save tiff
    for ii = 1:length(nameVector)
        filename = nameVector{ii};
        filepath = pathVector{ii};
        % Create a VideoReader object
        videoObj = VideoReader(fullfile(filepath, filename));
    
        if end_frame > videoObj.NumFrames
            error('End_frame exceeds the number of frames. Please check')
        end
        
        %convert bg tif to mat
        tifFileName = fullfile(bgPathVector{ii},bgNameVector{ii});
        bg = imread(tifFileName);
    
        % get tif saving path
        sep_path = strsplit(filepath,filesep);
        saveTiffPath = fullfile(aviRootPath, 'Comp_TIFF_Uint8', sep_path{end-1});
        if ~isfolder(saveTiffPath)
            mkdir(saveTiffPath)
        end
        return_tiff_path{ii} = saveTiffPath;
    
        % Loop through each frame and save it as a TIFF image
        frame_num = 1;
        frame_count = 1;
        count = 0;
        frame_vec = [];
        while hasFrame(videoObj)
            if frame_num >= (start_frame + bufferFramesForMask) && frame_num <= (end_frame - bufferFramesForMask)
                % Read the current frame
                frame = read(videoObj, frame_num);
                
                % substract bg and orig frame
                subs_bg_frame = imsubtract(bg, frame);
                subs_bg_frame = im2gray(subs_bg_frame);
                subs_bg_frame = imcomplement(subs_bg_frame);
                % binary_frame = imbinarize(subs_bg_frame);
                % binary_frame = imfill(binary_frame,'holes');
                % imshow(binary_frame)
                % imshow(subs_bg_frame)
    
                % find gray boundary
                boundaryThreshold = graythresh(subs_bg_frame);
    
                % get binary frame
                binary_frame = imbinarize(subs_bg_frame,boundaryThreshold+0.1);
                binary_frame = imfill(binary_frame,"holes");
                whiteBGImage = subs_bg_frame;
                whiteBGImage(binary_frame) = 255;
                blackBGImage = imcomplement(whiteBGImage);
                % imshow(whiteBGImage)
                
                % Generate the output TIFF filename (you can customize this)
                outputFileName = sprintf('frame_%d.tif', frame_num);
    
                imwrite(blackBGImage, fullfile(saveTiffPath, outputFileName));
                % save actual frame number for reference
                frame_vec{frame_count} = frame_num;
                frame_count = frame_count + 1;
            end
            
            if rem(frame_num, 100) == 0 || frame_num == (end_frame - bufferFramesForMask)
                fprintf(repmat('\b',1,count));
                count = fprintf('vid #: %d, frame number: %d, percentage: %.2f %%\n', ii, frame_num, (frame_num-start_frame)/(end_frame - start_frame - bufferFramesForMask)*100);
            end 
            if frame_num == end_frame - bufferFramesForMask
                break
            end
            frame_num = frame_num + 1;
        end
    
        % save vid info
        frameSize = [videoObj.Height, videoObj.Width];
        frameRate = cameraSpeed;
        startFrame = actionStartFrame;
        save(fullfile(saveTiffPath, 'metaData'), 'frameSize', 'frameRate', 'bg','isFlipped', 'startFrame','frame_vec');
        fprintf('%d. Conversion %s completed and saved to: %s.\n', ii, filename, saveTiffPath);
    end

    disp('Convert vid2tiff completed.');
end