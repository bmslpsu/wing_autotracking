% save video to tiff series. Where the tiff is substracted by background.
% The output is binary
clc;
clear;

% Specify the path to your AVI file
cameraSpeed = 8000;
fly_work_dir = 'G:\BMS Lab\wael chapter3';
fly_folder = 'fly_1';
aviRootPath = fullfile(fly_work_dir, fly_folder);
num_vid = 3;
bufferFramesForMask = 15;
start_frame = 3000;
end_frame = 7000;

% save to metaData
isFlipped = logical(0);
startFrame = 1;


%get directories of avis and bgs
nameVector = cell(num_vid,1);
pathVector = cell(num_vid,1);
bgPathVector = cell(num_vid,1);
bgNameVector = cell(num_vid,1);
for ii = 1:num_vid
    [filename, filepath] = utils.select_avi(aviRootPath, ii);
    nameVector{ii} = string(filename);
    pathVector{ii} = string(filepath);
end
for ii = 1:num_vid
    [filename, filepath] = utils.select_bg(aviRootPath, ii);
    bgPathVector{ii} = string(filepath);
    bgNameVector{ii} = string(filename);
end

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
    saveTiffPath = fullfile(aviRootPath, 'Comp_TIFF', sep_path{end-1});
    if ~isfolder(saveTiffPath)
        mkdir(saveTiffPath)
    end

    % Loop through each frame and save it as a TIFF image
    frame_num = 1;
    count = 0;
    while hasFrame(videoObj)
        if frame_num >= (start_frame + bufferFramesForMask) && frame_num <= (end_frame - bufferFramesForMask)
            % Read the current frame
            frame = read(videoObj, frame_num);
            
            % substract bg and orig frame
            subs_bg_frame = imsubtract(bg, frame);
            subs_bg_frame = im2gray(subs_bg_frame);
            % find gray boundary
            boundaryThreshold = graythresh(subs_bg_frame);
            binary_frame = imbinarize(subs_bg_frame,boundaryThreshold+0.1);
            
            % Generate the output TIFF filename (you can customize this)
            outputFileName = sprintf('frame_%d.tif', frame_num);
            
            % Save the frame as a TIFF image
            binary_frame = imfill(binary_frame,'holes'); % fill noise points

            % keep only the largest white area(fly) because the bottom view
            % is very dirty and the background does not match
            cc = bwconncomp(binary_frame);
            areas = regionprops(cc, 'Area'); % calculate area for each component
            [~, idx] = max([areas.Area]); % get largest part
            largestAreaImage = zeros(size(binary_frame));
            largestAreaImage(cc.PixelIdxList{idx}) = 1;
            largestAreaImage = imbinarize(largestAreaImage);

            % imshow(binary_frame)
            imwrite(largestAreaImage, fullfile(saveTiffPath, outputFileName));
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
    save(fullfile(saveTiffPath, 'metaData'), 'frameSize', 'frameRate', 'bg','isFlipped', 'startFrame');
    fprintf('%d. Conversion %s completed and saved to: %s.\n', ii, filename, saveTiffPath);
end

disp('Conversion completed.');