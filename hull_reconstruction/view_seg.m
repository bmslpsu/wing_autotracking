clc;
clear;

% init and get segmentation file
mov = 1;
mov_dir = fullfile(pwd,'videos',['mov' num2str(mov)],'Segmentation');
seg_name = ['mov' num2str(mov) '_seg.mat'];
seg_file_dir = fullfile(mov_dir,seg_name);
seg_file = load(seg_file_dir);
image_size = [512, 768];
% image_data = zeros(image_size);

% figure('Visible', 'on');
% size(seg_file.seg.wing1);

for c = 1:min(min(size(seg_file.seg.wing1,2),size(seg_file.seg.wing2,2)),size(seg_file.seg.body,2))
    save_dir = fullfile(mov_dir,'seg_pic',['cam' num2str(c)]);
    mkdir(save_dir);
    % read seg
    wing_1 = seg_file.seg.wing1{c};  % get wing 1 seg pixels
    wing_2 = seg_file.seg.wing2{c}; % get wing 2 seg pixels
    body = seg_file.seg.body{c}; % get body seg pixels
    max_frame = min(min(length(wing_1),length(wing_2)),length(body));
    fprintf('cam:%d',c);
    count = 0;
    for i = 1:max_frame
        % plot figures
        if ~isempty(wing_1(i).indIm) || ~isempty(wing_2(i).indIm) || ~isempty(body(i).indIm)
            if ~isempty(wing_1(i).indIm)
                wing_1i = wing_1(i).indIm;
                scatter(wing_1i(:,2), wing_1i(:, 1), 'r', 'filled');
                hold on
            end
            if ~isempty(wing_2(i).indIm)
                wing_2i = wing_2(i).indIm;
                scatter(wing_2i(:,2), wing_2i(:, 1), 'b', 'filled');
                hold on
            end
            if ~isempty(body(i).indIm)
                bodyi = body(i).indIm;
                scatter(bodyi(:,2), bodyi(:, 1), 'g', 'filled');
                hold on
            end

            hold off
            % gcf = flipud(gcf);
            axis([0 image_size(2)  0 image_size(1)]); 

            saveas(gcf,fullfile(save_dir,['frame_' num2str(i) '.tif']));
        end
        if rem(i,10)==0 || i == max_frame
            fprintf(repmat('\b',1,count));
            count=fprintf('current image:%d, percentage:%.2f%%\n',i, (i/max_frame)*100);
        end
    end
end
disp('conversion complete')

