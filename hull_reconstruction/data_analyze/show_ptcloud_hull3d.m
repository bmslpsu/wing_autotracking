clc;
clear;
close all;

show_ptcl = 1;
save_video_name = 'ptClVid_flipped_axisf_vec_AMT.mp4';

[parentPath, ~, ~] = fileparts(pwd);
storage_folder = fullfile(parentPath,'videos');
[ptcl_name, ptcl_path] = select_hull3d(storage_folder);
hull3d_file = load(fullfile(ptcl_path,ptcl_name));
[shull_name, shull_path] = select_shull(ptcl_path);
shull_file = load(fullfile(shull_path,shull_name));
[hull_name, hull_path] = select_hull(ptcl_path);
hull_file = load(fullfile(hull_path,hull_name));

% get point clouds
body_ptcl = hull3d_file.hull3d.body.hull;
rwing_ptcl = hull3d_file.hull3d.rightwing.hull.hull3d;
lwing_ptcl = hull3d_file.hull3d.leftwing.hull.hull3d;
frames_ptcl = hull3d_file.hull3d.frames;
% get vectors
body_vectors = shull_file.Shull.body.vectors; % contains x,y,z vectors
rwing_span = shull_file.Shull.rightwing.vectors.span;
lwing_span = shull_file.Shull.leftwing.vectors.span;
% get wing tips
rwing_tips = shull_file.Shull.rightwing.coords.tip;
lwing_tips = shull_file.Shull.leftwing.coords.tip;

% % % % % get LE TE (not correct)
% % % % rwing_LE = hull3d_file.hull3d.rightwing.hull.LE;
% % % % rwing_TE = hull3d_file.hull3d.rightwing.hull.TE;
% % % % lwing_LE = hull3d_file.hull3d.leftwing.hull.LE;
% % % % lwing_TE = hull3d_file.hull3d.leftwing.hull.TE;


%% save to video
% set output video
outputVideo = VideoWriter(fullfile(ptcl_path, save_video_name), 'MPEG-4');
outputVideo.FrameRate = 5; % set frame rate
% start writing
open(outputVideo);

% start generating ptcl
disp('save vid')
tot_time = 0;
count = 0;
if ~show_ptcl
    fig = figure('Visible', 'off');
end

for frame = 1:length(frames_ptcl)
    tic
    %
    try
        % get LE TE
        wingname= {'leftwing','rightwing'};
        [lwing_LE ,lwing_TE] = wing2LETE(hull3d_file.hull3d,hull_file.hull,wingname{1},frame);
        [rwing_LE ,rwing_TE] = wing2LETE(hull3d_file.hull3d,hull_file.hull,wingname{2},frame);

        % get translation
        fixed_point = [0, 0, 0];
        body_center = mean(body_ptcl{frame},1);
        rwing_center = mean(rwing_ptcl{frame},1);
        lwing_center = mean(lwing_ptcl{frame},1);
        translation_vector = fixed_point - body_center;
        mat_trans_vector_b = double(repmat((translation_vector), size(body_ptcl{frame}, 1), 1));
        mat_trans_vector_rw = double(repmat((translation_vector), size(rwing_ptcl{frame}, 1), 1));
        mat_trans_vector_lw = double(repmat((translation_vector), size(lwing_ptcl{frame}, 1), 1));
        mat_trans_rwing_LE = double(repmat((translation_vector), size(rwing_LE, 1), 1));
        mat_trans_rwing_TE = double(repmat((translation_vector), size(rwing_TE, 1), 1));
        mat_trans_lwing_LE = double(repmat((translation_vector), size(lwing_LE, 1), 1));
        mat_trans_lwing_TE = double(repmat((translation_vector), size(lwing_TE, 1), 1));
        mat_trans_rwing_center = double(repmat((translation_vector), size(rwing_center, 1), 1));
        mat_trans_lwing_center = double(repmat((translation_vector), size(lwing_center, 1), 1));
        mat_trans_rwing_tips = double(repmat((translation_vector), size(rwing_tips, 1), 1));
        mat_trans_lwing_tips = double(repmat((translation_vector), size(lwing_tips, 1), 1));
    
        % translate rwing
        rwing_ptcl_t{frame} = double(rwing_ptcl{frame}) + mat_trans_vector_rw;
        % translate body
        body_ptcl_t{frame} = double(body_ptcl{frame}) + mat_trans_vector_b;
        % translate lwing
        lwing_ptcl_t{frame} = double(lwing_ptcl{frame}) + mat_trans_vector_lw;
        % translate LE/TE
        rwing_LE = double(rwing_LE) + mat_trans_rwing_LE;
        rwing_TE = double(rwing_TE) + mat_trans_rwing_TE;
        lwing_LE = double(lwing_LE) + mat_trans_lwing_LE;
        lwing_TE = double(lwing_TE) + mat_trans_lwing_TE;
        % translate wing center
        rwing_center = double(rwing_center) + mat_trans_rwing_center;
        lwing_center = double(lwing_center) + mat_trans_lwing_center;
        % translate wing tips
        ori_rt = rwing_tips;
        rwing_tips = double(rwing_tips) + mat_trans_rwing_tips;
        lwing_tips = double(lwing_tips) + mat_trans_lwing_tips;

        % plot ptcl
        scatter3(rwing_ptcl_t{frame}(:, 1), rwing_ptcl_t{frame}(:, 2), rwing_ptcl_t{frame}(:, 3), '.', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8]);    
        hold on;
        h1 = scatter3(body_ptcl_t{frame}(:, 1), body_ptcl_t{frame}(:, 2), body_ptcl_t{frame}(:, 3), '.', 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0],'DisplayName', 'Body');
        scatter3(lwing_ptcl_t{frame}(:, 1), lwing_ptcl_t{frame}(:, 2), lwing_ptcl_t{frame}(:, 3), '.', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8]);
        % plot LE TE
        h2 = scatter3(lwing_LE(:, 1), lwing_LE(:, 2), lwing_LE(:, 3), '.', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1],'DisplayName', 'lwing_LE');
        h3 = scatter3(lwing_TE(:, 1), lwing_TE(:, 2), lwing_TE(:, 3), '.', 'MarkerFaceColor', [0 0 0.6], 'MarkerEdgeColor', [0 0 0.6],'DisplayName', 'lwing_TE');
        h4 = scatter3(rwing_LE(:, 1), rwing_LE(:, 2), rwing_LE(:, 3), '.', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0],'DisplayName', 'rwing_LE');
        h5 = scatter3(rwing_TE(:, 1), rwing_TE(:, 2), rwing_TE(:, 3), '.', 'MarkerFaceColor', [0.6 0 0], 'MarkerEdgeColor', [0.6 0 0],'DisplayName', 'rwing_TE');

        % % % plot body vector
        quiver3(0,0,0,body_vectors.X(frame,1), body_vectors.X(frame,2), body_vectors.X(frame,3),500,'k');
        quiver3(0,0,0,body_vectors.Y(frame,1), body_vectors.Y(frame,2), body_vectors.Y(frame,3),500,'k');
        quiver3(0,0,0,body_vectors.Z(frame,1), body_vectors.Z(frame,2), body_vectors.Z(frame,3),500,'k');

        % % % plot span
        quiver3(rwing_center(1),rwing_center(2),rwing_center(3), ...
            rwing_span(frame,1)+rwing_center(1), rwing_span(frame,2)+rwing_center(2), rwing_span(frame,3)+rwing_center(3),1200,'r');
        quiver3(lwing_center(1),lwing_center(2),lwing_center(3), ...
            lwing_span(frame,1)+lwing_center(1), lwing_span(frame,2)+lwing_center(2), lwing_span(frame,3)+lwing_center(3),1200,'b');
        % quiver3(rwing_center(1),rwing_center(2),rwing_center(3), ...
            % rwing_tips(frame,1)-rwing_center(1), rwing_tips(frame, 2)-rwing_center(2), rwing_tips(frame, 3)-rwing_center(3), ...
            % 0,'m');
        % quiver3(lwing_center(1),lwing_center(2),lwing_center(3), ...
            % lwing_tips(frame,1)-lwing_center(1), lwing_tips(frame, 2)-lwing_center(2), lwing_tips(frame, 3)-lwing_center(3), ...
            % 0,'m');
        % scatter3(lwing_tips(frame,1),lwing_tips(frame, 2),lwing_tips(frame, 3), '>', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
        % scatter3(lwing_center(1),lwing_center(2),lwing_center(3), '>', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
        % scatter3(rwing_span(frame,1), rwing_span(frame, 2), rwing_span(frame, 3), '.', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
        % scatter3(lwing_span(frame,1), lwing_span(frame, 2), lwing_span(frame, 3), '.', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1]);
        hold off;

        legend([h1,h2,h3,h4,h5], 'DisplayNames',{'Body','lwing_LE','lwing_TE', 'rwing_LE', 'rwing_TE'});
        % axis([xmin xmax ymin ymax zmin zmax]); 
        axis([-35 35 -35 35 -35 35]); 

        num_voxels = size(rwing_ptcl{frame},1)+size(body_ptcl{frame},1)+size(lwing_ptcl{frame},1);
        text(0.1, 0.8, ['Voxels: ', num2str(num_voxels)], 'Units', 'normalized', 'FontSize', 12, 'Color', 'red');
        text(0.8, 0.8, ['Frame: ', num2str(frames_ptcl(frame))], 'Units', 'normalized', 'FontSize', 12, 'Color', 'red');
        view(50, -20); % larger x, turns to left. larger y, turns counter-clockwise
        % if frame == 21
        %     pause('10')
        % end
        % pause(1)
        writeVideo(outputVideo, getframe(gcf));
    end

    cur_time = toc;

    % get info
    tot_time = tot_time + cur_time;
    avgtime = tot_time/frame;
    % if rem(frame, 10) == 0 || frame == (length(frames_ptcl))
    %     fprintf(repmat('\b',1,count));
    %     count = fprintf('cur_frame: %d, percentage: %2f%%,eta: %2fs\n',frame, frame/length(frames_ptcl)*100,avgtime*(length(frames_ptcl)-frame));
    % end 
    
end

% show ptcl or not
if ~show_ptcl
    close(fig);
end
close(outputVideo);
disp('done')


%%
function [filename, filepath] = select_hull3d(root)
    % select an avi file
    [filename, filepath] = uigetfile(fullfile(root,'hull3d*.*'), ['select a hull3d point cloud']);
    
    % check if selected a file
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('cancelled');
        error('not selected a file');
    else
        disp(['hull3d point cloud file selected: ', fullfile(filepath, filename)]);
    end
end

function [filename, filepath] = select_shull(root)
    % select an avi file
    [filename, filepath] = uigetfile(fullfile(root,'shull*.*'), ['select a shull file, not hull3d']);
    
    % check if selected a file
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('cancelled');
        error('not selected a file');
    else
        disp(['shull file selected: ', fullfile(filepath, filename)]);
    end
end

function [filename, filepath] = select_hull(root)
    % select an avi file
    [filename, filepath] = uigetfile(fullfile(root,'hull*.*'), ['select a hull file, not hull3d']);
    
    % check if selected a file
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('cancelled');
        error('not selected a file');
    else
        disp(['hull file selected: ', fullfile(filepath, filename)]);
    end
end

function [LE, TE] = wing2LETE(hull3d,hull,wingname,fr)
    % frm = find(fr == hull.frames);
    wing = hull3d.(wingname).hull.hull3d{fr};
    
    chrd = hull.(wingname).vectors.chord(fr,:);
    chord_dir = (hull.rotmat_EWtoL'* (chrd'))'; % rotate chord to EW axis
    Tip = hull.(wingname).coords.tip(fr,:);
    [LE ,TE] = hull.SplitWin_LETE(chord_dir,Tip,double(wing));
end

%%
function plot_fly_LETE_strk_fig2(hull,hull3d,fr,addSP)
    %UNTITLED6 Summary of this function goes here
    %   Detailed explanation goes here
    ind2mm = 50e-6*1000;
    % frm = find(fr == hull.frames);
    frm = fr;
    wingcell = {'rightwing','leftwing'};
    figure;
    
    colormat = {'r','m';'b','c'};
    hull.body.hull3d = hull3d.body.hull{frm};
    hull.body.vectors.Xtmp =hull.body.vectors.X;
    hull.body.vectors.X =hull.body.vectors.X(frm,:);
    bod = ind2mm*(hull.rotmat_EWtoL *( double([hull3d.body.hull{frm}])'))';
    CMbody = mean(bod);
    realC = hull.real_coord{frm};
    [head_tailCM,headVoxels,tailVoxels] =hull.Xaxis(realC,'plot',0,'save2hull',0); % calculate X body
    tailVoxels = ind2mm*(hull.rotmat_EWtoL * [double(tailVoxels)]')';
    headVoxels = ind2mm*(hull.rotmat_EWtoL * [double(headVoxels)]')';
    
    
    a = ismember(bod,[tailVoxels;headVoxels],'rows');
    middle = bod(a == 0,:) - CMbody;
    
    for kwing = 1:1:2
        wingname = wingcell{kwing};
        wing = double(hull3d.(wingname).hull.hull3d{frm});
        chrd = hull.(wingname).vectors.chord(frm,:);
        chord_dir = (hull.rotmat_EWtoL'* (chrd'))'; % rotate chord to EW axis
        Tip = hull.(wingname).coords.tip(frm,:);
        CMwing = mean(ind2mm*(hull.rotmat_EWtoL*wing')') - CMbody
        [LE ,TE] = hull.SplitWin_LETE(chord_dir,Tip,wing);
        
        hullRealLE = ind2mm*(hull.rotmat_EWtoL *LE')'-CMbody;
        hullRealTE = ind2mm*(hull.rotmat_EWtoL *TE')'-CMbody;
        hull.(wingname).hull3d = hull3d.(wingname).hull.hull3d{frm};
        [cm_coord,wingCM,TE_chrd,LE_chrd]  = hull.EstTip_calcSpanChord(wingname);
        RealLE = ind2mm*(hull.rotmat_EWtoL*LE_chrd')'-CMbody;
        RealTE = ind2mm*(hull.rotmat_EWtoL*TE_chrd')'-CMbody;
        
        hold on;plot3(hullRealLE(:,1),hullRealLE(:,2),hullRealLE(:,3),'.','color',colormat{kwing,1});
        hold on;plot3(hullRealTE(:,1),hullRealTE(:,2),hullRealTE(:,3),'.','color',colormat{kwing,2});
        hold on;plot3(RealLE(:,1),RealLE(:,2),RealLE(:,3),'.k');
        hold on;plot3(RealTE(:,1),RealTE(:,2),RealTE(:,3),'.','color',[0.5 0.5 0.5]);  
    %     hold on;quiver3(CMwing(1),CMwing(2),CMwing(3),chrd(1),chrd(2),chrd(3),1,'color',colormat{kwing,1},'linewidth',2);
    end
    
    tailmm = tailVoxels-CMbody;
    headmm = headVoxels-CMbody;
    cmtail = mean(tailmm);
    hold on;plot3(tailmm(:,1),tailmm(:,2),tailmm(:,3),'.','color',[0,0.9,0],'markersize',5);hold on;
    plot3(headmm(:,1),headmm(:,2),headmm(:,3),'.','color',[0,0.9,0],'markersize',5);
    plot3(middle(:,1),middle(:,2),middle(:,3),'.','color',[0.3 0.67 0.04],'markersize',5);
    
    % hold on;a0 = quiver3(cmtail(1),cmtail(2),cmtail(3),hull.body.vectors.X(1),hull.body.vectors.X(2),hull.body.vectors.X(3),3,'color','k','linewidth',3);
    xlabel('X [mm]');ylabel('Y [mm]');zlabel('Z [mm]');
    axis equal 
    grid on
    % set(gcf,'renderer','painters');
    
    if addSP == 1
    % add wing stroke, get X Y Z wing stroke and plot their plane
    pln =hull.body.vectors.strkPlan(frm,:);
    ybody = hull.body.vectors.Y(frm,:);
    % hold on;a1 = quiver3(0,0,0,ybody(1),ybody(2),ybody(3),1,'linewidth',3,'color','k','MaxHeadSize',0.8);
    Zbody = hull.body.vectors.Z(frm,:);
    % hold on;a2 = quiver3(0,0,0,Zbody(1),Zbody(2),Zbody(3),1,'linewidth',3,'color','k','MaxHeadSize',0.8);
    
    Ystrk = ybody;
    Xstrk = cross(Ystrk,pln);
    sqr = [Xstrk-Ystrk;Xstrk+Ystrk;-Xstrk+Ystrk;-Xstrk-Ystrk];
    X = sqr(:,1);
    Y = sqr(:,2);
    Z = -1/pln(3)*(pln(1)*X + pln(2)*Y);
    
    
    ptch = patch('XData',1.5*X,'YData',1.5*Y,'ZData',1.5*Z) ;
    ptch.FaceColor = [0.5 0 1];
    ptch.FaceAlpha = 0.1;
    hold on;quiver3(0,0,0,pln(1),pln(2),pln(3),2,'linewidth',2,'color',[0.8 0 1]);
    end
    legend off

end

