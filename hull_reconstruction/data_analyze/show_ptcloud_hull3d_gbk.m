clc;
clear;
close all;

show_ptcl = 1;
save_video_name = 'pointCloudVideo_nonflipped_axisfixed.avi';

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

% % % % % get LE TE not correct
% % % % rwing_LE = hull3d_file.hull3d.rightwing.hull.LE;
% % % % rwing_TE = hull3d_file.hull3d.rightwing.hull.TE;
% % % % lwing_LE = hull3d_file.hull3d.leftwing.hull.LE;
% % % % lwing_TE = hull3d_file.hull3d.leftwing.hull.TE;


%% save to video
% 创建 VideoWriter 对象
outputVideo = VideoWriter(fullfile(ptcl_path, save_video_name));
outputVideo.FrameRate = 5; % 设置帧率

% 打开 VideoWriter，准备�?存视频
% open(outputVideo);

% �?帧写入视频
disp('save vid')
tot_time = 0;
count = 0;
if ~show_ptcl
    fig = figure('Visible', 'off');
end

for frame = 1:length(frames_ptcl)
    tic
    % 在Figure中绘制三个点云
    try
        % 跳过包�?�NaN值的帧

        % get LE TE
        wingname= {'leftwing','rightwing'};
        [lwing_LE ,lwing_TE] = wing2LETE(hull3d_file.hull3d,hull_file.hull,wingname{1},frame);
        [rwing_LE ,rwing_TE] = wing2LETE(hull3d_file.hull3d,hull_file.hull,wingname{2},frame);

        % 计算身体部分的平移，使中心�?于固定点
        fixed_point = [0, 0, 0];
        body_center = mean(body_ptcl{frame},1);
        translation_vector = fixed_point - body_center;
        mat_trans_vector_b = double(repmat((translation_vector), size(body_ptcl{frame}, 1), 1));
        mat_trans_vector_rw = double(repmat((translation_vector), size(rwing_ptcl{frame}, 1), 1));
        mat_trans_vector_lw = double(repmat((translation_vector), size(lwing_ptcl{frame}, 1), 1));
        % % % mat_trans_vector_body_x = double(repmat((translation_vector), size(body_vectors.X(frame,:), 1), 1));
        % % % mat_trans_vector_body_y = double(repmat((translation_vector), size(body_vectors.Y(frame,:), 1), 1));
        % % % mat_trans_vector_body_z = double(repmat((translation_vector), size(body_vectors.Z(frame,:), 1), 1));
        % % % mat_trans_rwing_span = double(repmat((translation_vector), size(rwing_span(frame,:), 1), 1));
        % % % mat_trans_lwing_span = double(repmat((translation_vector), size(lwing_span(frame,:), 1), 1));
        mat_trans_rwing_LE = double(repmat((translation_vector), size(rwing_LE, 1), 1));
        mat_trans_rwing_TE = double(repmat((translation_vector), size(rwing_TE, 1), 1));
        mat_trans_lwing_LE = double(repmat((translation_vector), size(lwing_LE, 1), 1));
        mat_trans_lwing_TE = double(repmat((translation_vector), size(lwing_TE, 1), 1));
    
        % 平移�?�翅膀
        rwing_ptcl_t{frame} = double(rwing_ptcl{frame}) + mat_trans_vector_rw;
        % 平移身体
        body_ptcl_t{frame} = double(body_ptcl{frame}) + mat_trans_vector_b;
        % 平移左翅膀
        lwing_ptcl_t{frame} = double(lwing_ptcl{frame}) + mat_trans_vector_lw;
        
        % % % % % % 平移身体�?��?
        % % % % % body_vectors.X(frame,:) = double(body_vectors.X(frame,:)) + mat_trans_vector_body_x;
        % % % % % body_vectors.Y(frame,:) = double(body_vectors.Y(frame,:)) + mat_trans_vector_body_y;
        % % % % % body_vectors.Z(frame,:) = double(body_vectors.Z(frame,:)) + mat_trans_vector_body_z;
        % % % % % % 平移span
        % % % % % rwing_span(frame,:) = double(rwing_span(frame,:)) + mat_trans_rwing_span;
        % % % % % lwing_span(frame,:) = double(lwing_span(frame,:)) + mat_trans_lwing_span;
        % 平移LE/TE
        rwing_LE = double(rwing_LE) + mat_trans_rwing_LE;
        rwing_TE = double(rwing_TE) + mat_trans_rwing_TE;
        lwing_LE = double(lwing_LE) + mat_trans_lwing_LE;
        lwing_TE = double(lwing_TE) + mat_trans_lwing_TE;

        scatter3(rwing_ptcl_t{frame}(:, 1), rwing_ptcl_t{frame}(:, 2), rwing_ptcl_t{frame}(:, 3), '.', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8]);    
        hold on;
        scatter3(body_ptcl_t{frame}(:, 1), body_ptcl_t{frame}(:, 2), body_ptcl_t{frame}(:, 3), '.', 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', [0 1 0]);
        scatter3(lwing_ptcl_t{frame}(:, 1), lwing_ptcl_t{frame}(:, 2), lwing_ptcl_t{frame}(:, 3), '.', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8]);
        % plot LE TE
        scatter3(lwing_LE(:, 1), lwing_LE(:, 2), lwing_LE(:, 3), '.', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [0 0 1]);
        scatter3(lwing_TE(:, 1), lwing_TE(:, 2), lwing_TE(:, 3), '.', 'MarkerFaceColor', [0 0 0.6], 'MarkerEdgeColor', [0 0 0.6]);
        scatter3(rwing_LE(:, 1), rwing_LE(:, 2), rwing_LE(:, 3), '.', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
        scatter3(rwing_TE(:, 1), rwing_TE(:, 2), rwing_TE(:, 3), '.', 'MarkerFaceColor', [0.6 0 0], 'MarkerEdgeColor', [0.6 0 0]);

        % % % plot body vector
        %%%quiver3(meanRwin(1),meanRwin(2),meanRwin(3),chord_dir(1),chord_dir(2),chord_dir(3),10);
        % % scatter3(body_vectors.X(frame,1), body_vectors.X(frame,2), body_vectors.X(frame,3), '.', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        % % scatter3(body_vectors.Y(frame,1), body_vectors.Y(frame,2), body_vectors.Y(frame,3), '.', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        % % scatter3(body_vectors.Z(frame,1), body_vectors.Z(frame,2), body_vectors.Z(frame,3), '.', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);

        % % % plot span
        % % scatter3(rwing_span(frame,1), rwing_span(frame, 2), rwing_span(frame, 3), '.', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        % % scatter3(lwing_span(frame,1), lwing_span(frame, 2), lwing_span(frame, 3), '.', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
            
        hold off;

        legend('right wing', 'body', 'left wing');
        % axis([xmin xmax ymin ymax zmin zmax]); % 请替�?�xmin�?xmax�?ymin�?ymax�?zmin和zmax为您希望的范围
        axis([-25 25 -25 25 -25 25]); % 请替�?�xmin�?xmax�?ymin�?ymax�?zmin和zmax为您希望的范围

        num_voxels = size(rwing_ptcl{frame},1)+size(body_ptcl{frame},1)+size(lwing_ptcl{frame},1);
        text(0.1, 0.8, ['Voxels: ', num2str(num_voxels)], 'Units', 'normalized', 'FontSize', 12, 'Color', 'red');
        text(0.8, 0.8, ['Frame: ', num2str(frames_ptcl(frame))], 'Units', 'normalized', 'FontSize', 12, 'Color', 'red');
        view(45, -30); % x 越大越往左转，y 越大越往上转
        if frame == 36
            pause(20)
        end
%         pause(0.5)
        % 将当�?帧写入视频
        % writeVideo(outputVideo, getframe(gcf));
    end

    cur_time = toc;

    % get info
    tot_time = tot_time + cur_time;
    avgtime = tot_time/frame;
    if rem(frame, 10) == 0 || frame == (length(frames_ptcl))
        fprintf(repmat('\b',1,count));
        count = fprintf('cur_frame: %d, percentage: %2f%%,eta: %2fs\n',frame, frame/length(frames_ptcl)*100,avgtime*(length(frames_ptcl)-frame));
    end 
    
end

% 关闭 VideoWriter，完�?视频�?存
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

function [LE_srt,TE_srt,ind,LE_srt2d,TE_srt2d] = CreateBoundary(obj,realC,LE,TE,wingname)
    kcamind = 1;
    for kcam = obj.cameras.all.camvec
        TwoD = Functions.Create2DCoords(realC,obj.cameras.all.DLT_coefs(:,kcam),{LE...
            ,TE,obj.(wingname).hull3d,obj.(wingname).coords.tip,mean(obj.(wingname).hull3d),obj.body.hull3d},1,obj.cameras.all.size_image(1),obj.cameras.all.RotMat_vol);
        [TwoDwingIm,TwoDDownIm,TwoDUpIm] = Split2Dbound(obj,TwoD,obj.wingana.image2D_close);
        % project 3d LE and TE and define the intersecting
        % pixels as LE and TE accordingaly
        BodOnW=0;LETE_same=0;

        if size(intersect(TwoD{3},TwoD{6},'rows'),1)>size(unique(TwoD{3},'rows'),1)*3/4
            % if more than 75% of the wing intersects the body use
            % the LE hull (not boundary)
            wingLE = {LE,LE};
            wingTE = {TE,TE};BodOnW=1;
        else
            [wingLE,wingTE,wingLE2d,wingTE2d] = intersect2D(obj,TwoDwingIm,LE,TE,TwoD{1},TwoD{2},TwoDUpIm,TwoDDownIm);
        end
        % Each wingLE/wingTE has 2 options, each option is constructed
        % using one of the 2D boundary. (each boundary can
        % intersect with the LE as well as the TE. usually the
        % right option will contain more pixels)
        idxs_LETE = [cell2mat(cellfun(@(x) length(x), wingLE, 'UniformOutput', false));...
            cell2mat(cellfun(@(x) length(x), wingTE, 'UniformOutput', false))];
        [~,indM]=max([idxs_LETE(1,1)+idxs_LETE(2,2),idxs_LETE(1,2)+idxs_LETE(2,1)]);
        if indM==1
            op_indLE = 1;op_indTE=2;
        else
            op_indLE = 2;op_indTE=1;
        end
        % calculate the difference between amount of pixels for
        % each option. keep only the indices marked as LE or TE
        % in the 3D hull. (eventually having 3 hulls, one for
        % each camera)

        valDiff(kcamind)=abs(idxs_LETE(1,op_indLE)+idxs_LETE(2,op_indTE)-(idxs_LETE(1,op_indTE)+idxs_LETE(2,op_indLE)));
        LE_srt{kcamind}= [wingLE(op_indLE) wingLE(op_indTE)];
        TE_srt{kcamind}= [wingTE(op_indTE) wingTE(op_indLE)];
        LE_srt2d{kcamind}= [wingLE2d(op_indLE) wingLE2d(op_indTE)];
        TE_srt2d{kcamind}= [wingTE2d(op_indTE) wingTE2d(op_indLE)];
        [~,ind]=sort(valDiff(:),'descend');
        kcamind = kcamind + 1;
    end
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

