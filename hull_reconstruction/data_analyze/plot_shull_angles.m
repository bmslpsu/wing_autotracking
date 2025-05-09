clc;
clear;
close all;

%% angle defination
% phi == stroke
% theta == dev
% psi == rotation

%%
RAWANGLEFlag = 1; % if need the raw data
FILTEDANGLEFlag = 1; % if need the filtered data
FILLEDANGLEFlag = 1; % if need the filled data

[parentPath, ~, ~] = fileparts(pwd);
storage_folder = fullfile(parentPath,'videos');
[shull_name, shull_path] = select_shull(storage_folder);
shull_file = load(fullfile(shull_path,shull_name));
% save_name_pre = [shull_path(end-1) 'cam '];
save_name_pre = [];

body_angles = shull_file.Shull.body.angles;
rwing_angles = shull_file.Shull.rightwing.angles;
lwing_angles = shull_file.Shull.leftwing.angles;
frames = shull_file.Shull.frames;

% check number of nan per stroke
window_size = 40;
AvgPercentNAN(rwing_angles.phi, window_size, 'phi');
AvgPercentNAN(rwing_angles.psi, window_size, 'psi');
AvgPercentNAN(rwing_angles.theta, window_size, 'theta');

if FILLEDANGLEFlag
    figure;
    plot(frames,fillmissing(rwing_angles.phi,'linear'),'.-', 'color','r');
    hold on
    plot(frames,fillmissing(rwing_angles.psi,'linear'),'.-', 'color','g');
    plot(frames,fillmissing(rwing_angles.theta,'linear'),'.-', 'color','b');
    legend('phi==stroke', 'psi==rotation', 'theta==deviation')
    title('filled right wing angle')
    hold off 
    save_name = [save_name_pre 'filled right wing angle.fig'];
    savefig(fullfile(shull_path,save_name));
    fprintf('fig saved to %s\n', fullfile(shull_path,save_name))
    
    figure;
    plot(frames,fillmissing(lwing_angles.phi,'linear'),'.-', 'color','r');
    hold on
    plot(frames,fillmissing(lwing_angles.psi,'linear'),'.-', 'color','g');
    plot(frames,fillmissing(lwing_angles.theta,'linear'),'.-', 'color','b');
    legend('phi==stroke', 'psi==rotation', 'theta==deviation')
    title('filled left wing angle')
    hold off 
    save_name = [save_name_pre 'filled left wing angle.fig'];
    savefig(fullfile(shull_path,save_name));
    fprintf('fig saved to %s\n', fullfile(shull_path,save_name))
    
    figure;
    plot(frames,fillmissing(body_angles.pitch,'linear'),'.-', 'color','r');
    hold on
    plot(frames,fillmissing(body_angles.roll,'linear'),'.-', 'color','g');
    plot(frames,fillmissing(body_angles.yaw,'linear'),'.-', 'color','b');
    legend('pitch', 'roll', 'yaw')
    title('filled body angle')
    hold off 
    save_name = [save_name_pre 'filled body angle.fig'];
    savefig(fullfile(shull_path,save_name));
    fprintf('fig saved to %s\n', fullfile(shull_path,save_name))
end

if RAWANGLEFlag
    figure;
    plot(frames,rwing_angles.phi,'.-', 'color','r');
    hold on
    plot(frames,rwing_angles.psi,'.-', 'color','g');
    plot(frames,rwing_angles.theta,'.-', 'color','b');
    legend('phi==stroke', 'psi==rotation', 'theta==deviation')
    title('raw right wing angle')
    hold off 
    save_name = [save_name_pre 'raw right wing angle.fig'];
    savefig(fullfile(shull_path,save_name));
    fprintf('fig saved to %s\n', fullfile(shull_path,save_name))
    
    figure;
    plot(frames,lwing_angles.phi,'.-', 'color','r');
    hold on
    plot(frames,lwing_angles.psi,'.-', 'color','g');
    plot(frames,lwing_angles.theta,'.-', 'color','b');
    legend('phi==stroke', 'psi==rotation', 'theta==deviation')
    title('raw left wing angle')
    hold off 
    save_name = [save_name_pre 'raw left wing angle.fig'];
    savefig(fullfile(shull_path,save_name));
    fprintf('fig saved to %s\n', fullfile(shull_path,save_name))
    
    figure;
    plot(frames,body_angles.pitch,'.-', 'color','r');
    hold on
    plot(frames,body_angles.roll,'.-', 'color','g');
    plot(frames,body_angles.yaw,'.-', 'color','b');
    legend('pitch', 'roll', 'yaw')
    title('raw body angle')
    hold off 
    save_name = [save_name_pre 'raw body angle.fig'];
    savefig(fullfile(shull_path,save_name));
    fprintf('fig saved to %s\n', fullfile(shull_path,save_name))
end

if FILTEDANGLEFlag
    % phi == stroke
    % theta == dev
    % psi == rotation
    
    cutFreq=500; % used for Wael's video
    SamplingFreq=8000; % used for Wael's video

    % filt right wing angles
    filted_rwing_phi = FiltStrokeAngle(rwing_angles.phi,cutFreq,SamplingFreq);
    filted_rwing_psi = FiltRotationAngle(rwing_angles.psi,cutFreq,SamplingFreq);
    filted_rwing_theta = FiltDevAngle(rwing_angles.theta,cutFreq,SamplingFreq);
    % plot and save angles
    figure;
    plot(frames,filted_rwing_phi,'.-', 'color','r');
    hold on
    plot(frames,filted_rwing_psi,'.-', 'color','g');
    plot(frames,filted_rwing_theta,'.-', 'color','b');
    legend('phi==stroke', 'psi==rotation', 'theta==deviation')
    title('filtered right wing angle')
    hold off 
    save_name = [save_name_pre 'filtered right wing angle.fig'];
    savefig(fullfile(shull_path,save_name));
    fprintf('fig saved to %s\n', fullfile(shull_path,save_name))

    % filt left wing angles
    filted_lwing_phi = FiltStrokeAngle(lwing_angles.phi,cutFreq,SamplingFreq);
    filted_lwing_psi = FiltRotationAngle(lwing_angles.psi,cutFreq,SamplingFreq);
    filted_lwing_theta = FiltDevAngle(lwing_angles.theta,cutFreq,SamplingFreq);
    % plot and save angles
    figure;
    plot(frames,filted_lwing_phi,'.-', 'color','r');
    hold on
    plot(frames,filted_lwing_psi,'.-', 'color','g');
    plot(frames,filted_lwing_theta,'.-', 'color','b');
    legend('phi==stroke', 'psi==rotation', 'theta==deviation')
    title('filtered left wing angle')
    hold off 
    save_name = [save_name_pre 'filtered left wing angle.fig'];
    savefig(fullfile(shull_path,save_name));
    fprintf('fig saved to %s\n', fullfile(shull_path,save_name))

    % % plot and save angles, filt the body angles later
    % figure;
    % plot(frames,body_angles.pitch,'.-', 'color','r');
    % hold on
    % plot(frames,body_angles.roll,'.-', 'color','g');
    % plot(frames,body_angles.yaw,'.-', 'color','b');
    % legend('pitch', 'roll', 'yaw')
    % title('filtered body angle')
    % hold off 
    % save_name = [save_name_pre 'filtered body angle.fig'];
    % savefig(fullfile(shull_path,save_name));
    % fprintf('fig saved to %s\n', fullfile(shull_path,save_name))
end


%% functions
function [filename, filepath] = select_shull(root)
    % select an avi file
    [filename, filepath] = uigetfile(fullfile(root,'shull*.*'), ['select a shull']);
    
    % check if selected a file
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('cancelled');
        error('not selected a file');
    else
        disp(['shull file selected: ', fullfile(filepath, filename)]);
    end
end

%%
function FilteredAngle=FiltStrokeAngle(StrokeAngle,cutFreq,SamplingFreq)
% fill nan data
StrokeAngle = fillmissing(StrokeAngle,'linear');
% start with a hampel filter (using def settings)
FilteredAngle=hampel(StrokeAngle);
% low-pass filtering (variable cut-off freq)
[b_stroke, a_stroke]=butter(6,cutFreq/(SamplingFreq/2));
% filter data 
FilteredAngle=filtfilt(b_stroke,a_stroke,FilteredAngle);
% spline interpolation?
end

%%
function FilteredAngle=FiltRotationAngle(RotationAngle,cutFreq,SamplingFreq)
% fill nan data
RotationAngle = fillmissing(RotationAngle,'linear');
% start with a hampel filter (using def settings)
FilteredAngle=hampel(RotationAngle);
% low-pass filtering (variable cut-off freq)
[b_rotation, a_rotation]=butter(4,1.8*cutFreq/(SamplingFreq/2)); %rotation angle is faster hence a larger cut-off freq
% filter data 
FilteredAngle=filtfilt(b_rotation,a_rotation,FilteredAngle);
% spline interpolation?
end

%%
function FilteredAngle=FiltDevAngle(DevAngle,cutFreq,SamplingFreq)
% fill nan data
DevAngle = fillmissing(DevAngle,'linear');
% start with a hampel filter (using def settings)
FilteredAngle=hampel(DevAngle);
% low-pass filtering (variable cut-off freq)
[b_dev, a_dev]=butter(4,1*cutFreq/(SamplingFreq/2)); %smaller cut-off to reduce noise
% filter data 
FilteredAngle=filtfilt(b_dev,a_dev,FilteredAngle);
% spline interpolation?
end

%%
function AvgPercentNAN(angle_vector, window_size, angle_type)
    nan_count = zeros(1, round(length(angle_vector) / window_size));
    index = 1;
    for i = 1:window_size:length(angle_vector)
        endIndex = min(i + window_size - 1, length(angle_vector));
        if endIndex <= length(angle_vector)
            current_window = angle_vector(i:endIndex);
            nan_count(index) = sum(isnan(current_window));
            index = index + 1;
        end
    end
    avg_nan = mean(nan_count);
    fprintf('The average percentage of nan per wing beat is: %.2f%% for %s\n', (avg_nan/40)*100, angle_type);
end

