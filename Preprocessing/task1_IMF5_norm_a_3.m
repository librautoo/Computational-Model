close all
clear all

load('yu20150723a_Task-1.mat');

[M,N] = size(data);


%% =========================
%  VMD 参数设置
%% =========================

fs = 256;                 % 采样率 Hz

NumIMFs = 6;              % VMD 分解的 IMF 个数
PenaltyFactor = 2000;     % 带宽惩罚参数 alpha
MaxIterations = 500;      % 最大迭代次数

TargetIMF = 5;            % 这次只提取 IMF5 做平均

Len = 400;                % cue 出现之后 400 个采样点做 VMD
PlotLen = 400;            % 图中显示前 400 个采样点

t = (1:Len) / fs;         % 时间轴，单位：秒


%% =========================
%  人工筛选保留的 trial
%% =========================

keepLeft = [9,14,19,27,35,37,38,39,40,41,42,46,47,54];

keepRight = [8,17,18,22,23,24,29,33,37,39,40,42,43,44];


%% =========================
%  找 cue marker
%% =========================

k1 = 0; 
k2 = 0;

LeftwardsCue = [];
RightwardsCue = [];

for j = 2:N

    % 第8行：视觉提示 marker
    % -1 表示左向 cue，+1 表示右向 cue
    if data(8,j) == -1 && data(8,j-1) == 0
        k1 = k1 + 1;
        LeftwardsCue(k1) = j;
    end

    if data(8,j) == 1 && data(8,j-1) == 0
        k2 = k2 + 1;
        RightwardsCue(k2) = j;
    end

end

M1 = length(LeftwardsCue);   % 左向 cue 个数
M2 = length(RightwardsCue);  % 右向 cue 个数

fprintf('Task-1: Total left cue trials = %d\n', M1);
fprintf('Task-1: Total right cue trials = %d\n', M2);


%% =========================
%  Right cue trials: 提取每个 trial 的 IMF5
%% =========================

right_IMF5_Fz = [];
right_IMF5_F3 = [];
right_IMF5_F4 = [];

muR = 0;

for i = 1:M2

    % 当前右向 cue 的起点
    num = RightwardsCue(i);

    % 防止索引超过数据长度
    if num + Len > N
        continue;
    end

    % trial 计数
    muR = muR + 1;

    % 只分析人工筛选出来的 right trial
    if ~ismember(muR, keepRight)
        continue;
    end

    % cue 出现之后 400 samples 的 EEG
    A = data(1,num+1:num+Len);   % Fz
    B = data(2,num+1:num+Len);   % F3
    C = data(3,num+1:num+Len);   % F4

    fprintf('Right trial %d: Cue = %d, window = [%d, %d]\n', ...
            muR, num, num+1, num+Len);

    % =========================
    %  VMD 分解
    % =========================

    [imfA, residualA, infoA] = vmd_one_channel(A, NumIMFs, PenaltyFactor, MaxIterations);
    [imfB, residualB, infoB] = vmd_one_channel(B, NumIMFs, PenaltyFactor, MaxIterations);
    [imfC, residualC, infoC] = vmd_one_channel(C, NumIMFs, PenaltyFactor, MaxIterations);

    % =========================
    %  提取 IMF5
    % =========================

    right_IMF5_Fz = [right_IMF5_Fz; imfA(:,TargetIMF).'];
    right_IMF5_F3 = [right_IMF5_F3; imfB(:,TargetIMF).'];
    right_IMF5_F4 = [right_IMF5_F4; imfC(:,TargetIMF).'];

end

fprintf('Right trials used for IMF5 average = %d\n', size(right_IMF5_Fz,1));


%% =========================
%  Left cue trials: 提取每个 trial 的 IMF5
%% =========================

left_IMF5_Fz = [];
left_IMF5_F3 = [];
left_IMF5_F4 = [];

muL = 0;

for i = 1:M1

    % 当前左向 cue 的起点
    num = LeftwardsCue(i);

    % 防止索引超过数据长度
    if num + Len > N
        continue;
    end

    % trial 计数
    muL = muL + 1;

    % 只分析人工筛选出来的 left trial
    if ~ismember(muL, keepLeft)
        continue;
    end

    % cue 出现之后 400 samples 的 EEG
    A = data(1,num+1:num+Len);   % Fz
    B = data(2,num+1:num+Len);   % F3
    C = data(3,num+1:num+Len);   % F4

    fprintf('Left trial %d: Cue = %d, window = [%d, %d]\n', ...
            muL, num, num+1, num+Len);

    % =========================
    %  VMD 分解
    % =========================

    [imfA, residualA, infoA] = vmd_one_channel(A, NumIMFs, PenaltyFactor, MaxIterations);
    [imfB, residualB, infoB] = vmd_one_channel(B, NumIMFs, PenaltyFactor, MaxIterations);
    [imfC, residualC, infoC] = vmd_one_channel(C, NumIMFs, PenaltyFactor, MaxIterations);

    % =========================
    %  提取 IMF5
    % =========================

    left_IMF5_Fz = [left_IMF5_Fz; imfA(:,TargetIMF).'];
    left_IMF5_F3 = [left_IMF5_F3; imfB(:,TargetIMF).'];
    left_IMF5_F4 = [left_IMF5_F4; imfC(:,TargetIMF).'];

end

fprintf('Left trials used for IMF5 average = %d\n', size(left_IMF5_Fz,1));


%% =========================
%  对 IMF5 求平均
%% =========================

mean_right_IMF5_Fz = mean(right_IMF5_Fz, 1);
mean_right_IMF5_F3 = mean(right_IMF5_F3, 1);
mean_right_IMF5_F4 = mean(right_IMF5_F4, 1);

mean_left_IMF5_Fz = mean(left_IMF5_Fz, 1);
mean_left_IMF5_F3 = mean(left_IMF5_F3, 1);
mean_left_IMF5_F4 = mean(left_IMF5_F4, 1);


%% =========================
%  Figure 1: Right trials 平均 IMF5
%% =========================

figure;

subplot(3,1,1);
plot(t, mean_right_IMF5_Fz, 'LineWidth', 1.5);
title('Task-1 Right cue trials: mean IMF5 at Fz');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
xlim([0 PlotLen/fs]);
grid on;


subplot(3,1,2);
plot(t, mean_right_IMF5_F3, 'LineWidth', 1.5);
title('Task-1 Right cue trials: mean IMF5 at F3');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
xlim([0 PlotLen/fs]);
grid on;


subplot(3,1,3);
plot(t, mean_right_IMF5_F4, 'LineWidth', 1.5);
title('Task-1 Right cue trials: mean IMF5 at F4');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
xlim([0 PlotLen/fs]);
grid on;

sgtitle(['Task-1 Right cue trials | Mean IMF5 | kept trials = ', mat2str(keepRight)]);


%% =========================
%  Figure 2: Left trials 平均 IMF5
%% =========================

figure;

subplot(3,1,1);
plot(t, mean_left_IMF5_Fz, 'LineWidth', 1.5);
title('Task-1 Left cue trials: mean IMF5 at Fz');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
xlim([0 PlotLen/fs]);
grid on;


subplot(3,1,2);
plot(t, mean_left_IMF5_F3, 'LineWidth', 1.5);
title('Task-1 Left cue trials: mean IMF5 at F3');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
xlim([0 PlotLen/fs]);
grid on;


subplot(3,1,3);
plot(t, mean_left_IMF5_F4, 'LineWidth', 1.5);
title('Task-1 Left cue trials: mean IMF5 at F4');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
xlim([0 PlotLen/fs]);
grid on;

sgtitle(['Task-1 Left cue trials | Mean IMF5 | kept trials = ', mat2str(keepLeft)]);


%% =========================
%  Figure 3: Right vs Left 平均 IMF5 对比
%% =========================

figure;

subplot(3,1,1);
plot(t, mean_right_IMF5_Fz, 'LineWidth', 1.5); hold on;
plot(t, mean_left_IMF5_Fz, 'LineWidth', 1.5);
title('Mean IMF5 comparison at Fz');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
legend('Right cue','Left cue');
xlim([0 PlotLen/fs]);
grid on;


subplot(3,1,2);
plot(t, mean_right_IMF5_F3, 'LineWidth', 1.5); hold on;
plot(t, mean_left_IMF5_F3, 'LineWidth', 1.5);
title('Mean IMF5 comparison at F3');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
legend('Right cue','Left cue');
xlim([0 PlotLen/fs]);
grid on;


subplot(3,1,3);
plot(t, mean_right_IMF5_F4, 'LineWidth', 1.5); hold on;
plot(t, mean_left_IMF5_F4, 'LineWidth', 1.5);
title('Mean IMF5 comparison at F4');
xlabel('Time after cue onset (s)');
ylabel('Amplitude');
legend('Right cue','Left cue');
xlim([0 PlotLen/fs]);
grid on;


sgtitle('Task-1 Mean IMF5 comparison: Right cue vs Left cue');






%% =========================
%  对平均后的 IMF5 做峰值归一化
%  目标：每条平均曲线的最大正峰 = 1
%% =========================

norm_mean_right_IMF5_Fz = normalize_by_positive_peak(mean_right_IMF5_Fz);
norm_mean_right_IMF5_F3 = normalize_by_positive_peak(mean_right_IMF5_F3);
norm_mean_right_IMF5_F4 = normalize_by_positive_peak(mean_right_IMF5_F4);

norm_mean_left_IMF5_Fz = normalize_by_positive_peak(mean_left_IMF5_Fz);
norm_mean_left_IMF5_F3 = normalize_by_positive_peak(mean_left_IMF5_F3);
norm_mean_left_IMF5_F4 = normalize_by_positive_peak(mean_left_IMF5_F4);


%% =========================
%  Figure 4: Right trials 归一化平均 IMF5
%% =========================

figure;

subplot(3,1,1);
plot(t, norm_mean_right_IMF5_Fz, 'LineWidth', 1.5);
title('Task-1 Right cue trials: normalized mean IMF5 at Fz');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;


subplot(3,1,2);
plot(t, norm_mean_right_IMF5_F3, 'LineWidth', 1.5);
title('Task-1 Right cue trials: normalized mean IMF5 at F3');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;



subplot(3,1,3);
plot(t, norm_mean_right_IMF5_F4, 'LineWidth', 1.5);
title('Task-1 Right cue trials: normalized mean IMF5 at F4');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;


sgtitle(['Task-1 Right cue trials | Normalized mean IMF5 | kept trials = ', mat2str(keepRight)]);


%% =========================
%  Figure 5: Left trials 归一化平均 IMF5
%% =========================

figure;

subplot(3,1,1);
plot(t, norm_mean_left_IMF5_Fz, 'LineWidth', 1.5);
title('Task-1 Left cue trials: normalized mean IMF5 at Fz');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;


subplot(3,1,2);
plot(t, norm_mean_left_IMF5_F3, 'LineWidth', 1.5);
title('Task-1 Left cue trials: normalized mean IMF5 at F3');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;



subplot(3,1,3);
plot(t, norm_mean_left_IMF5_F4, 'LineWidth', 1.5);
title('Task-1 Left cue trials: normalized mean IMF5 at F4');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;



sgtitle(['Task-1 Left cue trials | Normalized mean IMF5 | kept trials = ', mat2str(keepLeft)]);


%% =========================
%  Figure 6: Right vs Left 归一化平均 IMF5 对比
%% =========================

figure;

subplot(3,1,1);
plot(t, norm_mean_right_IMF5_Fz, 'LineWidth', 1.5); hold on;
plot(t, norm_mean_left_IMF5_Fz, 'LineWidth', 1.5);
title('Normalized mean IMF5 comparison at Fz');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
legend('Right cue','Left cue');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;



subplot(3,1,2);
plot(t, norm_mean_right_IMF5_F3, 'LineWidth', 1.5); hold on;
plot(t, norm_mean_left_IMF5_F3, 'LineWidth', 1.5);
title('Normalized mean IMF5 comparison at F3');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
legend('Right cue','Left cue');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;



subplot(3,1,3);
plot(t, norm_mean_right_IMF5_F4, 'LineWidth', 1.5); hold on;
plot(t, norm_mean_left_IMF5_F4, 'LineWidth', 1.5);
title('Normalized mean IMF5 comparison at F4');
xlabel('Time after cue onset (s)');
ylabel('Normalized amplitude');
legend('Right cue','Left cue');
xlim([0 PlotLen/fs]);
ylim([-1.2 1.2]);
grid on;


sgtitle('Task-1 Normalized mean IMF5 comparison: Right cue vs Left cue');







%% =========================================================
%  本脚本用到的局部函数：单通道 VMD 分解
%% =========================================================

function [imf, residual, info] = ...
    vmd_one_channel(x, NumIMFs, PenaltyFactor, MaxIterations)

    % 转成列向量
    x = x(:);

    % VMD 分解
    [imf, residual, info] = vmd(x, ...
        NumIMFs = NumIMFs, ...
        PenaltyFactor = PenaltyFactor, ...
        MaxIterations = MaxIterations);

    % 确保 imf 的行数是时间点，列数是 IMF
    if size(imf,1) ~= length(x) && size(imf,2) == length(x)
        imf = imf.';
    end

end



%% =========================================================
%  本脚本用到的局部函数：按最大正峰归一化
%% =========================================================

function x_norm = normalize_by_positive_peak(x)

    peakVal = max(x);

    if peakVal == 0 || isnan(peakVal)
        x_norm = x;
    else
        x_norm = x / peakVal;
    end

end