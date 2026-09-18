close all
clear all

load yu20150723b_Task-1;

[M,N] = size(data);


%% =========================
%  VMD 参数设置
%% =========================

fs = 256;                 % 采样率 Hz

NumIMFs = 6;              % VMD 分解的 IMF 个数
PenaltyFactor = 2000;     % 带宽惩罚参数 alpha
MaxIterations = 500;      % 最大迭代次数

DropHighFreqIMFs = 2;     % 去掉中心频率最高的 2 个 IMF
AddResidual = false;      % 不加回 residual，使降噪更明显

Len = 400;                % cue 出现之后 400 个采样点做 VMD
PlotLen = 400;            % 图中显示前 400 个采样点

t = (1:Len) / fs;         % 时间轴，单位：秒


%% =========================
%  人工筛选保留的 trial
%% =========================

keepLeft = [7,14,19,28,33,36,38,41,42,46,47];

keepRight = [6,8,10,11,16,18,25,30,38,43,49,50];


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
%  Right cue only + VMD
%% =========================

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

    fprintf('Task-1 Right cue trial %d: Cue = %d, window = [%d, %d], duration = %.3f s\n', ...
            muR, num, num+1, num+Len, Len/fs);


    % =========================
    %  VMD 降噪
    % =========================

    [A_vmd, imfA, residualA, infoA, keepIdxA, dropIdxA] = ...
        vmd_denoise_one_channel(A, NumIMFs, PenaltyFactor, MaxIterations, DropHighFreqIMFs, AddResidual);

    [B_vmd, imfB, residualB, infoB, keepIdxB, dropIdxB] = ...
        vmd_denoise_one_channel(B, NumIMFs, PenaltyFactor, MaxIterations, DropHighFreqIMFs, AddResidual);

    [C_vmd, imfC, residualC, infoC, keepIdxC, dropIdxC] = ...
        vmd_denoise_one_channel(C, NumIMFs, PenaltyFactor, MaxIterations, DropHighFreqIMFs, AddResidual);


    % % =========================
    % %  原始 vs VMD 降噪对比图
    % % =========================
    % 
    % figure;
    % 
    % subplot(3,1,1);
    % plot(t, A, 'LineWidth', 0.8); hold on;
    % plot(t, A_vmd, 'LineWidth', 1.2);
    % title(['Task-1 Right cue trial ', num2str(muR), ' cue-aligned EEG at Fz']);
    % legend('Raw','VMD denoised');
    % xlabel('Time after cue onset (s)');
    % ylabel('Amplitude');
    % xlim([0 PlotLen/fs]);
    % grid on;
    % 
    % subplot(3,1,2);
    % plot(t, B, 'LineWidth', 0.8); hold on;
    % plot(t, B_vmd, 'LineWidth', 1.2);
    % title(['Task-1 Right cue trial ', num2str(muR), ' cue-aligned EEG at F3']);
    % legend('Raw','VMD denoised');
    % xlabel('Time after cue onset (s)');
    % ylabel('Amplitude');
    % xlim([0 PlotLen/fs]);
    % grid on;
    % 
    % subplot(3,1,3);
    % plot(t, C, 'LineWidth', 0.8); hold on;
    % plot(t, C_vmd, 'LineWidth', 1.2);
    % title(['Task-1 Right cue trial ', num2str(muR), ' cue-aligned EEG at F4']);
    % legend('Raw','VMD denoised');
    % xlabel('Time after cue onset (s)');
    % ylabel('Amplitude');
    % xlim([0 PlotLen/fs]);
    % grid on;
    % 
    % sgtitle(['Task-1 Right cue trial ', num2str(muR), ...
    %          ' | Cue sample = ', num2str(num)]);


    % =========================
    %  VMD IMF 分解图：Right cue trial
    % =========================

    plot_vmd_imfs(imfA, residualA, keepIdxA, dropIdxA, ...
        ['Task-1 Right cue trial ', num2str(muR), ' Fz'], PlotLen, fs);

    plot_vmd_imfs(imfB, residualB, keepIdxB, dropIdxB, ...
        ['Task-1 Right cue trial ', num2str(muR), ' F3'], PlotLen, fs);

    plot_vmd_imfs(imfC, residualC, keepIdxC, dropIdxC, ...
        ['Task-1 Right cue trial ', num2str(muR), ' F4'], PlotLen, fs);

end

fprintf('Task-1: Total right cue trials checked = %d\n', muR);
fprintf('Task-1: Right cue trials kept for VMD = %s\n', mat2str(keepRight));


%% =========================
%  Left cue only + VMD
%% =========================

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

    fprintf('Task-1 Left cue trial %d: Cue = %d, window = [%d, %d], duration = %.3f s\n', ...
            muL, num, num+1, num+Len, Len/fs);


    % =========================
    %  VMD 降噪
    % =========================

    [A_vmd, imfA, residualA, infoA, keepIdxA, dropIdxA] = ...
        vmd_denoise_one_channel(A, NumIMFs, PenaltyFactor, MaxIterations, DropHighFreqIMFs, AddResidual);

    [B_vmd, imfB, residualB, infoB, keepIdxB, dropIdxB] = ...
        vmd_denoise_one_channel(B, NumIMFs, PenaltyFactor, MaxIterations, DropHighFreqIMFs, AddResidual);

    [C_vmd, imfC, residualC, infoC, keepIdxC, dropIdxC] = ...
        vmd_denoise_one_channel(C, NumIMFs, PenaltyFactor, MaxIterations, DropHighFreqIMFs, AddResidual);


    % =========================
    %  原始 vs VMD 降噪对比图
    % =========================

    figure;

    subplot(3,1,1);
    plot(t, A, 'LineWidth', 0.8); hold on;
    plot(t, A_vmd, 'LineWidth', 1.2);
    title(['Task-1 Left cue trial ', num2str(muL), ' cue-aligned EEG at Fz']);
    legend('Raw','VMD denoised');
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;

    subplot(3,1,2);
    plot(t, B, 'LineWidth', 0.8); hold on;
    plot(t, B_vmd, 'LineWidth', 1.2);
    title(['Task-1 Left cue trial ', num2str(muL), ' cue-aligned EEG at F3']);
    legend('Raw','VMD denoised');
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;

    subplot(3,1,3);
    plot(t, C, 'LineWidth', 0.8); hold on;
    plot(t, C_vmd, 'LineWidth', 1.2);
    title(['Task-1 Left cue trial ', num2str(muL), ' cue-aligned EEG at F4']);
    legend('Raw','VMD denoised');
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;

    sgtitle(['Task-1 Left cue trial ', num2str(muL), ...
             ' | Cue sample = ', num2str(num)]);


    % =========================
    %  VMD IMF 分解图：Left cue trial
    % =========================

    plot_vmd_imfs(imfA, residualA, keepIdxA, dropIdxA, ...
        ['Task-1 Left cue trial ', num2str(muL), ' Fz'], PlotLen, fs);

    plot_vmd_imfs(imfB, residualB, keepIdxB, dropIdxB, ...
        ['Task-1 Left cue trial ', num2str(muL), ' F3'], PlotLen, fs);

    plot_vmd_imfs(imfC, residualC, keepIdxC, dropIdxC, ...
        ['Task-1 Left cue trial ', num2str(muL), ' F4'], PlotLen, fs);

end

fprintf('Task-1: Total left cue trials checked = %d\n', muL);
fprintf('Task-1: Left cue trials kept for VMD = %s\n', mat2str(keepLeft));










%% =========================================================
%  本脚本用到的局部函数 1：单通道 VMD 降噪
%% =========================================================

function [x_clean, imf, residual, info, keepIdx, dropIdx] = ...
    vmd_denoise_one_channel(x, NumIMFs, PenaltyFactor, MaxIterations, DropHighFreqIMFs, AddResidual)

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

    K = size(imf,2);

    allIdx = 1:K;

    % 尝试根据中心频率找出最高频 IMF
    cf = [];

    if isstruct(info) && isfield(info, 'CentralFrequencies')
        cf = info.CentralFrequencies;
    end

    % 有些 MATLAB 版本里 CentralFrequencies 可能是多行，取最后一行
    if ~isempty(cf)

        cf = squeeze(cf);

        if size(cf,1) > 1 && size(cf,2) == K
            cf = cf(end,:);
        elseif size(cf,2) > 1 && size(cf,1) == K
            cf = cf(:,end).';
        else
            cf = cf(:).';
        end

        % 按中心频率从高到低排序
        [~, sortIdx] = sort(cf, 'descend');

        % 删除中心频率最高的若干 IMF
        dropIdx = sortIdx(1:min(DropHighFreqIMFs,K));

    else

        % 如果拿不到中心频率，就默认删除前几个 IMF
        % 通常前几个 IMF 更偏高频
        dropIdx = 1:min(DropHighFreqIMFs,K);

    end

    % 保留剩余 IMF
    keepIdx = setdiff(allIdx, dropIdx);

    % 用保留的 IMF 重构降噪信号
    x_clean = sum(imf(:,keepIdx), 2);

    % 是否加回 residual
    if AddResidual
        residual = residual(:);

        if length(residual) == length(x_clean)
            x_clean = x_clean + residual;
        end
    end

    % 转回行向量，方便和原信号画在一起
    x_clean = x_clean.';

end


%% =========================================================
%  本脚本用到的局部函数 2：画 VMD IMF 分解图
%% =========================================================

function plot_vmd_imfs(imf, residual, keepIdx, dropIdx, figTitle, PlotLen, fs)

    K = size(imf,2);

    % 防止 PlotLen 超过数据长度
    L = min(PlotLen, size(imf,1));

    t = (1:L) / fs;

    figure;

    for kk = 1:K

        subplot(K+1,1,kk);

        plot(t, imf(1:L,kk), 'LineWidth', 0.9);
        grid on;

        if ismember(kk, dropIdx)
            title(['IMF ', num2str(kk), '  Removed']);
        elseif ismember(kk, keepIdx)
            title(['IMF ', num2str(kk), '  Kept']);
        else
            title(['IMF ', num2str(kk)]);
        end

        xlim([0 L/fs]);

    end

    % 最后一行画 residual
    subplot(K+1,1,K+1);

    residual = residual(:);

    if length(residual) >= L
        plot(t, residual(1:L), 'LineWidth', 0.9);
    else
        t_res = (1:length(residual)) / fs;
        plot(t_res, residual, 'LineWidth', 0.9);
    end

    grid on;
    title('Residual');
    xlim([0 L/fs]);
    xlabel('Time after cue onset (s)');

    sgtitle(['VMD decomposition: ', figTitle]);

end