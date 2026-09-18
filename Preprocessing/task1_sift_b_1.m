close all
clear all

load yu20150723b_Task-1;

[M,N] = size(data);


%% =========================
%  基本参数设置
%% =========================

fs = 256;          % 采样率 Hz
Len = 1000;         % cue 出现之后截取 1000 samples
PlotLen = 1000;     % 图中显示 1000 samples

t = (1:Len) / fs;  % 时间轴，单位：秒


%% =========================
%  找 cue 和 action marker
%% =========================

k1 = 0;
k2 = 0;

l1 = 0;
l2 = 0;

LeftwardsCue = [];
RightwardsCue = [];

LeftAct = [];
RightAct = [];

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

    % 第9行：动作反应 marker
    % Task-1: -1 表示左向 action，+1 表示右向 action
    if data(9,j) == -1 && data(9,j-1) == 0
        l1 = l1 + 1;
        LeftAct(l1) = j;
    end

    if data(9,j) == 1 && data(9,j-1) == 0
        l2 = l2 + 1;
        RightAct(l2) = j;
    end

end

M1 = length(LeftwardsCue);   % 左向 cue 个数
M2 = length(RightwardsCue);  % 右向 cue 个数

fprintf('Task-1: Total left cue trials = %d\n', M1);
fprintf('Task-1: Total right cue trials = %d\n', M2);
fprintf('Task-1: Total left actions = %d\n', length(LeftAct));
fprintf('Task-1: Total right actions = %d\n', length(RightAct));


%% =========================
%  Right cue only: 画原始 EEG + 标出 action
%% =========================

muR = 0;

for i = 1:M2

    % 当前右向 cue 的起点
    num = RightwardsCue(i);

    % 防止索引超过数据长度
    if num + Len > N
        continue;
    end

    % cue 出现之后 400 samples 的原始 EEG
    A = data(1,num+1:num+Len);   % Fz
    B = data(2,num+1:num+Len);   % F3
    C = data(3,num+1:num+Len);   % F4

    muR = muR + 1;

    % 找当前 right cue 后面的第一个 right action
    possibleActs = RightAct(RightAct > num);

    if ~isempty(possibleActs)
        APoi = possibleActs(1);
        RT_samples = APoi - num;
        RT_sec = RT_samples / fs;
    else
        APoi = NaN;
        RT_samples = NaN;
        RT_sec = NaN;
    end

    fprintf('Right cue trial %d: Cue = %d, Action = %d, RT = %d samples, RT = %.3f s\n', ...
            muR, num, APoi, RT_samples, RT_sec);

    figure;

    subplot(3,1,1);
    plot(t, A, 'LineWidth', 1.0);
    title(['Task-1 Right cue trial ', num2str(muR), ' raw EEG at Fz']);
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;
    hold on;
    if ~isnan(RT_samples) && RT_samples <= PlotLen
        xline(RT_sec, 'r--', 'Action', 'LineWidth', 1.2);
    end

    subplot(3,1,2);
    plot(t, B, 'LineWidth', 1.0);
    title(['Task-1 Right cue trial ', num2str(muR), ' raw EEG at F3']);
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;
    hold on;
    if ~isnan(RT_samples) && RT_samples <= PlotLen
        xline(RT_sec, 'r--', 'Action', 'LineWidth', 1.2);
    end

    subplot(3,1,3);
    plot(t, C, 'LineWidth', 1.0);
    title(['Task-1 Right cue trial ', num2str(muR), ' raw EEG at F4']);
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;
    hold on;
    if ~isnan(RT_samples) && RT_samples <= PlotLen
        xline(RT_sec, 'r--', 'Action', 'LineWidth', 1.2);
    end

    if ~isnan(RT_samples) && RT_samples > PlotLen
        sgtitle(['Task-1 Right cue trial ', num2str(muR), ...
                 ' | Cue sample = ', num2str(num), ...
                 ' | Action outside display: RT = ', num2str(RT_sec, '%.3f'), ' s']);
    else
        sgtitle(['Task-1 Right cue trial ', num2str(muR), ...
                 ' | Cue sample = ', num2str(num), ...
                 ' | RT = ', num2str(RT_sec, '%.3f'), ' s']);
    end

end

fprintf('Task-1: Total right cue trials plotted = %d\n', muR);


%% =========================
%  Left cue only: 画原始 EEG + 标出 action
%% =========================

muL = 0;

for i = 1:M1

    % 当前左向 cue 的起点
    num = LeftwardsCue(i);

    % 防止索引超过数据长度
    if num + Len > N
        continue;
    end

    % cue 出现之后 400 samples 的原始 EEG
    A = data(1,num+1:num+Len);   % Fz
    B = data(2,num+1:num+Len);   % F3
    C = data(3,num+1:num+Len);   % F4

    muL = muL + 1;

    % 找当前 left cue 后面的第一个 left action
    possibleActs = LeftAct(LeftAct > num);

    if ~isempty(possibleActs)
        APoi = possibleActs(1);
        RT_samples = APoi - num;
        RT_sec = RT_samples / fs;
    else
        APoi = NaN;
        RT_samples = NaN;
        RT_sec = NaN;
    end

    fprintf('Left cue trial %d: Cue = %d, Action = %d, RT = %d samples, RT = %.3f s\n', ...
            muL, num, APoi, RT_samples, RT_sec);

    figure;

    subplot(3,1,1);
    plot(t, A, 'LineWidth', 1.0);
    title(['Task-1 Left cue trial ', num2str(muL), ' raw EEG at Fz']);
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;
    hold on;
    if ~isnan(RT_samples) && RT_samples <= PlotLen
        xline(RT_sec, 'r--', 'Action', 'LineWidth', 1.2);
    end

    subplot(3,1,2);
    plot(t, B, 'LineWidth', 1.0);
    title(['Task-1 Left cue trial ', num2str(muL), ' raw EEG at F3']);
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;
    hold on;
    if ~isnan(RT_samples) && RT_samples <= PlotLen
        xline(RT_sec, 'r--', 'Action', 'LineWidth', 1.2);
    end

    subplot(3,1,3);
    plot(t, C, 'LineWidth', 1.0);
    title(['Task-1 Left cue trial ', num2str(muL), ' raw EEG at F4']);
    xlabel('Time after cue onset (s)');
    ylabel('Amplitude');
    xlim([0 PlotLen/fs]);
    grid on;
    hold on;
    if ~isnan(RT_samples) && RT_samples <= PlotLen
        xline(RT_sec, 'r--', 'Action', 'LineWidth', 1.2);
    end

    if ~isnan(RT_samples) && RT_samples > PlotLen
        sgtitle(['Task-1 Left cue trial ', num2str(muL), ...
                 ' | Cue sample = ', num2str(num), ...
                 ' | Action outside display: RT = ', num2str(RT_sec, '%.3f'), ' s']);
    else
        sgtitle(['Task-1 Left cue trial ', num2str(muL), ...
                 ' | Cue sample = ', num2str(num), ...
                 ' | RT = ', num2str(RT_sec, '%.3f'), ' s']);
    end

end

fprintf('Task-1: Total left cue trials plotted = %d\n', muL);