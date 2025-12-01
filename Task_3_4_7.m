%% Declaration of algorithm parameters

y_nominal = Innov_nominal(1:N,11); % measurement data
y = Innov_attack(1:N,11); % measurement data
y(1) = 0; % to remove the initial spike and focus on when the KF converges
% y = y_nominal;

theta_0 = [mean(y_nominal(1:N,1))];
sigma_0 = [std(y_nominal(1:N,1))];

% Thresholds for positive and negative tests
threshold_pos = [500]; % Example thresholds for positive tests for each measurement
threshold_neg = [-200]; % Example thresholds for negative tests for each measurement

% Initialize variables to count alarms
n_alarm_pos = zeros(1, size(y, 2)); % number of alarms for positive test
n_alarm_neg = zeros(1, size(y, 2)); % number of alarms for negative test

% Initialize arrays to store time indices of alarms
k_alarm_pos = cell(1, size(y, 2)); % time indices k for alarms of positive test
k_alarm_neg = cell(1, size(y, 2)); % time indices k for alarms of negative test

% Initialize CUSUM variables
g_pos = zeros(size(y));
g_neg = zeros(size(y));
s = (y - theta_0) ./ sigma_0; 

%% CUSUM Test
for m = 1:size(y, 2)
    for k = 1:size(y, 1)-1
        g_pos(k+1, m) = max(0, g_pos(k, m) + s(k, m))-9; % positive test
        g_neg(k+1, m) = min(0, g_neg(k, m) + s(k, m))+9; % negative test

        % Check for positive alarms
        if g_pos(k+1, m) > threshold_pos(m)
            n_alarm_pos(m) = n_alarm_pos(m) + 1;
            k_alarm_pos{m} = [k_alarm_pos{m} k+1];
            g_pos(k+1, m) = 0; % reset
        end

        % Check for negative alarms
        if g_neg(k+1, m) < threshold_neg(m)
            n_alarm_neg(m) = n_alarm_neg(m) + 1;
            k_alarm_neg{m} = [k_alarm_neg{m} k+1];
            g_neg(k+1, m) = 0; % reset
        end
    end

    % Plotting
    titles = {'x(13)=b_A_x', 'x(14)=b_A_y', 'x(15)=b_A_z', 'x(16)=b_p', 'x(17)=b_q', 'x(18)=b_r'};
    figure;
    
    % Subplot 1: Plotting g_pos and g_neg with alarms
    subplot(2, 1, 1);
    p1 = plot(g_pos(:, m));
    hold on;
    p2 = plot(g_neg(:, m));
    p3 = []; % Initialize p3
    p4 = []; % Initialize p4
    for i = 1:length(k_alarm_pos{m})
        p3=plot([k_alarm_pos{m}(i) k_alarm_pos{m}(i)], [-20 20], 'r--');
    end
    for i = 1:length(k_alarm_neg{m})
        p4=plot([k_alarm_neg{m}(i) k_alarm_neg{m}(i)], [-20 20], 'r-.');
    end
    legend([p1,p2,p3,p4], 'Positive Test', 'Negative Test', 'Alarms for positive test', 'Alarms for negative test')
    xlabel('Step')
    ylabel('g_t')
    title('CUSUM test for Alpha - Angle of attack)');
    hold off;
    
    % Subplot 2: Plotting y
    subplot(2, 1, 2);
    plot(y(:, m), 'b')
    title(['Measurement: ' titles(m)])
    grid on

end

%% plot
titles = {'x(1)=x', 'x(2)=y', 'x(3)=z', 'x(4)=u', 'x(5)=v', 'x(6)=w', 'x(7)=phi', 'x(8)=theta', 'x(9)=psi', 'x(10)=v_wxE', 'x(11)=v_wyE', 'x(12)=v_wzE', 'x(13)=b_A_x', 'x(14)=b_A_y', 'x(15)=b_A_z', 'x(16)=b_p', 'x(17)=b_q', 'x(18)=b_r'};

figure;
for i = 1:size(y, 2)
    plot(1:size(y, 1), y(:, i), 'b')
    hold on;
    % Plot vertical lines for positive alarms
    for j = 1:length(k_alarm_pos{i})
        line([k_alarm_pos{i}(j), k_alarm_pos{i}(j)], [-1, 2], 'Color', 'r');
    end
    % Plot vertical lines for negative alarms
    for j = 1:length(k_alarm_neg{i})
        line([k_alarm_neg{i}(j), k_alarm_neg{i}(j)], [-1, 2], 'Color', 'g');
    end
    hold off;
    title('Alpha - Angle of attack')
    legend('Measurement', 'Positive Alarms')
    grid on
end
