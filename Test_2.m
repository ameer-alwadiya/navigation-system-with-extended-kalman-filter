% clear; clc; close all;

% x2 = load('KF Result_Task3_2(innovation).mat');
% x3 = load('KF Result_Task3_4(innovation).mat');
% 
x2 = Innov_nominal;
x3 = Innov_attack; 

 X2(:,1)=x2(:,11);
X3(:,1)=x3(:,11);
leak=1;

    straingauge=X3(500:end);
    theta_0 =0;%mean(X2(100:500,num));
    sigma_0 =std(straingauge);
  
     thres=theta_0+2;

    threshold_pos = thres;
    threshold_neg =-thres;
    n_alarm_pos=0; %number of alarms for postive test
    n_alarm_neg=0; %number of alarms for negative test
    k_alarm_pos=[]; %time indices k for alarms of postive test
    k_alarm_neg=[]; %time indices k for alarms of negative test
    % Two-Sided CUSUM Test with threshold = 20
    g_pos = 0*straingauge;
    g_neg = 0*straingauge;
    s = (straingauge-theta_0)/sigma_0;
    for k = 1:size(straingauge,1)-1
     g_pos(k+1) = g_pos(k)+s(k)-leak; %neglect leakage term
     g_neg(k+1) = g_neg(k)+s(k)+leak; %neglect leakage term
     %positive test
     if g_pos(k+1) <0
         g_pos(k+1)=0;
     end
     if g_pos(k+1) > threshold_pos 
        n_alarm_pos=n_alarm_pos+1;
        k_alarm_pos=[k_alarm_pos k+1];
        g_pos(k+1)=0; %reset
     end
     %negative test
     if g_neg(k+1) >0
         g_neg(k+1)=0;
     end
     if g_neg(k+1) < threshold_neg
        n_alarm_neg=n_alarm_neg+1;
        k_alarm_neg=[k_alarm_neg k+1];
        g_neg(k+1)=0;%reset
     end
    end
  
    figure
    p1=plot(g_pos);
    hold on
    p2=plot(g_neg);
    for i=1:length(k_alarm_pos)
        p3=plot([k_alarm_pos(i) k_alarm_pos(i)],[-20 20],'b--');
    end
    for i=1:length(k_alarm_neg)
        p4=plot([k_alarm_neg(i) k_alarm_neg(i)],[-20 20],'r-.');
    end
    %legend([p1,p2,p3,p4],'Positive Test', 'Negative Test','Alarms for positive test','Alarms for negative test')
    yline(threshold_neg)
    yline(threshold_pos)
    ylim([threshold_neg-0.025 threshold_pos+0.025])
    xlabel('Step')
    ylabel('g_t')
