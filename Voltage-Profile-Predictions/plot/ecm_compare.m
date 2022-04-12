close all

figure(1)

l_sim = length(vterm_R);

%time=[0:0.1:1800];
%14000 = stop time from simulink
%time=[0:T_est*100:13000];
time=[0:T_est*100:sim_time];

%time=[1:size(vterm_2RC, 1)];

%Peaks around 2415 index at 1208

subplot(311);
%plot(time(250:1000),drive_profile(250:1000));
plot(drive_profile_index(250:1000),drive_profile(250:1000));
ylabel('I_{bat} (A)');
%xlim([25 100]);

subplot(312);
hold on;
%plot(time(250:1000),vterm_2RC(250:1000));
%plot(time(250:1000),vterm_R(250:1000));

% plot(time(1:18001),vterm_2RC(1:18001));
% plot(time(1:18001),vterm_R(1:18001));
% plot(time(1:18001),exp_vterm(1:18001));

plot(time,vterm_2RC);
plot(time,vterm_R);

%plot((exp_vterm_index(1:3001)-1)/10,exp_vterm(1:3001));
plot(exp_vterm_index,exp_vterm);

ylabel('V_{term} (V)');
legend('2RC','R_{total}', 'exp vterm');
%xlim([25 100]);

subplot(313);
plot(time(250:1000),100*(vterm_2RC(250:1000)-vterm_R(250:1000))./vterm_2RC(250:1000));
ylabel('R vs RC Error (%)');
xlabel('Time (s)');
% xlim([25 100]);

%Analysis, Errors only work if same sized araays
%find max index of vterm, compare discharge and charge

[maxv, maxi] = max(exp_vterm);

maxv_time = exp_vterm_index(maxi);

%vq1 = interp1(x,v,xq); %default method is 'linear'
%x: known x values (voltages)
%v: known y vales y = v(x) (SOCs)
%xq: query points (interpolated x point with data)

maxv_simtime_index = interp1(time, 1:length(time),maxv_time); %index at 2417s


figure(2)
subplot(211);
hold on
plot(time(1:maxv_simtime_index), vterm_2RC(1:maxv_simtime_index))
plot(time(1:maxv_simtime_index), vterm_R(1:maxv_simtime_index))
plot(exp_vterm_index(1:maxi), exp_vterm(1:maxi))
ylabel('V_{term} (V)');
legend('2RC','R_{total}', 'exp vterm', 'Location', 'Best');

subplot(212);
hold on
plot(time(1:maxi), 100*(vterm_2RC(1:maxi) - exp_vterm(1:maxi))./exp_vterm(1:maxi))
plot(time(1:maxi), 100*(vterm_R(1:maxi) - exp_vterm(1:maxi))./exp_vterm(1:maxi))
ylabel('V_{term} Error (%)');
legend('2RC','R_{total}', 'Location', 'Best');
xlabel('Time (s)');
sgtitle('Charging Profile')

figure(3)
subplot(211);
hold on
plot(time(maxv_simtime_index+1:end), vterm_2RC(maxv_simtime_index+1:end))
plot(time(maxv_simtime_index+1:end), vterm_R(maxv_simtime_index+1:end))
plot(exp_vterm_index(maxi+1:end), exp_vterm(maxi+1:end))
ylabel('V_{term} (V)');
legend('2RC','R_{total}', 'exp vterm', 'Location', 'Best');

subplot(212);
hold on
plot(time(maxi+1:end), 100*(vterm_2RC(maxi+1:end) - exp_vterm(maxi+1:l_sim))./exp_vterm(maxi+1:l_sim))
plot(time(maxi+1:end), 100*(vterm_R(maxi+1:end) - exp_vterm(maxi+1:l_sim))./exp_vterm(maxi+1:l_sim))
ylabel('V_{term} Error (%)');
legend('2RC','R_{total}', 'Location', 'Best');
xlabel('Time (s)');
sgtitle('Discharging Profile')


RMSE_vterm_2RC = sqrt(mean((exp_vterm(1:l_sim) - vterm_2RC).^2));
RMSE_vterm_R = sqrt(mean((exp_vterm(1:l_sim) - vterm_R).^2));
disp('RMSE vterm_2RC')
disp(RMSE_vterm_2RC)
disp('RMSE vterm_R')
disp(RMSE_vterm_R)
