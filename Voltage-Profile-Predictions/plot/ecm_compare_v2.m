close all

%Plotting scripting: Run load_sim_data.m then battery.slx
%Check indices, times, and coloumb count to ensure functionalit
%Change T_est in load_sim_data.m to match testing time step /100

l_sim = length(vterm_R);

%sim_time = stop time from simulink
%steps = T_est*100 (decimation)
time=[0:T_est*100:sim_time];

figure(1)
subplot(311)
plot(drive_profile_index(250:1000),drive_profile(250:1000));
ylabel('I_{bat} (A)');
%xlim([25 100]);

subplot(312);
hold on;
plot(time,vterm_2RC);
plot(time,vterm_R);
plot(exp_vterm_index,exp_vterm);
ylabel('V_{term} (V)');
legend('2RC','R_{0}', 'exp vterm');

subplot(313);
plot(time(250:1000),100*(vterm_2RC(250:1000)-vterm_R(250:1000))./vterm_2RC(250:1000));
ylabel('R vs RC Error (%)');
xlabel('Time (s)');

%%
%Analysis, Errors only work if same sized araays
%find max index of vterm, compare discharge and charge

%Peaks around 2415 index at 1208
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
legend('2RC','R_{0}', 'Location', 'Best');
xlabel('Time (s)');
sgtitle('Charging Profile')

figure(3)
subplot(211);
hold on
plot(time(maxv_simtime_index+1:end), vterm_2RC(maxv_simtime_index+1:end))
plot(time(maxv_simtime_index+1:end), vterm_R(maxv_simtime_index+1:end))
plot(exp_vterm_index(maxi+1:end), exp_vterm(maxi+1:end))
ylabel('V_{term} (V)');
legend('2RC','R_{0}', 'exp vterm', 'Location', 'Best');

subplot(212);
hold on
plot(time(maxi+1:end), 100*(vterm_2RC(maxi+1:end) - exp_vterm(maxi+1:l_sim))./exp_vterm(maxi+1:l_sim))
plot(time(maxi+1:end), 100*(vterm_R(maxi+1:end) - exp_vterm(maxi+1:l_sim))./exp_vterm(maxi+1:l_sim))
ylabel('V_{term} Error (%)');
legend('2RC','R_{0}', 'Location', 'Best');
xlabel('Time (s)');
sgtitle('Discharging Profile')

RMSE_vterm_2RC = sqrt(mean((exp_vterm(1:l_sim) - vterm_2RC).^2));
RMSE_vterm_R = sqrt(mean((exp_vterm(1:l_sim) - vterm_R).^2));
disp('RMSE vterm_2RC')
disp(RMSE_vterm_2RC)
disp('RMSE vterm_R')
disp(RMSE_vterm_R)

%% Discharge error

figure(10);
plot(time(maxi+1:end), 100*(vterm_2RC(maxi+1:end) - exp_vterm(maxi+1:l_sim))./exp_vterm(maxi+1:l_sim))
ylabel('V_{term} Error (%)');
%legend('2RC','R_{total}', 'Location', 'Best');
xlabel('Time (s)');
sgtitle('2RC Error Circuit Discharging Profile')
ylim([-5, 5])

figure(11);
plot(time(maxi+1:end), 100*(vterm_R(maxi+1:end) - exp_vterm(maxi+1:l_sim))./exp_vterm(maxi+1:l_sim))
ylabel('V_{term} Error (%)');
%legend('2RC','R_{total}', 'Location', 'Best');
xlabel('Time (s)');
sgtitle('R_{0} Error Circuit Discharging Profile')
ylim([-5, 5])

RMSE_vterm_2RC_d = sqrt(mean((exp_vterm(maxi+1:end) - vterm_2RC(maxi+1:end)).^2));
RMSE_vterm_R_d = sqrt(mean((exp_vterm(maxi+1:end) - vterm_R(maxi+1:end)).^2));
disp('discharge')
disp('RMSE vterm_2RC')
disp(RMSE_vterm_2RC_d)
disp('RMSE vterm_R')
disp(RMSE_vterm_R_d)

%% 

figure(4);
hold on;
plot(time,vterm_2RC);
plot(time,vterm_R);
plot(exp_vterm_index,exp_vterm);
title("Model S Module 5 Voltage Profile");
ylabel('V_{term} (V)');
legend('2RC','R_{0}', 'xperimental vterm');
xlabel('Time (s)');

%% FFT vterm, Rtotal and 2RC

var_y = 0.005;

%vterm
T_vterm = 2;             % Sampling period  
Fs_vterm = 1/T_vterm;            % Sampling frequency                    
     
%already have time vector 
L_vterm = length(exp_vterm);   % Length of signal
%exp_vterm_index    % Time vector

Y_vterm = fft(exp_vterm);
P2_vterm = abs(Y_vterm/L_vterm);
P1_vterm = P2_vterm(1:L_vterm/2+1);

f_vterm = Fs_vterm*(0:(L_vterm/2))/L_vterm;
figure(5);
plot(f_vterm, P1_vterm)
title('FFT of Experimental Battery Voltage');
xlabel('f (Hz)')
ylabel('|V_{term}(f)| (V)')
ylim([0 var_y])

%Rtotal and 2RC
T_p = 100*T_est;             % Sampling period  
Fs_p = 1/T_p;            % Sampling frequency                    
     
%already have time vector 
L_p = length(time);   % Length of signal
%time    % Time vector

%Rtotal
Y_R = fft(vterm_R);
P2_R = abs(Y_R/L_p);
P1_R = P2_R(1:L_p/2+1);

f_p = Fs_p*(0:(L_p/2))/L_p;
figure(6);
plot(f_p, P1_R)
title('FFT of Predicted Battery Voltage R');
xlabel('f (Hz)')
ylabel('|V_{Rtotal}(f)| (V)')
ylim([0 var_y])
xlim([0 0.25])

%2RC
Y_2RC = fft(vterm_2RC);
P2_2RC = abs(Y_2RC/L_p);
P1_2RC = P2_2RC(1:L_p/2+1);

figure(7);
plot(f_p, P1_2RC)
title('FFT of Predicted Battery Voltage 2RC');
xlabel('f (Hz)')
ylabel('|V_{Rtotal}(f)| (V)')
ylim([0 var_y])
xlim([0 0.25])


%% 2RC FFT error for same sampling period
var_ye = 0.001;
figure(8);
plot(f_vterm, P1_vterm-P1_R)
title('FFT Error of Predicted Battery Voltage R');
xlabel('f (Hz)')
ylabel('|V(f)| (V)')
ylim([-var_ye var_ye])

figure(9);
plot(f_vterm, P1_vterm-P1_2RC)
title('FFT Error of Predicted Battery Voltage 2RC');
xlabel('f (Hz)')
ylabel('|V(f)| (V)')
ylim([-var_ye var_ye])

%% Submodule 5 results

%lecm
%Q = 236.8
%ext_v2
% RMSE vterm_2RC
%     0.1108
% 
% RMSE vterm_R
%     0.1109
%ext_c_v1
% RMSE vterm_2RC
%     0.1100
% 
% RMSE vterm_R
%     0.1100

%baseline
%ext_v2
% RMSE vterm_2RC
%     0.1326
% 
% RMSE vterm_R
%     0.1326
%ext_c_v1
% RMSE vterm_2RC
%     0.1317
% 
% RMSE vterm_R
%     0.1317

%lecm
%Q = 350
%ext_c_v1
% RMSE vterm_2RC
%     0.0841
% 
% RMSE vterm_R
%     0.0842

%baseline
%Q = 350
%ext_c_v1
% RMSE vterm_2RC
%     0.1023
% 
% RMSE vterm_R
%     0.1024

%baseline
%Q = 350
%ext_c_v2
% RMSE vterm_2RC
%     0.0905
% 
% RMSE vterm_R
%     0.0915
%adjusted vTerm = -0.1
% RMSE vterm_2RC
%     0.0238
% 
% RMSE vterm_R
%     0.0275
%ext_c_v1, vterm = -0.01
% RMSE vterm_2RC
%     0.1023
% 
% RMSE vterm_R
%     0.1024