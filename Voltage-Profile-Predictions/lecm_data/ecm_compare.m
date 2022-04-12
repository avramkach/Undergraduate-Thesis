figure;

%time=[0:0.1:1800];
time=[0:T_est*100:600]; %600 sec

subplot(311);
plot(time(250:1000),drive_profile(250:1000));
ylabel('I_bat (V)');
xlim([25 100]);


%need to make sure times align with indices of exp_vterm

subplot(312);
hold on;
plot(time,vterm_2RC);
plot(time,vterm_R);

%plot(time,exp_vterm(1:18001));

%times need to line up, dt = 200ms and starts at 100ms?, dt = 2s
plot((exp_vterm_index(1:3001)-1)/10,exp_vterm(1:3001)); %run for 600 sec and plot to match 6001 samples generated for vterms

ylabel('V_{term} (V)');
legend('2RC','R_{total}', 'exp');
%xlim([0 100]);

subplot(313);
plot(time(250:1000),100*(vterm_2RC(250:1000)-vterm_R(250:1000))./vterm_2RC(250:1000));
ylabel('Error (%)');
xlabel('Time (s)');
xlim([25 100]);

