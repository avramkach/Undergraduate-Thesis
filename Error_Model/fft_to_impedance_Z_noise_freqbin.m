clear all
close all

%10 cycles is the same as buffer if statement

%% FFT to Impedance 10 cycles/periods for each test frequency (hacked without second fft)
%Run through normal frequencies just plot nyquist for now
close all

%test frequencies
w = [1952;1708;1464;1220;976;732;488;244;122;61;30.50000;15.25000;7.62500;3.81250;1.90625;0.95300;0.47600;0.40000;0.30000;0.20000;0.10000]*2*pi';
%w = [244;100]*2*pi;
n = length(w);
f = w/(2*pi);
%make last frequency 1.90625 hz? too long to simulate

%transfer function R + C circuit
syms wi

%ECM parameters from extraction paper without SEI
% L = 0.175*10^(-6);
% R0 = 0.550;
% R1 = 0.119;
% C1 = 1.146;
% sig = 0.0346;

%Cell5_0 extraction
L = 0;
R0 = 0.26424;
R1 = 0.21616;
C1 = 0.42183;
sig = 0.11111;

% L = 0.00005
% R0 = 0.5;
% R1 = 1;
% C1 = 1;

Zc = 1./(j*wi*C1);
ZL = j*wi*L;
Zw = sig./(sqrt(wi))*(1-j);

%R + R//C Circuit
%Tau = 1/(RC)
%Z = R0 + 1/(1/R1+1/Zc);

%L + R + R//C Circuit
%Z = R0 + 1/(1/R1+1/Zc) + ZL;

Z = R0 + 1/(1/(R1+Zw)+1/Zc);

%Z = 1 %1ohm

Zf = subs(Z, w); %output impedance at test frequency wi

%Change to cell for fft computation
Zf = sym2cell(Zf);
for i = 1:n
    Zo(i) = double(Zf{i});
end

%Sampling parameters
showp = 10; %periods to run
fs = 250000; %sampling frequency (samples/s)
Ts = 1/fs; %sample increments 
%10 cycles Ts

disp('FFT to Impedance 10 cycles/periods')
for i = 1:n
    %input signal I(t)
    fi = w(i)/(2*pi)
    
    t1 = Ts:Ts:showp*2*pi/w(i); %make it a fast fft for 2^N points
    x = sin(w(i)*t1);
    N = length(t1);
    
    %FFT input signal
    I = fft(x);
    %extract complex
    
    IP1 = 0;
    IP2 = 0;
    IP2 = I/N; %Take P2 = y0/L where y0 is the fft output
    IP1 = IP2(1:N/2+1);
    IP1(2:end-1) = 2*IP1(2:end-1);
    fbinsI = fs*(0:(N/2))/N; %fs is the set samping rate, fbins are the frequency bins
    [zeroI, nI] = min(abs(fbinsI - fi)); %Find closest frequency index to known input frequency
    
    Icom(i) = I(nI); 
    
    %Compute output signal
    %V(f) = I(f)*Z(f)
    %V = I*Zo(i);
    Vcom(i) = Icom(i)*Zo(i); %hack for computation
     
end


Zouth = Vcom./Icom;
figure(1);
plot(real(Zouth),-imag(Zouth),'o-'); %,'MarkerSize',14);
hold on

% FFT to Impedance with Buffer If Statement 2^N (1024) points makes it much
% faster
disp('FFT to Impedance with Buffer')
SNR = 20;
for i = 1:n
    %input signal I(t)
    
    fi = w(i)/(2*pi)

    if (fi > 244)
        %No downsampling, frequencies are always integer*244
        %Periods run for are fi/244 (i.e. 1952/244 = 8, runs for periods/cycles)

        fs = 250000; %%sampling rate (samples/s) (Can't set higher)
        Ts = 1/fs;

        %spacing is (x2-x1)/(n-1) = Ts.
        %n = 1024, x1 = 0
        %starting time x1 and ending time x2 set on 1 cycles using Ts b
        x1 = 0;
        x2 = Ts*(1024-1) + x1; %just set to rate/sample to fill in buffer
        t1 = linspace(x1, x2, 1024);
        N = length(t1);
        
    elseif (fi < 244 && fi > 1) %run downsampling with 3 cycles (or 2)
        %Frequencies are ~/2
        %Run 1 cycle with reduced fsample
        %Implement specific set sampling? Yes
        %Implement Decimation for higher frequency range? Yes do another if statement  
        cycles = 1;

        fs = fi/cycles * 1024; %Set reduced sampling, sampling constraint: Integer value they usually round down need to comfirm
        %Ts = 1/fs; ideal 
        
        nd = round(250000/fs);
        
        fs = 250000/nd;
        Ts = 1/fs;
        %fADC = floor(fs/250000)*1000;
        %fADC
        %spacing is (x2-x1)/(n-1) = Ts.
        %n = 1024, x1 = 0
        %starting time x1 and ending time x2 set on 1 cycles using Ts b
        x1 = 0;
        x2 = Ts*(1024-1) + x1;
        t1 = linspace(x1, x2, 1024);

        N = length(t1);

    elseif (fi < 1) %run downsampling with 1
        %Frequencies are ~/2
        %Run 1 cycle with reduced fsample
        %Implement specific set sampling? Yes
        %Implement Decimation for higher frequency range? Yes do another if statement  
        fs = floor(fi * 1024); %Set reduced sampling, sampling constraint: Integer value they usually round down need to comfirm
        %Ts = 1/fs;
        
        nd = round(250000/fs);
        
        fs = 250000/nd;
        Ts = 1/fs;

        %spacing is (x2-x1)/(n-1) = Ts.
        %n = 1024, x1 = 0
        %starting time x1 and ending time x2 set on 1 cycles using Ts b
        x1 = 0;
        x2 = Ts*(1024-1) + x1;
        t1 = linspace(x1, x2, 1024);

        N = length(t1);
        disp('<10Hz input')

    else %fi == 244Hz
        fs = 250000; %sampling rate (samples/s)
        Ts = 1/fs; %sampling increment (s/samples)

        %spacing is (x2-x1)/(n-1) = Ts.
        %n = 1024, x1 = 0
        %starting time x1 and ending time x2 set on 1 cycles using Ts b
        x1 = 0;
        x2 = Ts*(1024-1) + x1;
        t1 = linspace(x1, x2, 1024);
        N = length(t1);
    end
    
    
    x = sin(w(i)*t1);
    
    %SNR = SNR + 10; starting at 100 works make custom array
    %xn = awgn(x,SNR,'measured');
    xn = x; %no noise
    
    %FFT input signal
    %I = fft(x);
    
    %FFT input signal with noise
    In = fft(xn);
    %extract complex
    
    IP1 = 0;
    IP2 = 0;
    IP2 = In/N; %Take P2 = y0/L where y0 is the fft output
    IP1 = IP2(1:N/2+1);
    IP1(2:end-1) = 2*IP1(2:end-1);
    fbinsI = fs*(0:(N/2))/N; %fs is the set samping rate, fbins are the frequency bins
    [zeroI, nI] = min(abs(fbinsI - fi)); %Find closest frequency index to known input frequency
    %Single sided spectrum is fbins vs P1
    
    %Icom(i) = max(In); 
    
    %Compute output signal
    %V(f) = I(f)*Z(f)
    %V = I*Zo(i);
    
    %Output signal V(t) from V(f)
%     V = sym2cell(V);
%     V = V{1:end}
%     
%     y = ifft(V); %same number of samples and frequencies
%     
%     %Take V(f) from V(t)
%     V = fft(y);
    %y = ifft(V);
    
    %y = awgn(y,SNR,'measured');
    
    
    y = x*Zo(i); %V(t) = I(t)Z(fi) ?? works
    
    %y = awgn(y,SNR,'measured');
    
    V = fft(y);
    
    VP1 = 0;
    VP2 = 0;
    VP2 = V/N; %Take P2 = y0/L where y0 is the fft output
    VP1 = VP2(1:N/2+1);
    VP1(2:end-1) = 2*VP1(2:end-1);
    fbinsV = fs*(0:(N/2))/N; %fs is the set samping rate, fbins are the frequency bins
    [zeroV, nV] = min(abs(fbinsV - fi)); %Find closest frequency index to known input frequency
    %Single sided spectrum is fbins vs P1
    
    %Z in loop
    Zout(i) = VP1(nV)/IP1(nI);
end

plot(real(Zout),-imag(Zout),'x-','MarkerSize',14);
xlabel('Z_{real} (\Omega)');
ylabel('-Z_{imag} (\Omega)');
title('Nyquist Plot using NLLS Fitted Submodule 5 Parameters');
legend('10 cycles', 'System Buffer');

%Set axis limits
%xlim([0.5 inf]);
%ylim([0 inf]);

mag_Zouth = abs(Zouth);
phase_Zouth = angle(Zouth);

mag_Zout = abs(Zout);
phase_Zout = angle(Zout);

figure(4);
sgtitle('Phase and Mag Bode Plot Submodule 5 _ 0')
subplot(312);

r1 = plot(f,mag_Zouth,'o-');
hold on
r2 = plot(f,mag_Zout,'x-');
hold on

set(gca,'xscale','log');
ylabel('|Z| (\Omega)');

subplot(313);
i1 = plot(f,phase_Zouth,'o-');
hold on
i2 = plot(f,phase_Zout,'x-');
hold on

set(gca,'xscale','log');
xlabel('Frequency (Hz)');
ylabel('\phi(Z) (rad)');

hSub = subplot(311); plot(1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');
legend(hSub, [i1, i2], {'10 cycles', 'System Buffer with Noise'}, 'Location', 'South');
