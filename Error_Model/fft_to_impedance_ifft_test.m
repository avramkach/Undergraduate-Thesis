clear all
close all

%10 cycles is the same as buffer if statement

%% FFT to Impedance 10 cycles/periods for each test frequency (hacked without second fft)
%Run through normal frequencies just plost nyquist for now

%test frequencies
%w = [1952;1708;1464;1220;976;732;488;244;122;61;30.50000;15.25000;7.62500;3.81250;1.90625;0.95300;0.47600;0.40000;0.30000;0.20000;0.10000]*2*pi';

w = [15.25000;7.62500]*2*pi';
n = length(w);



%transfer function R + C circuit
syms wi

%ECM parameters from extraction paper without SEI
L = 0.175*10^(-6);
R0 = 0.550;
R1 = 0.119;
C1 = 1.146;
sig = 0.0346;

% R0 = 0.5;
% R1 = 1;
% C1 = 1;

Zc = 1./(j*wi*C1);
ZL = j*wi*L;
Zw = sig./(sqrt(wi))*(1-j);

%R + R//C Circuit
%Tau = 1/(RC)
%Z = R0 + 1/(1/R1+1/Zc)
%Z = R0 + 1/(1/R1+1/Zc) + Zw;
Z = R0 + 1/(1/R1+1/Zc) + ZL + Zw;

Zf = subs(Z, w); %output impedance at test frequency wi

%Change to cell for fft computation
Zf = sym2cell(Zf);
for i = 1:n
    Zo(i) = double(Zf{i});
end

%Sampling parameters
showp = 2; %periods to run
fs = 250000; %sampling frequency (samples/s)
Ts = 1/fs; %sample increments 

disp('FFT to Impedance 2 cycles/periods')
for i = 1:n
    %input signal I(t)
    fi = w(i)/(2*pi)
    
    if fi > 10
        cycles = 3;
    else 
        cycles = 1;
    end
    
    %force cycles 3
    cycles = 3;
    
    fs = floor(fi/cycles * 1024) %Set reduced sampling, sampling constraint: Integer value they usually round down need to comfirm
    Ts = 1/fs;

    %spacing is (x2-x1)/(n-1) = Ts.
    %n = 1024, x1 = 0
    %starting time x1 and ending time x2 set on 1 cycles using Ts b
    x1 = 0;
    x2 = Ts*(1024-1) + x1;
    t1 = linspace(x1, x2, 1024);
    
    %t1 = Ts:Ts:showp*2*pi/w(i); %make it a fast fft for 2^N points
    x = sin(w(i)*t1);
    
    xt = x;
    t1t = t1;
    
    %FFT input signal
    I = fft(x);
    %extract complex
    Icom(i) = max(I); 
    
    %Compute output signal
    %V(f) = I(f)*Z(f)
    %V = I*Zo(i);
    V = Icom(i)*Zo(i); %hack for computation
    
    %extract complex
    Vcom(i) = max(V);
    
    V = I*Zo(i);
    y = ifft(V);
    V = fft(y);
    
    %Z in loop
    Z = V./I;
    
    [maxz, nz] = max(abs(Z));
    
    Zouth(i) = Z(nz);
    
    N = length(t1);
    f = (0:N-1)*fs/N; %0, fs/N, 2*fs/N,...,(N-1)*fs/N
    
    figure(1);
    sgtitle('3 Cycles for Each Test Frequency')
    hSub = subplot(411); plot(1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');
    lgd = legend(hSub, string(w/(2*pi)), 'Location', 'south');
    title(lgd, 'Frequency (Hz)');
    
    subplot(412);
    plot(t1,x)
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    title('Output Measurement Signal V(t) = sin(wt)')
    hold on

    subplot(413);
    plot(f,real(Z));
    xlabel('Frequency (Hz)');
    title('Z Real Component from FFT');
    hold on
    
    subplot(414);
    plot(f,imag(Z));
    xlabel('Frequency (Hz)');
    title('Z Imaginary Component from FFT');
    hold on
    
    figure(5);
    sgtitle('V and I: 3 Cycles for Each Test Frequency')
    hSub = subplot(511); plot(1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');
    lgd = legend(hSub, string(w/(2*pi)), 'Location', 'south');
    title(lgd, 'Frequency (Hz)');
    
    subplot(512);
    plot(f,real(I));
    xlabel('Frequency (Hz)');
    title('I Real Component from FFT');
    hold on
    
    subplot(513);
    plot(f,imag(I));
    xlabel('Frequency (Hz)');
    title('I Real Component from FFT');
    hold on
    
    subplot(514);
    plot(f,real(V));
    xlabel('Frequency (Hz)');
    title('V Real Component from FFT');
    hold on
    
    subplot(515);
    plot(f,imag(V));
    xlabel('Frequency (Hz)');
    title('V Real Component from FFT');
    hold on
end

Zout = Vcom./Icom;
figure(3);
plot(real(Zout),-imag(Zout),'.-','MarkerSize',14);
hold on
plot(real(Zouth),-imag(Zouth),'o-','MarkerSize',14);
hold on

% FFT to Impedance with Buffer If Statement 2^N (1024) points makes it much
% faster
disp('FFT to Impedance with Buffer')

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
        
    elseif (fi < 244 && fi > 10) %run downsampling with 3 cycles (or 2)
        %Frequencies are ~/2
        %Run 1 cycle with reduced fsample
        %Implement specific set sampling? Yes
        %Implement Decimation for higher frequency range? Yes do another if statement  
        cycles = 3;

        fs = floor(fi/cycles * 1024) %Set reduced sampling, sampling constraint: Integer value they usually round down need to comfirm
        Ts = 1/fs;

        %spacing is (x2-x1)/(n-1) = Ts.
        %n = 1024, x1 = 0
        %starting time x1 and ending time x2 set on 1 cycles using Ts b
        x1 = 0;
        x2 = Ts*(1024-1) + x1;
        t1 = linspace(x1, x2, 1024);

        N = length(t1);

    elseif (fi < 10) %run downsampling with 1
        %Frequencies are ~/2
        %Run 1 cycle with reduced fsample
        %Implement specific set sampling? Yes
        %Implement Decimation for higher frequency range? Yes do another if statement  
        fs = floor(fi * 1024) %Set reduced sampling, sampling constraint: Integer value they usually round down need to comfirm
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
        fi = w(i)/(2*pi);
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
    
    
    %FFT input signal
    I = fft(x);
    
    %extract complex
    Icom(i) = max(I); 
    
    %Compute output signal
    %V(f) = I(f)*Z(f)
    V = I*Zo(i);
    
    y = ifft(V);
    
    V = fft(y);
    
    %Z in loop
    Z = V./I;
    
    [maxz, nz] = max(abs(Z));
    
    Zout(i) = Z(nz);
    
    N = length(t1);
    f = (0:N-1)*fs/N; %0, fs/N, 2*fs/N,...,(N-1)*fs/N
    
    figure(2);
    sgtitle('Buffer: 3 Cycles for 15Hz and 1 cycle for 7Hz')
    hSub = subplot(411); plot(1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');
    lgd = legend(hSub, string(w/(2*pi)), 'Location', 'south');
    title(lgd, 'Frequency (Hz)');
    
    subplot(412);
    plot(t1,x)
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    title('Output Measurement Signal V(t) = sin(wt)')
    hold on

    subplot(413);
    plot(f,real(Z));
    xlabel('Frequency (Hz)');
    title('Z Real Component from FFT');
    hold on
    
    subplot(414);
    plot(f,imag(Z));
    xlabel('Frequency (Hz)');
    title('Z Imaginary Component from FFT');
    hold on
    
    figure(4);
    sgtitle('V and I Buffer: 3 Cycles for 15Hz and 1 cycle for 7Hz')
    hSub = subplot(511); plot(1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');
    lgd = legend(hSub, string(w/(2*pi)), 'Location', 'south');
    title(lgd, 'Frequency (Hz)');
    
    subplot(512);
    plot(f,real(I));
    xlabel('Frequency (Hz)');
    title('I Real Component from FFT');
    hold on
    
    subplot(513);
    plot(f,imag(I));
    xlabel('Frequency (Hz)');
    title('I Real Component from FFT');
    hold on
    
    subplot(514);
    plot(f,real(V));
    xlabel('Frequency (Hz)');
    title('V Real Component from FFT');
    hold on
    
    subplot(515);
    plot(f,imag(V));
    xlabel('Frequency (Hz)');
    title('V Real Component from FFT');
    hold on
end

figure(3);
plot(real(Zout),-imag(Zout),'x-','MarkerSize',14);
xlabel('Z_{real} (\Omega)');
ylabel('-Z_{imag} (\Omega)');
title('Nyquist');
legend('3 Cycles Correct', '3 Cycles IFFT', 'System Buffer');

%Set axis limits
%xlim([0.5 inf]);
%ylim([0 inf]);