%% EN2570 Digital Signal Processing Project 
clc; 
clear all; 
close all;

%% Calculate the paramteres

indexNo = str2double('180039');
%indexNo = str2double(inputdlg('Enter Index Number','Index Number',1));

%Get A,B,C in 180ABC ( index )
A = mod(floor(indexNo/100),10);
B = mod(floor(indexNo/10),10);
C = mod(indexNo,10);
%A = 0; B = 3; C = 9;

A_tilda_p = 0.03 + (0.01 * A) ; %dB Maximum passband ripple, A˜p
A_tilda_a = 45 + B;             %dB Minimum stopband attenuation, A˜a
Omega_p1 = C*100 + 300;         %rad/s Lower passband edge
Omega_p2 = C*100 + 700;         %rad/s Upper passband edge 
Omega_a1 = C*100 + 150;         %rad/s Lower stopband edge
Omega_a2 = C*100 + 800;         %rad/s Upper stopband edge
Omega_s = 2*(C*100 + 1200);     %rad/s Sampling frequency

%% Derived Filter Specifications

B_t1 = Omega_p1-Omega_a1;       %rad/s lower transition width
B_t2 = Omega_a2-Omega_p2;       %rad/s upper transisiton width
B_t = min(B_t1,B_t2);           %rad/s critical transition width
Omega_c1 = Omega_p1-B_t/2;      %rad/s lower cutoff frequency
Omega_c2 = Omega_p2+B_t/2;      %rad/s upper cutoff frequency
T = 2*pi/Omega_s;               %s sampling period

%% Kaiser Window Parameters

%Choose Delta
delta_p = ((10^(0.05*A_tilda_p))-1)/((10^(0.05*A_tilda_p))+1); 
delta_a = 10^(-1*(0.05*A_tilda_a)); 
delta = min(delta_p,delta_a);

%Actual Stopband attenuation 
A_a = -20*log10(delta);

%Choose Alpha
if A_a<=21
    alpha = 0; 
elseif A_a>21 && A_a<= 50 
    alpha = 0.5842*(A_a-21)^0.4 + 0.07886*(A_a-21); 
else
    alpha = 0.1102*(A_a-8.7); 
end

%Choose D
if A_a <= 21
    D = 0.9222;
else
    D = (A_a-7.95)/14.36;
end

% Calculating order of the filter N
N = ceil((Omega_s*D/B_t) +1);
if mod(N,2) == 0
    N = N+1;
end

%Range of N
n_lim = floor(N-1)/2;

% Length of the filter
n = -n_lim:1:n_lim;

% Calculating beta
beta = alpha*sqrt(1-(2*n/(N-1)).^2);

%% Generating I_0_alpha 

bessellimit =50; 
I_0_alpha = 1; 
for k = 1:bessellimit 
    termk = (1/factorial(k)*(alpha/2).^k).^2; 
    I_0_alpha = I_0_alpha + termk; 
end

%% Generating I_0_beta
I_0_beta = 1; 
for k = 1:bessellimit 
    termk = (1/factorial(k)*(beta/2).^k).^2; 
    I_0_beta = I_0_beta + termk; 
end

%% Obtaining Kaiser Window w_k(nT) 
wk_nt = I_0_beta/I_0_alpha; 
figure 
stem(n,wk_nt)
xlim([-70 70])
xlabel('n') 
ylabel('Amplitude') 
title('Kaiser Window - Time Domain');

%% Generating Impulse Response h(nT)
n_left = -(N-1)/2:-1;
hnt_left = 1./(n_left*pi).*(sin(Omega_c2*n_left*T)-sin(Omega_c1*n_left*T));

n_right = 1:(N-1)/2; 
hnt_right = 1./(n_right*pi).*(sin(Omega_c2*n_right*T)-sin(Omega_c1*n_right*T));

hnt_0 = 2/Omega_s*(Omega_c2-Omega_c1);

hnt = [hnt_left,hnt_0,hnt_right]; 
figure 
stem(n,hnt)
xlim([-70 70])
xlabel('n') 
ylabel('Amplitude') 
title(strcat(['Filter Response - Rectangular window - Time Domain']));

figure 
[h,w] = freqz(hnt); 
w = w/T; 
h = 20*log10(abs(h)); 
plot(w,(h))
xlim([0 2100])
xlabel('Frequency (rad/s)') 
ylabel('Magnitude (dB)') 
title(strcat(['Filter Response - Rectangular Window - Frequency Domain']));

%% Applying the Kaiser Window to the filter

Hw_nT = hnt.*wk_nt;

figure
stem(n,Hw_nT);
xlim([-70 70])
xlabel('n')
ylabel('Amplitude')
title(strcat(['Filter Response - Kaiser Window - Time Domain']));
%fvtool(Hw_nT);

figure 
[h,w] = freqz(Hw_nT); 
w = w/T; 
h = 20*log10(abs(h)); 
plot(w,h)
xlim([0 2100])
xlabel('Frequency (rad/s)') 
ylabel('Magnitude (dB)') 
title(strcat(['Filter Response - Kaiser Window - Frequency Domain']));

%% Plotting the Pass band

figure 
start = round(length(w)/(Omega_s/2)*Omega_c1); 
finish = round(length(w)/(Omega_s/2)*Omega_c2); 
w_pass = w(start:finish); 
h_pass = h(start:finish); 
plot(w_pass,h_pass) 
axis([-inf, inf, -0.1, 0.1]); 
xlabel('Frequency (rad/s)') 
ylabel('Magnitude (dB)') 
title('Passband - Frequency Domain');

%% Input signal generation

% component frequencies of the input
Omega_1 = Omega_c1/2;
Omega_2 = Omega_c1 + (Omega_c2-Omega_c1)/2;
Omega_3 = Omega_c2 + (Omega_s/2-Omega_c2)/2;

% generate discrete signal and evelope
samples = 500;
n1 = 0:1:samples;
n2 = 0:0.1:samples;
X_nT = sin(Omega_1.*n1.*T)+sin(Omega_2.*n1.*T)+sin(Omega_3.*n1.*T);
X_env = sin(Omega_1.*n2.*T)+sin(Omega_2.*n2.*T)+sin(Omega_3.*n2.*T);

%% Using DFT to check the filtering

% Filtering using frequency domain multiplication
len_fft = length(X_nT)+length(Hw_nT)-1; % length for fft in x dimension
x_fft = fft(X_nT,len_fft);
Hw_nT_fft = fft(Hw_nT,len_fft);
out_fft = Hw_nT_fft.*x_fft; % A shift in time is added here
out = ifft(out_fft,len_fft);
rec_out = out(floor(N/2)+1:length(out)-floor(N/2)); % account for shifting delay

% Ideal Output Signal
ideal_out = sin(Omega_1.*n2.*T)+sin(Omega_3.*n2.*T);

% Frequency domain representation of input signal before filtering 
figure
subplot(2,1,1)
len_fft = 2^nextpow2(numel(n1))-1;
x_fft = fft(X_nT,len_fft);
x_fft_plot = [abs([x_fft(floor(len_fft/2)+1:len_fft)]),abs(x_fft(1)),abs(x_fft(2:floor(len_fft/2)))];
f = Omega_s*linspace(0,1,len_fft)-Omega_s/2;
plot(f,x_fft_plot); 
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Input signal',' ','- Frequency Domain']));
% Time domain representation of input signal before filtering
subplot(2,1,2)
stem(n1,X_nT)
%xlim([200 250])  %Uncomment for zoomed plot
xlabel('n')
ylabel('Amplitude')
title(strcat(['Input signal',' ','- Time Domain']));
% Frequency domain representation of output signal after filtering
figure
subplot(3,1,1)
len_fft = 2^nextpow2(numel(n1))-1;
xfft_out = fft(rec_out,len_fft);
x_fft_out_plot =[abs([xfft_out(floor(len_fft/2)+1:len_fft)]),abs(xfft_out(1)),abs(xfft_out(2:floor(len_fft/2)))];
f = Omega_s*linspace(0,1,len_fft)-Omega_s/2;
plot(f,x_fft_out_plot); 
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Output signal',' ','- Frequency Domain']));
% Time domain representation of output signal after filtering
subplot(3,1,2)
stem(n1,out(1:samples+1))
%xlim([200 250])    %Uncomment for zoomed plot
xlabel('n')
ylabel('Amplitude')
title(strcat(['Output signal from designed filter',' ','- Time Domain']));
subplot(3,1,3)
stem(n1,rec_out)
%xlim([200 250])    %Uncomment for zoomed plot
xlabel('n')
ylabel('Amplitude')
title(strcat(['Expected Output signal',' ','- Time Domain']));
