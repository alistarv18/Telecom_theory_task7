clc, clear variables
pkg load signal;

hChan = [1 0.25 0.25];
L = length(hChan);

L_FF = 21;
L_FB = 3;

SNR = 10;

N_delay = 0;



% Generate information symbols, PAM-4
N = 1e4;
M = 4; N_train=2e3;

a = -2*(randi(M,1,N)-1)+M-1;
x_chan = filter(hChan,1,a);

Es = mean(abs(x_chan).^2);                  % Average Energy of transmitted symbol
Eb = Es/log2(M);                                        % Bit energy
EbN0 = 10.^(SNR/10);                                      % Convert dB -> Linear scale
N0 = Eb./EbN0;
sNoi =  sqrt(N0/2)*randn(1,N);

xRx = x_chan + sNoi;

[y,y_train, coeHist] = DFEequ(xRx,a,L_FF,L_FB, N_train);
y = y(1+N_delay:end);
a = a(1:N-N_delay);

##y = y_train(1+N_delay:end);
##a = a(1:N_train-N_delay);

figure(1)
subplot(211),plot(coeHist), grid on
xlabel('Sample no');ylabel('Filter coef');

% Calculate Cross-Correlation
[val, lag] = xcorr(y,a,'unbiased');
subplot(212),stem(lag,val), grid on
xlabel('Time Shift Tag'),ylabel('Cross-correlation');

NumErr = sum(y~=a)
ErrVec = y~=a
n = 1:length(ErrVec);n(ErrVec)
