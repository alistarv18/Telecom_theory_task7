clear
SNR = 10:0.5:20;

N = 1e4;                    % Number of transmitted QAM-16 symbols (recommend 1e4 for debugging)
M = 16;                     % Size of QAM-16 alphabet (each orthogonal component of QAM-16 is PAM-4)
mBits = log2(M);            % Size of QAM-16 symbol binary payload (4 bits)

% Impulse response of channel
hChan = [1 0.25 0.25];
L = length(hChan);                        % Length of channel's IR

% Viterbi algorithm related parameters
NdelayVit =  5*(L-1);               % Delay during traceback procedure before starting to read received data

% MSE equalizer related parameters
Leq = 45;	                    % Order of MSE equalizer (length of equalizer IR)
NdelayMSE = (Leq-1)/2;		        % The maximum delay of MSE equalizer (half of the equalizer IR length)

% DFE equalizer related parameters
% Add them as needed by custom made function
L_FF = 51;
L_FB = 1;
N_train = 5000;
NdelayDFE = 0;
% Generate QAM16 signal
% Specify QAM16 constellation map as a vector-column
QAM16map = [-1-1i; -1-3i; -1+1i; -1+3i;
            -3-1i; -3-3i; -3+1i; -3+3i;
            +1-1i; +1-3i; +1+1i; +1+3i;
            +3-1i; +3-3i; +3+1i; +3+3i];

% Generate QAM-16 symbols
symDt =  randi(M,1,N);                   % Indices of the QAM16 symbols

% QAM-16 transmitted symbols
a = QAM16map(symDt).';

% Channel output adding ISI
sChanOu = filter(hChan,1,a);

% Evaluate Noise spectral power Density
Es = mean(abs(sChanOu).^2);                  % Average Energy of transmitted symbol
Eb = Es/mBits;                                        % Bit energy
EbN0 = 10.^(SNR/10);                                      % Convert dB -> Linear scale
N0 = Eb./EbN0;                                        % Noise power spectral density

for k = 1:length(SNR)
  fprintf('SNR = %.1f dB\n\n',SNR(k));

  % Add noise to received signal
  sNoi =  sqrt(N0(k)/2)*(randn(1,N) + 1i*randn(1,N));                                             % Generated complex noise, noise power density is N0
  sRx  =  sChanOu + sNoi;                                             % Received signal = Channel Output signal + Noise


  % =========================== Viterbi Algorithm part =========================
  % Apply Viterbi algorithm for Real and Imag parts separately Use provided ViterbiDetectFun() function.
  % If You have access to MATLAB with Communications Toolbox You can use MLSEEQ() function here.
  fprintf('Real part: '); rcvDataRe = ViterbiDetectFun(hChan, real(sRx),NdelayVit);
  fprintf('Imag part: '); rcvDataIm =ViterbiDetectFun(hChan, imag(sRx),NdelayVit);
  fprintf('\n');

  % Form a complex received symbols sequence
  rcvData = rcvDataRe + rcvDataIm * 1i;

  % Demodulate, evaluate bit error rate
  binRx = QAM16demod(rcvData, QAM16map);                                                          % Use previous Practical Exercises QAM16demod() function here!
  binDt =  (dec2bin(symDt(1:end-NdelayVit-1)-1) == '1')';
  binDt = binDt(:); % Truncate last L symbols before conversion!
  BER_Viterbi(k) = sum(binRx ~= binDt)/length(binRx);                                                % Number of mismatches (bit errors)


  % =========================== MSE Equalizer part =============================
  % Estimate equalizer, separately Real and Imag parts.
  % You need to create Your own function for this - see Practical Exercise 5.
  fprintf('Processing MSE Equalizer...\n');
  fprintf('Real part: ');  hEqRe = MSEEequFun(real(sRx), real(a), Leq,hChan);
  fprintf('Imag part: ');  hEqIm = MSEEequFun(imag(sRx), imag(a), Leq,hChan);

  % Apply equalizer to each QAM-16 sub-channel
  QkRe = filter(hEqRe,1,real(sRx));
  QkIm = filter(hEqIm,1,imag(sRx));

  % Form a complex received symbols sequence
  Qk = QkRe + 1i*QkIm;

  % Apply delay of Leq/2 caused by equalizer's filter
  Qk = Qk(NdelayMSE+1:end);
  symk = symDt(1:end-NdelayMSE);

  % Demodulate, evaluate bit error rate
  binRx = QAM16demod(Qk,QAM16map);
  binDt = (dec2bin(symk-1) == '1')'; binDt = binDt(:);
  BER_MSE(k) = sum(binRx ~= binDt)/length(binRx);

  % =========================== DFE Equalizer part =============================
  % Estimate equalizer, separately Real and Imag parts.
  % You need to create Your own function for this - see PE07 Part 2 record.
  % Apply equalizer to each QAM-16 sub-channel
  fprintf('\nProcessing DFE Equalizer...\n');
  fprintf('Real part: '); dtRxRe = DFEequ(real(sRx),real(a),L_FF,L_FB, N_train);
  fprintf('Imag part: '); dtRxIm = DFEequ(imag(sRx),imag(a),L_FF,L_FB, N_train);


  % Form a complex received symbols sequence
  dtRx = dtRxRe + dtRxIm*1i;
  % Apply delay of Leq/2 caused by equalizer's filter
  dtRx = dtRx(1+NdelayDFE:end);
  dtTx = symDt(1:end-NdelayDFE);
  % Demodulate, evaluate bit error rate
    binRx = QAM16demod(dtRx,QAM16map);
  binDt = (dec2bin(dtTx-1) == '1')'; binDt = binDt(:);
  BER_DFE(k) = sum(binRx ~= binDt)/length(binRx);

  fprintf('\n------------------------------------------------------------\n');
end

% Calculate theoretical BER for QAM-16 (Gray Mapping assumed)
alfa = sqrt(3*EbN0*mBits/(M-1)) ;                 % Argument of Q function
Perr = 4/mBits * 0.5 * erfc(alfa/sqrt(2));        % Q function calculation. Note: mBits = log2(M)
                                                  % Q(a) = 1 - (0.5 + 0.5*erf(a/sqrt(2)))


% Plot results. See Provided reference in Practical Exercise guidelines
figure(1)
semilogy(SNR, BER_Viterbi,'-rx'), grid on;      % Viterbi results
hold on, semilogy(SNR, BER_MSE,'-bx'), hold off; % MSE results
hold on, semilogy(SNR, BER_DFE,'-gx'), hold off;           	      % DFE results
hold on, semilogy(SNR,Perr,'-*'), hold off;       % Theoretical results

% Add plot annotations (axis, title, limit)
xlabel('E_b/N_0'); ylabel('Bit Error Rate');
title('BER for QAM-16 with different equalizers')
ylim([1e-8 1e-1])
legend('Viteerbi','MMSE','DFE','Theory');
