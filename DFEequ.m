function [yd,yd_train, coeffsHist] = DFEequ(x,y_train,L_FF,L_FB, N_train)

  stepsize = 0.001;

  U = [+3; +1; -1; -3];
  N = length(x); M = length(U);

  StateFF = zeros(L_FF,1); StateFB = zeros(L_FB,1);
  hFF = zeros(L_FF,1);hFF((L_FF-1)/2)=1;
  hFB = zeros(L_FB,1);

  %Pre-allocate vectors
  eqSym = zeros(1,N_train); error = eqSym; yd_train = eqSym;
  coeffsHist = zeros(N_train,length([hFF;hFB]));
  fprintf('Training');

  # Training Sequence Driven EQU-----------------------------------------------
  for n = 1:N_train
    if rem(n,N/50) == 0, fprintf('-');end

    StateFF(2:end) = StateFF(1:end-1);StateFF(1)=x(n);

    eqSym(n) = hFF'*StateFF + hFB'*StateFB; %y

    % Perform detection of a receeived symbol via mapping
    dist = (ones(M,1)*eqSym(n) - U).^2;
    [~,idx] = min(dist); yd_train(n) = U(idx);

    %Find error of equalization
    %error(n) = yd(n) - eqSym(n);
     error(n) = y_train(n) - eqSym(n);

    %Update filter coefficients
    hFF = hFF + stepsize * error(n) * StateFF;
    hFB = hFB + stepsize * error(n) * StateFB;

    %Update the state of FB filter memory
    StateFB(2:end) = StateFB(1:end-1);StateFB(1)=y_train(n);

    %Record coefficients of both filters for output_precision
    coeffsHist(n,:) = [hFF;hFB];

  end

  fprintf('\n');

  # Decisiom Driven EQU---------------------------------------------


  %Pre-allocate vectors
  eqSym = zeros(1,N); error = eqSym; yd = eqSym;
  coeffsHistE = zeros(N_train,length([hFF;hFB]));
fprintf('Equalizing');

for n = 1:N
    if rem(n,N/50) == 0, fprintf('-');end

    StateFF(2:end) = StateFF(1:end-1);StateFF(1)=x(n);

    eqSym(n) = hFF'*StateFF + hFB'*StateFB; %y

    % Perform detection of a receeived symbol via mapping
    dist = (ones(M,1)*eqSym(n) - U).^2;
    [~,idx] = min(dist); yd(n) = U(idx);

    %Find error of equalization
     error(n) = yd(n) - eqSym(n);


    %Update filter coefficients
    hFF = hFF + stepsize * error(n) * StateFF;
    hFB = hFB + stepsize * error(n) * StateFB;

    %Update the state of FB filter memory
    StateFB(2:end) = StateFB(1:end-1);StateFB(1)=yd(n);

    %Record coefficients of both filters for output_precision
    coeffsHistE(n,:) = [hFF;hFB];

  end

  fprintf('\n');
  coeffsHistE = [coeffsHist; coeffsHistE];
