function hEq = MSEEequFun(Rk, Ak, Leq,hChan)
##  %Performs MSE equalization with iterations methods
##  % Rk - channel distorted output samples
##  % Ak - origianl data symbol samples
##  % hChan -
##
##
  LeqH = (Leq-1)/2;
  beta = 0.001;
  L = length(hChan);
  N = length(Rk);
  hEq = zeros(1,Leq);
  for i = LeqH + 1 : length(Rk) - LeqH
    if rem(i,N/50) == 0, fprintf('.');end
     rk = flipud(Rk(i-LeqH:i+LeqH)');
     Ek(i) = Ak(i)-hEq*rk;
     hEq = hEq + beta*Ek(i)*rk';
end
fprintf('\n');
