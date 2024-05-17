function [rcvData, WgtHist, IdxHist] = ViterbiDetectFun(hChan,RxSym,Ndelay)

% FUNCTION: Viterbi MLSE Equalizer detection
% Performs MLSE equalizing by finding the most probable path on Trellis
% diagram via Viterbi algorithm for a specified Impulse Response of channel.
% Implements Viterbi algorithm with Euclidean metric calculation.
%
% Input arguments:
%   hChan - the coefficients vector for Impulse Response of channel
%   RxSym - the sequence of received symbols with Intersymbol Interference (ISI)
%
% Output arguments:
%   rcvData - the received sequence of symbols. The sequence is truncated by
%             L final symbols, where L is the length of channel Impulse Response
%   WgtHist - a matrix of all weights (metrics) history calculated by Viterbi algorithm for each state
%   IdxHist - a matrix of all ORIGIN states to track transitions of Viterbi algorithm for each path


% All values of alphabet: starting from MAX going to MIN
U = [+3 +1 -1 -3]';         % The source symbols alphabet [max ... min]
M = length(U);              % Base of alphabet (number of unique symbols)
N = length(RxSym);          % The length of input data sequence

% Generate state flow diagram for specified channel
L = length(hChan);              % Impulse Response
NumStates = M^(L-1);            % Number of possible states of Trellis Diagram
% Find all possible states and transitions based on specified Depth of IR and alphabet symbols
% LenSym is a number of repeating symbols in a column for each of alphabet symbols (i.e., +1 or -1)
% NumRep is a number of replications for repeating symbols sequence
for k = 1:L
  NumRep = M^(k-1); LenSym = NumStates/NumRep;      % Find numbers of repetitions for a column with number k
  seq = kron(U,ones(LenSym,1));                     % Repeat each Alphabet Symbol LenSym times
  States(:,k) = repmat(seq,NumRep,1);               % Replicate symbols sequence NumRep times and assign a column
end
          
% Calculate according channel output values for state diagram
ChanSym = States * hChan(:);              % Convert channel impulse response vector into column-vector

% Make a matrix corresponding to the Trellis of State Diagram
% Each row of the matrix corresponds to one of all possible states (+1 +1, +1 -1, etc.)
% Each row contains a number of finite values matching the number of alphabet symbols
% Specifically, these values are positioned in column numbers according to the ORIGIN state
WgMtx = Inf*ones(NumStates,NumStates);    % Pre-allocate matrix filled with Inf values
                                          % The state diagram weights will be filled in by each row
for k = 1:NumStates
  idxOutVal = (k-1)*M + (1:M);              % Specify indices for a group of consecutive M values of state channel outputs
                                            % Each k iteration next M values should be read from ChanSym vector
  idxOrigSt = rem(idxOutVal-1,NumStates)+1; % Find indices of ORIGIN states numbers, which periodically repeat
                                            % The periodicy is achieved by using REM function, however, since
                                            % Matlab indexing starts at 1, an extra subtracting and adding 1 is needed
  WgMtx(k,idxOrigSt) = ChanSym(idxOutVal);  % Assign OUTPUT values in Trellis Matrix at positions marked by ORIGIN state numbers
end


% Set initial weight values
W = zeros(NumStates,1);                     % Current symbol's state minimum weights
WgtHist = zeros(NumStates,N);               % A history of weight values (untruncated) for traceback procedure
IdxHist = zeros(NumStates,N);               % A history of ORIGIN states (untruncated) for each transition, for traceback procedure
rcvData = zeros(1,length(RxSym));           % Pre-allocated outuput data array
         
% Update metrics for each received symbol
% Determine number of symbols to be processed 
N = length(RxSym);
for k = 0:N-1
  % Add some activity progress for visualization purpose only
  if rem(k,N/50) == 0, fprintf('.'); end
  
  % Calculate weights (metrics) for all possible transitions on trellis
  % diagram, considering the ORIGIN state weight as well
  Wgt = W' + (RxSym(k+1) - WgMtx).^2;
  
  % For each State, find the minimum weight (metric) and correspanding path, that
  % enters this state.
  [W, idxW] = min(Wgt,[],2);        % W - min weigth, idxW - number of ORIGIN state
                                    % from which the path enters current, k-th state
  % Store updated weight values in history, as well as previous states
  % numbers, from which the path with minimum weight transitions
  WgtHist(:,k+1) = W;               % History of weights (metrics)
  IdxHist(:,k+1) = idxW;            % History of ORIGIN states
end

% TRACEBACK procedure starts here
% So far, all received symbols were procesessed and all paths have been stored
% in a form of history a) weights, b) ORIGIN states

% Calculate the total weight (metric) for each of NumStates possible paths
% This is performed in reverse order, starting with the final symbol's
% metrics

% Initial values of total accumulated weights are taken as the weights of
% final states for each path
WgtSUM = WgtHist(:,N);             % A column vector, size is NumStates (number of candidate paths)

% Iteratively move to the first state in the recorded history, accumulating
% weight for each of the candidate paths
for j = N:-1:2
    % Update path: k state weight + k-1 state weight, corresponding to ORIGIN state
    % Here IdxHist is used as an index to find the ORIGIN state number for each path
    WgtSUM = WgtSUM + WgtHist(IdxHist(:,j),j-1);     
end

% Once all NumStates path accumulated weights (metrics) calculated,
% find the path with MINIMUM accumulated weight (metric)
[~, minIdx] = min(WgtSUM);          % We only want to know the STATE number of final symbol

% Prepare a "map" for detecting symbol values by repeating each symbol of
% alphabet NumStates/M (number of input paths for each state) times
U = kron(U,ones(NumStates/M,1));

% Read from ORIGIN states history the state of symbol, preceding to the
% last one (used as initial value for loop below
Idx = IdxHist(minIdx);

% Iteratively read previous states and traceback to the ORIGIN state of the
% first processed symbol (in reverse order)
for j = N:-1:2                       % Go from final state to the zero state (j=1)
    rcvData(j) = U(Idx);             % MAP symbol value corresponding to ORIGIN state
    Idx = IdxHist(Idx,j-1);          % Read the number of previous symbol's ORIGIN state
end

% Finally, truncate the last L-1 symbols as too unreliable (branching still occurs for final symbols)
rcvData = rcvData(2:end-Ndelay);     % We discard the first symbol too, as we are interested only in
                                     % transitions between states, thus the number of symbols
                                     % is by 1 less than number of ORIGIN states

fprintf('\n');


         
