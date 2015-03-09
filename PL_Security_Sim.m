%input # of scatterers 
%input number of samples for E and distance
%for each scatterer
%% This is where we define the variables
S = 10; % # of Scatterers
N = 10; % # of Eve samples
d = .1 % wavelength 
q = 20 % Number of samples ahead we attempt to predict
Lamb = 1 
%%
phiAS = 2*pi*rand(S, 1) %where S is each scatterer (creates a vertical 
                        % array with unique angles for each value in the
                       % array)

phiSE = 2*pi*rand(S, 1)                       


ThetaSE = 2*pi*rand(S, 1); %for each angle required. !!Needs association , S = scat

ThetaAS = 2*pi*rand(S, 1);

%--------------------------
%PSA = exp(j*(phiAS) + (phiSE)) %multiplies phasor values 
  %by exponent and imag 
%PSB = exp(j*(phiBS) + (phiSE))
%--------------------------

% d is specified in terms of lambda, therefore, lambda = 1
% This defines the point on the wave as determined by distance from center
% of E 
%This is the phase term not including any doppler shenanigans
ArrPhase = repmat(phiSE, [1, N+q]) + repmat(phiAS, [1, N+q]) + kron(((2*pi)/Lamb)*d*((1:(N+q))-(N+1)/2),cos(ThetaSE))

ArrScatCont = exp(j*ArrPhase);

%AS = sum(ALL AS) %places phasors into array !!Not sure how this is going to work!!

 
%BS = sum(ALL BS)
