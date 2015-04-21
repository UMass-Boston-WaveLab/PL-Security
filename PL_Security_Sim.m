function [] = PL_Security_Sim(S, N, d, q, p)
%%UMass Boston Physical Layer Security Channel Model
%Authors: Eric Brown, Clara Gamboa, Dr. K.C. Kerby-Patel
%
%       (SECTION 1) Variables
%       (SECTION 2) Parameters
%       (SECTION 3) Generated Signal
%       (SECTION 4) LPC Estimation 
%       (SECTION 5) ARYule Estimation
%--------------------------------------------

%%  (SECTION 1)
%This is where we define the variables
% S = 10;               % # of Scatterers
% N = 100;               % # of Eve samples
% d = 0.10;             % Spacing between eavesdropper samples in wavelengths
% q = 20;               % Number of samples ahead we attempt to predict
Lamb = 1;             % Wavelength = 1 (distances are normalized to the wavelength)
t = N + q;              % is the total number of readings
%---------------------------------------------

%% (SECTION 2)
%In this section we define our wave characteristics. Phase angles and
%thetas

pAS = 2*pi*rand(S, 1);    %where S is each scatterer (creates a vertical
                          %array with unique angles for each value in the
                          %array)
pSE = 2*pi*rand(S, 1);
tSE = 2*pi*rand(S, 1);    %for each angle required.
tAS = 2*pi*rand(S, 1);
%---------------------------------------------

%% (SECTION 3)
% Generated signal for testing, generate what sensors N through N+q should
% be seeing.
% d is specified in terms of lambda, therefore, lambda = 1
% This defines the point on the wave as determined by distance from center
% of E. NOTE: This is the phase term not including any doppler effect

%AP is the Array Phase 
%ASC is the Array of Scatterers Continued???
%H is the full channel as seen by eavesdroper

AP = repmat(pSE, [1,t]) + repmat(pAS, [1,t]) + kron(((2*pi)/Lamb)*d*((1:t)-(N+1)/2),cos(tSE));
ASC = exp(1i*AP);
H = sum(ASC, 1);   
%---------------------------------------------

%% (SECTION 4)
% LPC Estimation for N sensors  
% Variables
    x = H(1:N);        %sensor readings from our listening array
   % p = 10;            %p is the order of the linear preaditions(FIF filter)
                       %that predicts value of x

a = lpc(x, p);
%afilt = filter([0 -a(2:end)],1,x);
afilt = a;
estimates = x;
%a = lpc(x);
%estimates = x;

for Q = 1:q
    temp = predict(afilt, estimates);
    estimates = [estimates temp]; 
end



%PC is the Predictor coefficients
%est is the estimated signal
%e is the prediction error 
%[arcs, lags] is the autocorrelation sequence of the prediction error
%    a = lpc(x,p);
%    est = filter([0 -a(2:end)],1,x);
%    err = x-est;
%    [acs,lags] = xcorr(err,'coeff');
    
R =abs(xcorr(H));
ind = max(find(R>R(length(H))/2)); %max of R always occurs at 0 offset 
clen = ind-length(H); 
figure;
plot((1:length(R))-floor(length(R)/2),R)
title(sprintf('Autocorrelation Function Estimate from Samples; Correlation length = %.0f samples', clen))


% Plot of the Original signal vs Estimated signal  
figure;
    plot(1:t,abs(H(1:t)),1:t,abs(estimates),'--'), grid
    title 'Original Signal vs. LPC Estimate'
    xlabel 'Sensors 1 through N+q', ylabel 'Readings'
    legend('Original signal','LPC estimate')
% Plot of the Autocorrelation Prediction error
%    plot(lags,acs), grid
%    title 'Autocorrelation of the Prediction Error'
%    xlabel 'Lags', ylabel 'Normalized value'
%---------------------------------------------


%% (SECTION 5)
% Here we use the pmusic spectrum analysis method 




end