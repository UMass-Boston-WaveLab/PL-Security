%%UMass Boston Physical Layer Security Channel Model
%Authors: Eric Brown, Clara Gamboa, Dr. K.C. Kirby-Patel
%
%       (SECTION 1) Variables
%       (SECTION 2) Parameters
%       (SECTION 3) Generated Signal
%       (SECTION 4) LPC Estimation 
%       (SECTION 5) ARYule Estimation
%--------------------------------------------

%%  (SECTION 1)
%This is where we define the variables
S = 10;               % # of Scatterers
N = 10;               % # of Eve samples
d = 0.10;             % wavelength
q = 20;               % Number of samples ahead we attempt to predict
Lamb = 1;             %WHAT IS LAMB
t = N+q;              %is the total number of readings
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
%SR is the sensor readings as seen by easdroper

AP = repmat(pSE, [1,t]) + repmat(pAS, [1,t]) + kron(((2*pi)/Lamb)*d*((1:t)-(N+1)/2),cos(tSE));
ASC = exp(1i*AP);
SR = sum(ASC, 1);   
%---------------------------------------------

%% (SECTION 4)
% LPC Estimation for N+q sensors  
% Variables
    x = SR;        %sensor readings
    p = 4;         %p is the order of teh linear preaditionc(FIF filter)
                   %that predicts value of x
                   
%PC is the Predictor coefficients
%est is the estimated signal
%e is the prediction error 
%[arcs, lags] is the autocorrelation sequence of the prediction error
    atest = lpc(x,p);
    est = filter([0 -atest(2:end)],1,x);
    e = x-est;
    [acs,lags] = xcorr(e,'coeff');
    
% Plot of the Original signal vs Estimated signal    
    plot(1:t,x(1:t),1:t,est(1:t),'--'), grid
    title 'Original Signal vs. LPC Estimate'
    xlabel 'Sensors 1 through N+q', ylabel 'Readings'
    legend('Original signal','LPC estimate')
% Plot of the Autocorrelation Prediction error
    plot(lags,acs), grid
    title 'Autocorrelation of the Prediction Error'
    xlabel 'Lags', ylabel 'Normalized value'
%---------------------------------------------

%% (SECTION 5)
% Here we use the ARYule estimation method to predict the N to N+q sensor
% values. 
