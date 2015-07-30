function [PLSIMvsS, PLSIMvsN, PLSIMvsd, PLSIMvsP, PLSIMvsSNR,PLSIMvsf_d,PLSIMvsf_c] = PLSIMmontecarlo(reps)

%PLSIMMONTECARLO. This is built off of the CRLBMONTECARLO built by Prof.
%Kirby-Patel. I have used a fairly good base version of the channel for
%estimation as a default value set. If there is a really long run time then
%I will split each variable into it's own program and make it more
%parallel. -Eric 

%default values
S    = 6;                % # of Scatterers
N    = 200;               % # of sensor array samples
d    = 0.1;              % Spacing between eavesdropper samples in wavelengths
q    = 200;               % Number of samples ahead we attempt to predict
P    = 30;                % Number of complex sinusoids that make up the wireless channel
f_d  = 11000;             % doppler frequency
f_c  = 2400000;           % carrier frequency
SNR = 20;                   % Signal to Noise Ratio in dB.

%test ranges
testS = 1:20;
testN = 50:10:400;
testd = 0.1:0.1:0.5;
testP = 7:30;
testSNR = 10:30; %in dB
%testf_d = 5000:4000:200000;
%testf_c = 1200000:5000:3600000;

%PL_Sim = PL_Security_Sim_pmusic();
q_maxS = zeros(size(testS));
q_maxN = zeros(size(testN));
q_maxd = zeros(size(testd));
% q_maxq = zeros(size(testq));
q_maxP = zeros(size(testP));
q_maxSNR = zeros(size(testSNR));
% q_maxf_c = zeros(size(testf_c));
% q_maxf_d = zeros(size(testf_d));
for kk = 1:reps   
    
    for ii=1:length(testS)
  
        [H, H_hat] = PL_Security_Sim_pmusic(testS(ii), N, d, P,SNR); %tests the scenario
        err = abs((H-H_hat)/rms(H)); % finds error percentage
        q_maxS(ii) = q_maxS(ii)+find(err>0.05, 1)/reps;
    end
    
    for ii=1:length(testN)
        [H, H_hat] = PL_Security_Sim_pmusic(S, testN(ii), d, P,SNR); %tests the scenario
        err = abs((H-H_hat)/rms(H)); % finds error percentage
        q_maxN(ii) = q_maxN(ii) + find(err>0.05, 1)/reps;
    end
    
    for ii=1:length(testd)
        [H, H_hat] = PL_Security_Sim_pmusic(S, N, testd(ii), P,SNR); %tests the scenario
        err = abs((H-H_hat)/rms(H)); % finds error percentage
        q_maxd(ii) = q_maxd(ii) + find(err>0.05, 1)/reps;
    end
    
%    for ii=1:length(testq)
%         [H, H_hat] = PL_Sim(S, N, d, testq(ii), P,f_d,f_c,SNR); %tests the scenario
%         err = abs((H-H_hat)/rms(H)); % finds error percentage
%         q_maxq(ii) = find(err>0.05, 'first');
%     end
    
    for ii=1:length(testP)
        [H, H_hat] = PL_Security_Sim_pmusic(S, N, d, testP(ii),SNR); %tests the scenario
        err = abs((H-H_hat)/rms(H)); % finds error percentage
        q_maxP(ii) = q_maxP(ii) + find(err>0.05, 1)/reps;
    end
    
%     for ii=1:length(testf_d)
%         [H, H_hat] = PL_Sim(testS(ii), N, d, q, P,testf_d(ii),f_c,SNR); %tests the scenario
%         err = abs((H-H_hat)/rms(H)); % finds error percentage
%         q_maxf_d(ii) = find(err>0.05, 'first');
%     end
%     
%     for ii=1:length(testf_c)
%         [H, H_hat] = PL_Sim(testS(ii), N, d, q, P,f_d,testf_c(ii),SNR); %tests the scenario
%         err = abs((H-H_hat)/rms(H)); % finds error percentage
%         q_maxf_c(ii) = find(err>0.05, 'first');
%     end
    for ii=1:length(testSNR)
        [H, H_hat] = PL_Security_Sim_pmusic(S, N, d, P,testSNR(ii)); %tests the scenario
        err = abs((H-H_hat)/rms(H)); % finds error percentage
        q_maxSNR(ii) = q_maxSNR(ii) + find(err>0.05, 1)/reps;
    end
end
figure; 
plot(testS, q_maxS)
xlabel('Number of Scatterers')
ylabel('Average Max Prediction Length')

figure; 
plot(testN, q_maxN)
xlabel('Number of Measurements')
ylabel('Average Max Prediction Length')

figure; 
plot(testd, q_maxd)
xlabel('Space Between Samples (Wavelengths)')
ylabel('Average Max Prediction Length')

figure; 
plot(testP, q_maxP)
xlabel('Number of component Sinusoids in rootMUSIC Estimate')
ylabel('Average Max Prediction Length')

% figure; 
% plot(testf_d, real(PLSIMvsf_d))
% xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
% ylabel('Estimator Minimum Variance')
% 
% figure; 
% plot(testf_c, real(PLSIMvsf_c))
% xlabel(sprintf('Number of Samples Predicted Ahead (sample spacing = %.2f wavelengths)', d))
% ylabel('Estimator Minimum Variance')

figure; 
plot(testSNR, q_maxSNR)
xlabel('Signal to Noise Ratio (dB)')
ylabel('Average Max Prediction Length')
