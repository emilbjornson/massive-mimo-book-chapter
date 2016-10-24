%This Matlab script can be used to generate Figure 2 in the book chapter:
%
%Trinh Van Chien, Emil Bjornson, ?Massive MIMO Communications,? in 5G
%Mobile Communications, W. Xiang et al. (eds.), pp. 77-116, Springer, 2017.
%
%Download article: -
%
%This is version 1.0 (Last edited: 2016-10-24)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.



%Initialization
close all;
clear;


%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));

%Range of number of BS antennas
Mantennas = 10:10:100;

%Number of users
K = 10;

%Should the optimal beamforming be calculated? (true /false)
%(Note that it is the calculation of optimal beamforming is the main source
%of computational complexity. The running time with optimal beamforming
%takes hours or days, while it only takes a few minutes when only the
%heuristic beamforming schemes are computed. This shows clearly the need
%for simple heuristic beamforming!
computeCapacity = false;

%Number of realizations in the Monte Carlo simulations
monteCarloRealizations = 100;

%Combined channel matrix will be (K x K*M). This matrix gives the
%normalized variance of each channel element
channelVariances = ones(1,K);

%User weights for (unweighted) sum rate computation
weights = ones(K,1);

%Range of SNR values
SNRdB = 5; %dB scale
SNR = K*10.^(SNRdB/10); %Linear scale



%%Pre-allocation of matrices

%Matrices for saving sum rates with different beamforming stategies
sumRateZF = zeros(monteCarloRealizations,length(Mantennas));
sumRateFP = zeros(monteCarloRealizations,length(Mantennas));
sumrateCAPACITY = zeros(monteCarloRealizations,length(Mantennas));



%Go through the different number of BS antennas
for n = 1:length(Mantennas)
    
    %Extract the current number of antennas
    M = Mantennas(n);
    
    %Pre-generation of Rayleigh fading channel realizations (unit variance)
    Hall = (randn(K,M,monteCarloRealizations)+1i*randn(K,M,monteCarloRealizations))/sqrt(2);
    
    %Output the progress
    disp(['Progress: M = ' num2str(M)]);
    
    %Go through all channel realizations
    for m = 1:monteCarloRealizations
        
        
        
        %Generate channel matrix for m:th realization
        H = repmat(sqrt(channelVariances)',[1 M]) .* Hall(:,:,m);
        
        
        %Compute normalized beamforming vectors for MR
        wMRT = functionMRT(H);
        
        %Compute normalized beamforming vectors for ZF
        wZF = functionZFBF(H);
        
        

        %Calculate power allocation with MR (using Theorem 3.5 in [3])
        rhos = diag(abs(H*wMRT).^2)';
        powerAllocationFP = functionHeuristicPowerAllocation(rhos,SNR,weights);
        
        %Calculate sum rate without interference (by removing interference
        %from MR)
        W = kron(sqrt(powerAllocationFP),ones(M,1)).*wMRT;
        channelGains = abs(H*W).^2;
        signalGains = diag(channelGains);
        
        rates = log2(1+signalGains);
        sumRateFP(m,n) = weights'*rates;
        
        
        
        %Calculate power allocation with ZF (using Theorem 3.5 in [3])
        rhos = diag(abs(H*wZF).^2)';
        powerAllocationwZFBF = functionHeuristicPowerAllocation(rhos,SNR,weights);
        
        %Calculate sum rate with ZFBF
        W = kron(sqrt(powerAllocationwZFBF),ones(M,1)).*wZF;
        channelGains = abs(H*W).^2;
        signalGains = diag(channelGains);
        interferenceGains = sum(channelGains,2)-signalGains;
        rates = log2(1+signalGains./(interferenceGains+1));
        sumRateZF(m,n) = weights'*rates;
        

        %Compute sum rate using the capacity
        if computeCapacity == true
            
            sumrateCAPACITY(m,n) = real(function_capacity_broadcast(H,K,SNR));
            
        end
        
        
    end
    
    
    
end


%Plot simulation results
figure; hold on; box on;

plot(Mantennas,mean(sumRateFP,1),'k--','LineWidth',1);
plot(Mantennas,mean(sumrateCAPACITY,1),'ro-','LineWidth',1);
plot(Mantennas,mean(sumRateZF,1),'b*--','LineWidth',1);

legend('No interference','Non-linear: Sum capacity','Linear: ZF','Location','SouthEast');

xlabel('Number of BS Antennas (M) ');
ylabel('Spectral Efficiency (bit/s/Hz/cell)');
ylim([0 90]);
