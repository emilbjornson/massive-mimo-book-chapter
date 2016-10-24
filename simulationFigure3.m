%This Matlab script can be used to generate Figure 3 in the book chapter:
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


%Range of SNR values
SNRdB = 5; %dB scale
SNR = 10.^(SNRdB/10); %Linear scale

%Number of realizations in the Monte Carlo simulations
monteCarloRealizations = 10000;


%Number of users
K = 10;

%Range of number of BS antennas
M = 10:2:100;

%Extract maximum number of antennas
Mmax = max(M);


%%Pre-allocation of matrices for saving simulation results
rateZF_perfect = zeros(length(M),monteCarloRealizations);
rateZF_imperfect_TDD = zeros(length(M),monteCarloRealizations);
rateZF_imperfect_FDD50 = zeros(length(M),monteCarloRealizations);

rateMR_perfect = zeros(length(M),monteCarloRealizations);
rateMR_imperfect_TDD = zeros(length(M),monteCarloRealizations);
rateMR_imperfect_FDD50 = zeros(length(M),monteCarloRealizations);



%Go through all Monte Carlo realizations
for r = 1:monteCarloRealizations
    
    %Generate channel matrix realization
    H = (randn(K,Mmax)+1i*randn(K,Mmax))/sqrt(2);
    
    %Generate noise matrix in channel estimation
    N = (randn(K,Mmax)+1i*randn(K,Mmax))/sqrt(2);
    
    %Generate estimated channel matrix realization 
    %(assuming K-length pilot sequence)
    Hhat = sqrt(K*SNR)/(K*SNR+1)*(sqrt(K*SNR)*H+N);
    
    %Go through all number of BS antennas
    for ind = 1:length(M)
        
        
        %%MR precoding

        %Compute unit-norm MR precoding vectors, with perfect CSI
        MR_perfect = H(:,1:M(ind))'./repmat(sqrt(sum(abs(H(:,1:M(ind))').^2,1)),[M(ind) 1]);
        
        %Compute rate with MR and perfect CSI
        channelGains = abs(H(:,1:M(ind))*MR_perfect).^2;
        rateMR_perfect(ind,r) = sum(log2(1+diag(channelGains)./(sum(channelGains,2)-diag(channelGains)+1/SNR) ));
        
        
        %Compute unit-norm MR precoding vectors, with imperfect CSI
        MR_imperfect = Hhat(:,1:M(ind))'./repmat(sqrt(sum(abs(Hhat(:,1:M(ind))').^2,1)),[M(ind) 1]);
        channelGains = abs(H(:,1:M(ind))*MR_imperfect).^2;
        
        %Compute rate with MR and imperfect CSI
        rateMR_imperfect_TDD(ind,r) = sum(log2(1+diag(channelGains)./(sum(channelGains,2)-diag(channelGains) +1/SNR ) ));
        
        
        
        %%ZF precoding
        
        %Compute unit-norm ZF precoding vectors, with perfect CSI
        ZF_perfect = H(:,1:M(ind))'/(H(:,1:M(ind))*H(:,1:M(ind))');
        ZF_perfect = ZF_perfect./repmat(sqrt(sum(abs(ZF_perfect).^2,1)),[M(ind) 1]);
        
        %Compute rate with ZF and perfect CSI
        channelGains = abs(H(:,1:M(ind))*ZF_perfect).^2;
        rateZF_perfect(ind,r) = sum(log2(1+diag(channelGains)./(sum(channelGains,2)-diag(channelGains)+1/SNR) ));
        

        %Compute unit-norm ZF precoding vectors, with imperfect CSI
        ZF_imperfect = Hhat(:,1:M(ind))'/(Hhat(:,1:M(ind))*Hhat(:,1:M(ind))');
        ZF_imperfect = ZF_imperfect./repmat(sqrt(sum(abs(ZF_imperfect).^2,1)),[M(ind) 1]);

        %Compute rate with ZF and imperfect CSI
        channelGains = abs(H(:,1:M(ind))*ZF_imperfect).^2;
        rateZF_imperfect_TDD(ind,r) = sum(log2(1+diag(channelGains)./(sum(channelGains,2)-diag(channelGains)+1/SNR) ));
        
        
        %Check if the up to 50 FDD downlink pilots are sufficient to
        %estimate all M channel directions
        if M(ind) < 50
            
            %All channel dimensions are estimated
            rateZF_imperfect_FDD50(ind,r) = rateZF_imperfect_TDD(ind,r); 
            rateMR_imperfect_FDD50(ind,r) = rateMR_imperfect_TDD(ind,r); 
           
        elseif M(ind) == 50
            
            %Only the first 50 channel dimensions are estimated
            rateZF_imperfect_FDD50(ind:end,r) = rateZF_imperfect_TDD(ind,r); 
            rateMR_imperfect_FDD50(ind:end,r) = rateMR_imperfect_TDD(ind,r); 
           
        end
        
        
    end
    
end



%Plot Figure 3a
figure; hold on; box on;

plot(M,mean(rateMR_perfect,2),'b','LineWidth',1);
plot(M,mean(rateMR_imperfect_TDD,2),'b--','LineWidth',1);
plot(M,mean(rateMR_imperfect_FDD50,2),'b-.','LineWidth',1);
plot(M,mean(rateMR_imperfect_TDD(1,:),2)*ones(length(M),1),'b:','LineWidth',1);

xlabel('Number of BS Antennas (M)');
ylabel('Spectral Efficiency (bit/s/Hz)');
ylim([0 40]);
legend('Perfect CSI','TDD (\tau_p=10) or FDD (\tau_p=M)','FDD (\tau_p=min(M,50))','FDD (\tau_p=10)','Location','NorthWest');


%Plot Figure 3b
figure; hold on; box on;
plot(M,mean(rateZF_perfect,2),'k','LineWidth',1);
plot(M,mean(rateZF_imperfect_TDD,2),'k--','LineWidth',1);
plot(M,mean(rateZF_imperfect_FDD50,2),'k-.','LineWidth',1);
plot(M,mean(rateZF_imperfect_TDD(1,:),2)*ones(length(M),1),'k:','LineWidth',1);

xlabel('Number of BS Antennas (M)');
ylabel('Spectral Efficiency (bit/s/Hz)');
ylim([0 90]);
legend('Perfect CSI','TDD (\tau_p=10) or FDD (\tau_p=M)','FDD (\tau_p=min(M,50))','FDD (\tau_p=10)','Location','NorthWest');
