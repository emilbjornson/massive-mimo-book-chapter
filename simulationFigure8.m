%This Matlab script can be used to generate Figure 8 in the book chapter:
%
%Trinh Van Chien, Emil Björnson, "Massive MIMO Communications," in 5G
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


%Pathloss exponent
kappa = 3.7;

%Define range of SNR values
SNRvalues = [10^(0/10)  10^(20/10)]; %0 dB and 20 dB


%Define the coverage area (as a square with wrap-around). This dimension
%has no impact on  the Monte-Carlo results since we compute metrics that
%are scale invariant.
squareLength = 1000;

%Maximal number of Monte-Carlo setups
monteCarloRealizations = 1000;


%Maximal number of users per cell
K = 100;

%Compute range of number of users
Kvalues = 1:K;

%Select number of BS antennas for SNR value (same here)
Mvalues = [200 200];


%Number of tiers of hexagonals that are simulated, around the desired cell
tiers = 5;

%Coherence interval length
tau_c = 400;

%Percentage of the radius inside the cell where no UEs are allowed
forbiddenRegion = 0.1;

%Define intersite distance in a normalized scale
intersiteDistance = 2;
intersiteDistanceHalf = intersiteDistance/2;

dmax = intersiteDistanceHalf; %Normalized cell radius
dmin = dmax * forbiddenRegion; %Normalized shortest distance from a BS






%%Begin Monte-Carlo simulations

%Placeholders for storing spectral efficiencies from Monte-Carlo
%simulations
averageSEsReuse1_MR = zeros(monteCarloRealizations,length(Kvalues),length(Mvalues));
averageSEsReuse3_MR = zeros(monteCarloRealizations,length(Kvalues),length(Mvalues));
averageSEsReuse4_MR = zeros(monteCarloRealizations,length(Kvalues),length(Mvalues));

averageSEsReuse1_ZF = zeros(monteCarloRealizations,length(Kvalues),length(Mvalues));
averageSEsReuse3_ZF = zeros(monteCarloRealizations,length(Kvalues),length(Mvalues));
averageSEsReuse4_ZF = zeros(monteCarloRealizations,length(Kvalues),length(Mvalues));



%Go through all setups
for n = 1:monteCarloRealizations

    %Display simulation progress
    disp(['Realization ' num2str(n) ' out of ' num2str(monteCarloRealizations)]);
    
    
    %Vector with BS positions
    BSpositions = zeros(4,1);
    
    %Vectors that store which set of pilots that each BS uses (when there
    %is non-universal pilot reuse)
    reusePattern3 = zeros(4,1); %Pilot reuse 3
    reusePattern4 = zeros(4,1); %Pilot reuse 4
    
    
    %BS index of the next BS to be deployed
    itr = 1;
    
    %Deploy hexagonal cells in "tiers" number of tiers
    for alpha1 = 0:tiers
        for alpha2 = 0:tiers-alpha1
            
            if (alpha1 == 0) || (alpha2>=1)
                
                %Compute a BS location according to Eq. (30)
                BSloc = sqrt(3)*alpha2*intersiteDistanceHalf*1i + sqrt(3)*alpha1*intersiteDistanceHalf*exp(1i*pi*(30/180));
                
                
                %Special: The first BS is placed in the origin (this is
                %where the performance is computed)
                if (alpha1 == 0) && (alpha2 == 0)
                    BSpositions(itr) = BSloc;
                    reusePattern3(itr) = 1;
                    reusePattern4(itr) = 1;
                    itr = itr+1;
                    
                else
                    
                    %Compute the current BS location
                    basis = [3*intersiteDistanceHalf/2 0; sqrt(3)*intersiteDistanceHalf/2 sqrt(3)*intersiteDistanceHalf];
                    rotation = [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)];
                    
                    %Deploy BSs in all six directions from the center cell
                    for m = 1:6
                        
                        %Compute the reuse pattern identity for reuse 3
                        if mod(2*alpha1+alpha2,3)==0
                            reusePattern3(itr) = 1;
                        elseif mod(mod(2*alpha1+alpha2,3)+m,2)==1
                            reusePattern3(itr) = 2;
                        elseif mod(mod(2*alpha1+alpha2,3)+m,2)==0
                            reusePattern3(itr) = 3;
                        end
                        
                        %Compute the reuse pattern identity for reuse 4
                        if mod(alpha1,2)==0 && mod(alpha2,2)==0
                            reusePattern4(itr) = 1;
                        elseif m == 1
                            reusePattern4(itr) = 1+mod(2*alpha1+mod(alpha2,2),4);
                        else
                            dims = round(basis\(rotation^(m-1)*basis*[alpha1; alpha2]));
                            reusePattern4(itr) = 1+mod(2*dims(1)+mod(dims(2),2),4);
                        end
                        
                        %Deploy a BS
                        BSpositions(itr) = BSloc;
                        itr = itr+1;
                        
                        %Rotate the BS location to consider the next
                        %direction of the six directions from the center
                        %cell
                        BSloc = BSloc .* exp(1i*2*pi*60/360);
                        
                    end
                end
                
            end
            
        end
    end
    
    %Compute the final number of BSs
    nbrBSs = itr-1;
    
    %Prepare to put out UEs in the cells
    UEpositions = zeros(K,nbrBSs);
    
    %Initiate matrices where first and second order interference are computed
    interference1reuse1 = zeros(K,1);
    interference2reuse1 = zeros(K,1);
    interference1reuse3 = zeros(K,3);
    interference2reuse3 = zeros(K,3);
    interference1reuse4 = zeros(K,4);
    interference2reuse4 = zeros(K,4);
    
    
    %Go through all the cells
    for j = 1:nbrBSs
        
        %Generate UE locations randomly with uniform distribution inside the cells
        nbrToGenerate = K; %Number of UE locations left to generate
        notFinished = true(K,1); %Indices of the UE locations that are left to generate
        
        
        %Iterate the generation of UE locations until all of them are inside a
        %hexagonal cell
        while nbrToGenerate>0
            
            %Generate new UE locations uniformly at random in a circle of radius dmax
            UEpositions(notFinished,j) = sqrt( rand(nbrToGenerate,1)*(dmax^2-dmin^2)+ dmin^2 ) .* exp(1i*2*pi*rand(nbrToGenerate,1));
            
            %Check which UEs that are inside a hexagonal and declare as finished
            finished = checkHexagonal(UEpositions(:,j)',dmax);
            
            %Update which UEs that are left to generate
            notFinished = (finished==false);
            
            %Update how many UEs that are left to generate
            nbrToGenerate = sum(notFinished);
            
        end
        
        %Finalize UE locations by translating them around the serving BS
        UEpositions(:,j) = UEpositions(:,j) + BSpositions(j);
        
        %Compute the distance from the users in cell j to BS j
        distancesSquaredBSj = abs(UEpositions(:,j) - BSpositions(j));
        
        
        %Focus on interference caused to the first BS (in the center)
        l = 1;
        
        %Compute the distance from the users in cell j to BS l
        distancesSquaredBSl = abs(UEpositions(:,j) - BSpositions(l));
        
        
        %Compute inteference terms of the types that show up Eq. (7) in [9]
        interference1reuse1(:,l) = interference1reuse1(:,l) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse1(:,l) = interference2reuse1(:,l) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
        interference1reuse3(:,reusePattern3(j)) = interference1reuse3(:,reusePattern3(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse3(:,reusePattern3(j)) = interference2reuse3(:,reusePattern3(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
        
        interference1reuse4(:,reusePattern4(j)) = interference1reuse4(:,reusePattern4(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(kappa);
        interference2reuse4(:,reusePattern4(j)) = interference2reuse4(:,reusePattern4(j)) + (distancesSquaredBSj./distancesSquaredBSl).^(2*kappa);
        
    end
    
    
    
    %Go through all different number of BS antennas
    for mind = 1:length(Mvalues)
        
        m = Mvalues(mind);
        
        SNR = SNRvalues(mind);
        
        
        %Go through all different number of users
        for ind = 1:length(Kvalues)
            
            k = Kvalues(ind);
            
            
            %Compute different pilot lengths (depending on pilot reuse factor)
            taup1 = k;
            taup3 = 3*k;
            taup4 = 4*k;
            
            %Prepare to store the SINRs
            SINRreuse1_MR = zeros(k,1);
            SINRreuse3_MR = zeros(k,1);
            SINRreuse4_MR = zeros(k,1);

            SINRreuse1_ZF = zeros(k,1);
            SINRreuse3_ZF = zeros(k,1);
            SINRreuse4_ZF = zeros(k,1);
            
            
            %Compute the SINRs in accordance to Lemma 2 in [9] for MR combining
            SINRreuse1_MR(:,1) = m *ones(k,1) ./ ( (sum(interference1reuse1(1:k,1)) + 1/SNR)*(interference1reuse1(1:k,1)+1/(taup1*SNR)) + m*(interference2reuse1(1:k,1)-1));
            SINRreuse3_MR(:,1) = m *ones(k,1) ./ ( (sum(interference1reuse1(1:k,1)) + 1/SNR)*(interference1reuse3(1:k,1)+1/(taup3*SNR)) + m*(interference2reuse3(1:k,1)-1));
            SINRreuse4_MR(:,1) = m *ones(k,1) ./ ( (sum(interference1reuse1(1:k,1)) + 1/SNR)*(interference1reuse4(1:k,1)+1/(taup4*SNR)) + m*(interference2reuse4(1:k,1)-1));
            
            %Compute the SINRs in accordance to Lemma 2 in [9] for ZF combining
            if (m-k)>0 %Check if ZF exists
                SINRreuse1_ZF(:,1) = (m-k) *ones(k,1) ./ ( (sum(interference1reuse1(1:k,1)) - sum(interference2reuse1(1:k,1)./(interference1reuse1(1:k,1)+1/(taup1*SNR))) + 1/SNR)*( interference1reuse1(1:k,1)  +1/(taup1*SNR)) + (m-k)*(interference2reuse1(1:k,1)-1));
                SINRreuse3_ZF(:,1) = (m-k) *ones(k,1) ./ ( (sum(interference1reuse1(1:k,1)) - sum(interference2reuse3(1:k,1)./(interference1reuse3(1:k,1)+1/(taup3*SNR))) + 1/SNR)*( interference1reuse3(1:k,1) +1/(taup3*SNR)) + (m-k)*(interference2reuse3(1:k,1)-1));
                SINRreuse4_ZF(:,1) = (m-k) *ones(k,1) ./ ( (sum(interference1reuse1(1:k,1)) - sum(interference2reuse4(1:k,1)./(interference1reuse4(1:k,1)+1/(taup4*SNR))) + 1/SNR)*( interference1reuse4(1:k,1) +1/(taup4*SNR)) + (m-k)*(interference2reuse4(1:k,1)-1));
            end
            
            
            %Compute the average SEs over the network by Monte-Carlo simulation
            averageSEsReuse1_MR(n,ind,mind) = mean(log2(1+SINRreuse1_MR(:)));
            averageSEsReuse3_MR(n,ind,mind) = mean(log2(1+SINRreuse3_MR(:)));
            averageSEsReuse4_MR(n,ind,mind) = mean(log2(1+SINRreuse4_MR(:)));
            
            averageSEsReuse1_ZF(n,ind,mind) = mean(log2(1+SINRreuse1_ZF(:)));
            averageSEsReuse3_ZF(n,ind,mind) = mean(log2(1+SINRreuse3_ZF(:)));
            averageSEsReuse4_ZF(n,ind,mind) = mean(log2(1+SINRreuse4_ZF(:)));
            
        end
        
    end
    
end




%Plot Figure 8a
figure;
hold on; box on;

for mind = 1:length(Mvalues)
    
    upper1 = Kvalues.*max([(1-Kvalues/tau_c); zeros(size(Kvalues))]).*mean(averageSEsReuse1_MR(:,:,mind),1);
    upper3 = Kvalues.*max([(1-3*Kvalues/tau_c); zeros(size(Kvalues))]).*mean(averageSEsReuse3_MR(:,:,mind),1);
    upper4 = Kvalues.*max([(1-4*Kvalues/tau_c); zeros(size(Kvalues))]).*mean(averageSEsReuse4_MR(:,:,mind),1);
    
    upperMR = max([upper1; upper3; upper4],[],1);
    
    plot(Kvalues,upper1,'r-');
    plot(Kvalues,upper3,'k-.');
    plot(Kvalues,upper4,'b--');
    
end

legend('Pilot reuse: f=1','Pilot reuse: f=3','Pilot reuse: f=4','Location','NorthWest');
xlabel('Number of Users (K)');
ylabel('Spectral Efficiency (bit/s/Hz/cell)');


%Plot Figure 8b
figure;
hold on; box on;

for mind = 1:length(Mvalues)
    
    upper1 = Kvalues.*max([(1-Kvalues/tau_c); zeros(size(Kvalues))]).*mean(averageSEsReuse1_ZF(:,:,mind),1);
    upper3 = Kvalues.*max([(1-3*Kvalues/tau_c); zeros(size(Kvalues))]).*mean(averageSEsReuse3_ZF(:,:,mind),1);
    upper4 = Kvalues.*max([(1-4*Kvalues/tau_c); zeros(size(Kvalues))]).*mean(averageSEsReuse4_ZF(:,:,mind),1);
    
    upperZF = max([upper1; upper3; upper4],[],1);
    
    plot(Kvalues,upper1,'r-');
    plot(Kvalues,upper3,'k-.');
    plot(Kvalues,upper4,'b--');
    
end

legend('Pilot reuse: f=1','Pilot reuse: f=3','Pilot reuse: f=4','Location','NorthWest');
xlabel('Number of Users (K)');
ylabel('Spectral Efficiency (bit/s/Hz/cell)');
