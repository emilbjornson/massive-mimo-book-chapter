function [muValues1Mean,muValues2Mean,reuseMu1Mean,reuseMu1Mean2,reuseMu1MeanNext,reuseMu1Mean2Next,reuseMu2Mean,reuseMuMeanVariance,muValues1Worst,muValues2Worst,reuseMu1Worst,reuseMu1Worst2,reuseMu1WorstNext,reuseMu1Worst2Next,reuseMu2Worst,reuseMuWorstVariance,muValues1Best,muValues2Best,reuseMu1Best,reuseMu1Best2,reuseMu1BestNext,reuseMu1Best2Next,reuseMu2Best,reuseMuBestVariance,reuseFactor] = computeEnvironment(kappa,forbiddenRegion,monteCarloUEs)
%This function performs Monte Carlo simulations to compute various sums of
%the mu-parameters from Eq. (18) in the article:
% 
%Emil Bjornson, Erik G. Larsson, Merouane Debbah, "Massive MIMO for Maximal
%Spectral Efficiency: How Many Users and Pilots Should Be Allocated?,"
%vol. 15, no. 2, pp. 1293-1308, February 2016.
%
%Download article: http://arxiv.org/pdf/1412.7102
%
%This is version 1.1 (Last edited: 2016-08-22)
% 
%License: This code is licensed under the GPLv2 license. If you in any way 
%use this code for research that results in publications, please cite our 
%original article listed above. 
%
%INPUT
%kappa           = Pathloss exponent (typically in the range 2 - 6).
%forbiddenRegion = Percentage of the radius inside the cell where no UEs are allowed
%monteCarloUEs   = Number of random UE locations per cell (Optional)
%
%OUTPUT
%muValues1Mean        = Sum of mu^(1) over the whole network in mean case
%muValues2Mean        = Sum of mu^(2) over the whole network in mean case
%reuseMu1Mean         = Sum of mu^(1) for cells with same pilots as in origin in mean case (for different reuse factors)
%reuseMu1Mean2        = Sum of (mu^(1))^2 for cells with same pilots as in origin in mean case (for different reuse factors)
%reuseMu1MeanNext     = Sum of mu^(1) for cells with same pilots as in neighboring cell to origin in mean case (for different reuse factors)
%reuseMu1Mean2Next    = Sum of (mu^(1))^2 for cells with same pilots as in neighboring cell to origin in mean case (for different reuse factors)
%reuseMu2Mean         = Sum of mu^(2) for cells with same pilots as in origin in mean case (for different reuse factors)
%reuseMuMeanVariance  = Sum of mu^(2)-(mu^(1))^2 for cells with same pilots as in origin in mean case (for different reuse factors)
%muValues1Worst       = Sum of mu^(1) over the whole network in worst case
%muValues2Worst       = Sum of mu^(2) over the whole network in worst case
%reuseMu1Worst        = Sum of mu^(1) for cells with same pilots in worst case (for different reuse factors)
%reuseMu1Worst2       = Sum of (mu^(1))^2 for cells with same pilots in worst case (for different reuse factors)
%reuseMu1WorstNext    = Sum of mu^(1) for cells with same pilots as in neighboring cell to origin in worst case (for different reuse factors)
%reuseMu1Worst2Next   = Sum of (mu^(1))^2 for cells with same pilots as in neighboring cell to origin in worst case (for different reuse factors)
%reuseMu2Worst        = Sum of mu^(2) for cells with same pilots as in origin in worst case (for different reuse factors)
%reuseMuWorstVariance = Sum of mu^(2)-(mu^(1))^2 for cells with same pilots as in origin in worst case (for different reuse factors)
%muValues1Best        = Sum of mu^(1) over the whole network in best case
%muValues2Best        = Sum of mu^(2) over the whole network in best case
%reuseMu1Best         = Sum of mu^(1) for cells with same pilots as in origin in best case (for different reuse factors)
%reuseMu1Best2        = Sum of (mu^(1))^2 for cells with same pilots as in origin in best case (for different reuse factors)
%reuseMu1BestNext     = Sum of mu^(1) for cells with same pilots as in neighboring cell to origin in best case (for different reuse factors)
%reuseMu1Best2Next    = Sum of (mu^(1))^2 for cells with same pilots as in neighboring cell to origin in best case (for different reuse factors)
%reuseMu2Best         = Sum of mu^(2) for cells with same pilots as in origin in best case (for different reuse factors)
%reuseMuBestVariance  = Sum of mu^(2)-(mu^(1))^2 for cells with same pilots as in origin in best case (for different reuse factors)
%reuseFactor          = Matrix with the considered reuse factors


%Set number of UE locations in the Monte Carlo simulations, if this number
%is not set as an input parameter
if nargin<3
    monteCarloUEs = 1000000;
end


%Define matrix for storing UE locations
UElocations = zeros(1,monteCarloUEs);


%Define cell dimensions (the unit or exact size doesn't matter since
%everything is the mu-parameters are scale invariant)
intersiteDistance = 0.5; %Distance between neighboring BSs

dmax = intersiteDistance/2; %Cell radius
dmin = dmax * forbiddenRegion; %Shortest distance from a BS


%Generate UE locations randomly with uniform distribution inside the cells
nbrToGenerate = monteCarloUEs; %Number of UE locations left to generate
notFinished = true(monteCarloUEs,1); %Indices of the UE locations that are left to generate


%Iterate the generation of UE locations until all of them are inside a
%hexagonal cell
while nbrToGenerate>0
    
    %Generate new UE locations uniformly at random in a circle of radius dmax
    UElocations(1,notFinished) = sqrt( rand(1,nbrToGenerate)*(dmax^2-dmin^2)+ dmin^2 ) .* exp(1i*2*pi*rand(1,nbrToGenerate));
    
    %Check which UEs that are inside a hexagonal and declare as finished
    finished = checkHexagonal(UElocations(1,:)',dmax);
    
    %Update which UEs that are left to generate
    notFinished = (finished==false);
    
    %Update how many UEs that are left to generate
    nbrToGenerate = sum(notFinished);
    
end


%Angle between each edge point (360/6 = 60)
baseAngle = 60; 

%Select how many tiers of BSs should be considered around the cell of
%interest
howFar = 5;


%Placeholders for storing results for the mean interference case
muValues1Mean = zeros(howFar+1,howFar+1);
muValues2Mean = zeros(howFar+1,howFar+1);
muValues1Mean(1,1) = 1; %The mu-values in Eq. (18) are one within a cell
muValues2Mean(1,1) = 1; %The mu-values in Eq. (18) are one within a cell

reuseMu1Mean = zeros(howFar+1,howFar+1);
reuseMu1Mean2 = zeros(howFar+1,howFar+1);
reuseMu1MeanNext = zeros(howFar+1,howFar+1);
reuseMu1Mean2Next = zeros(howFar+1,howFar+1);

reuseMu2Mean = zeros(howFar+1,howFar+1);
reuseMuMeanVariance = zeros(howFar+1,howFar+1);


%Placeholders for storing results for the worst interference case
muValues1Worst = zeros(howFar+1,howFar+1);
muValues2Worst = zeros(howFar+1,howFar+1);
muValues1Worst(1,1) = 1; %The mu-values in Eq. (18) are one within a cell
muValues2Worst(1,1) = 1; %The mu-values in Eq. (18) are one within a cell

reuseMu1Worst = zeros(howFar+1,howFar+1);
reuseMu1Worst2 = zeros(howFar+1,howFar+1);
reuseMu1WorstNext = zeros(howFar+1,howFar+1);
reuseMu1Worst2Next = zeros(howFar+1,howFar+1);

reuseMu2Worst = zeros(howFar+1,howFar+1);
reuseMuWorstVariance = zeros(howFar+1,howFar+1);


%Placeholders for storing results for the best interference case
muValues1Best = zeros(howFar+1,howFar+1);
muValues2Best = zeros(howFar+1,howFar+1);
muValues1Best(1,1) = 1; %The mu-values in Eq. (18) are one within a cell
muValues2Best(1,1) = 1; %The mu-values in Eq. (18) are one within a cell

reuseMu1Best = zeros(howFar+1,howFar+1);
reuseMu1Best2 = zeros(howFar+1,howFar+1);
reuseMu1BestNext = zeros(howFar+1,howFar+1);
reuseMu1Best2Next = zeros(howFar+1,howFar+1);

reuseMu2Best = zeros(howFar+1,howFar+1);
reuseMuBestVariance = zeros(howFar+1,howFar+1);


%Placeholder for storing reuse factors
reuseFactor = zeros(howFar+1,howFar+1);

%Define the position of one of the neighboring cells, as seen from the
%origin.
nextNeighbor = sqrt(3)*dmax*exp(1i*pi*(30/180));


%Go through neighboring cells at different distances using the
%parameterization in Eq. (31). Only one search direction is considered, but
%there are six neighboring cells at the same distance.
for alpha1 = 1:1:howFar
    for alpha2 = 0:1:howFar
        
        %Put out another BS using coordinates u and v
        BSlocations = sqrt(3)*alpha1*dmax*exp(1i*pi*(30/180)) + sqrt(3)*alpha2*dmax*1i; %Coordinates based on Eq. (31)
        

        %Compute the reuse factor for the hexagonal topology when the
        %current neighboring cell is the first one that reuses the same
        %pilot sequences
        reuseFactor(alpha1+1,alpha2+1) = alpha1^2+alpha2^2+alpha1*alpha2;
        
        
        %Mean interference
        
        %Compute the mu^(1) and mu^(2) values for mean interference
        %seen from the base station in the origin
        muValues1Mean(alpha1+1,alpha2+1) = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocations)).^kappa);
        muValues2Mean(alpha1+1,alpha2+1) = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocations)).^(2*kappa));
        
        %Store the sum of interference for the cells that have the same
        %pilot sequences as the cell in the origin, when the pilot reuse
        %factor is alpha1^2+alpha2^2+alpha1*alpha2
        reuseMu1Mean(alpha1+1,alpha2+1) = reuseMu1Mean(alpha1+1,alpha2+1) + muValues1Mean(alpha1+1,alpha2+1); %Sum of mu^(1)
        reuseMu1Mean2(alpha1+1,alpha2+1) = reuseMu1Mean2(alpha1+1,alpha2+1) + muValues1Mean(alpha1+1,alpha2+1)^2; %Sum of (mu^(1))^2
        reuseMu2Mean(alpha1+1,alpha2+1) = reuseMu2Mean(alpha1+1,alpha2+1) + muValues2Mean(alpha1+1,alpha2+1); %Sum of mu^(2) 
        reuseMuMeanVariance(alpha1+1,alpha2+1) = reuseMuMeanVariance(alpha1+1,alpha2+1) + muValues2Mean(alpha1+1,alpha2+1) - muValues1Mean(alpha1+1,alpha2+1)^2; %Sum of mu^(2)-(mu^(1))^2
        
        %Store the sum of interference for the cells that have the same
        %pilot sequences as one of the neighbors of the cell in the origin, 
        %when the pilot reuse factor is alpha1^2+alpha2^2+alpha1*alpha2
        newMu1ReuseNext = mean((abs(UElocations(:))./abs(UElocations(:)+nextNeighbor)).^kappa);
        newMu1ReuseNextOneStep = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocations+nextNeighbor)).^kappa);
        reuseMu1MeanNext(alpha1+1,alpha2+1) = reuseMu1MeanNext(alpha1+1,alpha2+1) + newMu1ReuseNext + newMu1ReuseNextOneStep;
        reuseMu1Mean2Next(alpha1+1,alpha2+1) = reuseMu1Mean2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2 + newMu1ReuseNextOneStep.^2;
        
        
        %Worst-case interference
        
        %Compute the mu^(1) and mu^(2) values for worst-case interference,
        %seen from the base station in the origin
        muValues1Worst(alpha1+1,alpha2+1) = max((abs(UElocations(:))./abs(UElocations(:)+BSlocations)).^kappa);
        muValues2Worst(alpha1+1,alpha2+1) = max((abs(UElocations(:))./abs(UElocations(:)+BSlocations)).^(2*kappa));
        
        %Store the sum of interference for the cells that have the same
        %pilot sequences as the cell in the origin, when the pilot reuse
        %factor is alpha1^2+alpha2^2+alpha1*alpha2
        reuseMu1Worst(alpha1+1,alpha2+1) = reuseMu1Worst(alpha1+1,alpha2+1) + muValues1Worst(alpha1+1,alpha2+1); %Sum of mu^(1)
        reuseMu1Worst2(alpha1+1,alpha2+1) = reuseMu1Worst2(alpha1+1,alpha2+1) + muValues1Worst(alpha1+1,alpha2+1).^2;  %Sum of (mu^(1))^2
        reuseMu2Worst(alpha1+1,alpha2+1) = reuseMu2Worst(alpha1+1,alpha2+1) + muValues2Worst(alpha1+1,alpha2+1); %Sum of mu^(2)
        reuseMuWorstVariance(alpha1+1,alpha2+1) = reuseMuWorstVariance(alpha1+1,alpha2+1) + muValues2Worst(alpha1+1,alpha2+1) - muValues1Worst(alpha1+1,alpha2+1)^2; %Sum of mu^(2)-(mu^(1))^2
        
        %Store the sum of interference for the cells that have the same
        %pilot sequences as one of the neighbors of the cell in the origin, 
        %when the pilot reuse factor is alpha1^2+alpha2^2+alpha1*alpha2
        newMu1ReuseNext = max((abs(UElocations(:))./abs(UElocations(:)+nextNeighbor)).^kappa);
        newMu1ReuseNextOneStep = max((abs(UElocations(:))./abs(UElocations(:)+BSlocations+nextNeighbor)).^kappa);
        reuseMu1WorstNext(alpha1+1,alpha2+1) = reuseMu1WorstNext(alpha1+1,alpha2+1) + newMu1ReuseNext + newMu1ReuseNextOneStep;
        reuseMu1Worst2Next(alpha1+1,alpha2+1) = reuseMu1Worst2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2 + newMu1ReuseNextOneStep.^2;
        
        
        
        %Best-case interference
        
        %Compute the mu^(1) and mu^(2) values for best-case interference
        %seen from the base station in the origin
        muValues1Best(alpha1+1,alpha2+1) = min((abs(UElocations(:))./abs(UElocations(:)+BSlocations)).^kappa);
        muValues2Best(alpha1+1,alpha2+1) = min((abs(UElocations(:))./abs(UElocations(:)+BSlocations)).^(2*kappa));
        
        %Store the sum of interference for the cells that have the same
        %pilot sequences as the cell in the origin, when the pilot reuse
        %factor is alpha1^2+alpha2^2+alpha1*alpha2
        reuseMu1Best(alpha1+1,alpha2+1) = reuseMu1Best(alpha1+1,alpha2+1) + muValues1Best(alpha1+1,alpha2+1); %Sum of mu^(1)
        reuseMu1Best2(alpha1+1,alpha2+1) = reuseMu1Best2(alpha1+1,alpha2+1) + muValues1Best(alpha1+1,alpha2+1).^2; %Sum of (mu^(1))^2
        reuseMu2Best(alpha1+1,alpha2+1) = reuseMu2Best(alpha1+1,alpha2+1) + muValues2Best(alpha1+1,alpha2+1); %Sum of mu^(2)
        reuseMuBestVariance(alpha1+1,alpha2+1) = reuseMuBestVariance(alpha1+1,alpha2+1) + muValues2Best(alpha1+1,alpha2+1) - muValues1Best(alpha1+1,alpha2+1)^2; %Sum of mu^(2)-(mu^(1))^2
        
        %Store the sum of interference for the cells that have the same
        %pilot sequences as one of the neighbors of the cell in the origin, 
        %when the pilot reuse factor is alpha1^2+alpha2^2+alpha1*alpha2
        newMu1ReuseNext = min((abs(UElocations(:))./abs(UElocations(:)+nextNeighbor)).^kappa);
        newMu1ReuseNextOneStep = min((abs(UElocations(:))./abs(UElocations(:)+BSlocations+nextNeighbor)).^kappa);
        reuseMu1BestNext(alpha1+1,alpha2+1) = reuseMu1BestNext(alpha1+1,alpha2+1) + newMu1ReuseNext + newMu1ReuseNextOneStep;
        reuseMu1Best2Next(alpha1+1,alpha2+1) = reuseMu1Best2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2 + newMu1ReuseNextOneStep.^2;
        
        
        %Consider the next two cells with the same reuse factor (there are
        %two neigbors instead of one in the second interfering tier)
        for index = 0:1
            
            %Compute location of the next BS that use the same reuse factor
            BSlocation2 = BSlocations + BSlocations*exp(1i*pi*((index*baseAngle)/180));
            
            
            %Mean interference
            
            %Compute the mu^(1) and mu^(2) values for mean interference
            %seen from the base station in the origin
            newMu1Reuse = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocation2)).^kappa);
            newMu2Reuse = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocation2)).^(2*kappa));
            newMu1ReuseNext = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocation2+nextNeighbor)).^kappa); 
            
            reuseMu1Mean(alpha1+1,alpha2+1) = reuseMu1Mean(alpha1+1,alpha2+1) + newMu1Reuse; %Add to sum of mu^(1)
            reuseMu1Mean2(alpha1+1,alpha2+1) = reuseMu1Mean2(alpha1+1,alpha2+1) + newMu1Reuse.^2; %Add to sum of (mu^(1))^2
            reuseMu2Mean(alpha1+1,alpha2+1) = reuseMu2Mean(alpha1+1,alpha2+1) + newMu2Reuse; %Add to sum of mu^(2)
            reuseMuMeanVariance(alpha1+1,alpha2+1) = reuseMuMeanVariance(alpha1+1,alpha2+1) + newMu2Reuse - newMu1Reuse^2; %Add to sum of mu^(2)-(mu^(1))^2
            reuseMu1MeanNext(alpha1+1,alpha2+1) = reuseMu1MeanNext(alpha1+1,alpha2+1) + newMu1ReuseNext; %Add to sum of mu^(1) using pilot of neighboring cell
            reuseMu1Mean2Next(alpha1+1,alpha2+1) = reuseMu1Mean2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2; %Add to sum of (mu^(1))^2 using pilot of neighboring cell
            
            
            %Worst-case interference
            
            %Compute the mu^(1) and mu^(2) for the next BS, for worst-case
            %interference seen from the base station in the origin
            newMu1Reuse = max((abs(UElocations(:))./abs(UElocations(:)+BSlocation2)).^kappa);
            newMu2Reuse = max((abs(UElocations(:))./abs(UElocations(:)+BSlocation2)).^(2*kappa));
            newMu1ReuseNext = max((abs(UElocations(:))./abs(UElocations(:)+BSlocation2+nextNeighbor)).^kappa);
            
            %Store the results
            reuseMu1Worst(alpha1+1,alpha2+1) = reuseMu1Worst(alpha1+1,alpha2+1) + newMu1Reuse; %Add to sum of mu^(1)
            reuseMu1Worst2(alpha1+1,alpha2+1) = reuseMu1Worst2(alpha1+1,alpha2+1) + newMu1Reuse.^2; %Add to sum of (mu^(1))^2
            reuseMu2Worst(alpha1+1,alpha2+1) = reuseMu2Worst(alpha1+1,alpha2+1) + newMu2Reuse; %Add to sum of mu^(2)
            reuseMuWorstVariance(alpha1+1,alpha2+1) = reuseMuWorstVariance(alpha1+1,alpha2+1) + newMu2Reuse - newMu1Reuse^2; %Add to sum of mu^(2)-(mu^(1))^2
            reuseMu1WorstNext(alpha1+1,alpha2+1) = reuseMu1WorstNext(alpha1+1,alpha2+1) + newMu1ReuseNext; %Add to sum of mu^(1) using pilot of neighboring cell
            reuseMu1Worst2Next(alpha1+1,alpha2+1) = reuseMu1Worst2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2; %Add to sum of (mu^(1))^2 using pilot of neighboring cell
            
            
            %Best-case interference
            
            %Compute the mu^(1) and mu^(2) for the next BS, for best-case
            %interference seen from the base station in the origin
            newMu1Reuse = min((abs(UElocations(:))./abs(UElocations(:)+BSlocation2)).^kappa);
            newMu2Reuse = min((abs(UElocations(:))./abs(UElocations(:)+BSlocation2)).^(2*kappa));
            newMu1ReuseNext = min((abs(UElocations(:))./abs(UElocations(:)+BSlocation2+nextNeighbor)).^kappa);
            
            %Store the results
            reuseMu1Best(alpha1+1,alpha2+1) = reuseMu1Best(alpha1+1,alpha2+1) + newMu1Reuse; %Add to sum of mu^(1)
            reuseMu1Best2(alpha1+1,alpha2+1) = reuseMu1Best2(alpha1+1,alpha2+1) + newMu1Reuse.^2; %Add to sum of (mu^(1))^2
            reuseMu2Best(alpha1+1,alpha2+1) = reuseMu2Best(alpha1+1,alpha2+1) + newMu2Reuse; %Add to sum of mu^(2)
            reuseMuBestVariance(alpha1+1,alpha2+1) = reuseMuBestVariance(alpha1+1,alpha2+1) + newMu2Reuse - newMu1Reuse^2; %Add to sum of mu^(2)-(mu^(1))^2
            reuseMu1BestNext(alpha1+1,alpha2+1) = reuseMu1BestNext(alpha1+1,alpha2+1) + newMu1ReuseNext; %Add to sum of mu^(1) using pilot of neighboring cell
            reuseMu1Best2Next(alpha1+1,alpha2+1) = reuseMu1Best2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2; %Add to sum of (mu^(1))^2 using pilot of neighboring cell
            
            
            
            %Consider the next three cells with the same reuse factor
            %(there are three neigbors instead of two in the third interfering tier)
            for index2 = index:1
                
                %Compute location of the next BS that use the same reuse factor
                BSlocation3 = BSlocation2 + BSlocations*exp(1i*pi*((index2*baseAngle)/180));
                
                
                %Mean interference
                
                %Compute the mu^(1) and mu^(2) for the next BS, for mean
                %interference seen from the base station in the origin
                newMu1Reuse = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocation3)).^kappa);
                newMu2Reuse = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocation3)).^(2*kappa));
                newMu1ReuseNext = mean((abs(UElocations(:))./abs(UElocations(:)+BSlocation3+nextNeighbor)).^kappa);
                
                reuseMu1Mean(alpha1+1,alpha2+1) = reuseMu1Mean(alpha1+1,alpha2+1) + newMu1Reuse; %Add to sum of mu^(1)
                reuseMu1Mean2(alpha1+1,alpha2+1) = reuseMu1Mean2(alpha1+1,alpha2+1) + newMu1Reuse.^2; %Add to sum of (mu^(1))^2
                reuseMu2Mean(alpha1+1,alpha2+1) = reuseMu2Mean(alpha1+1,alpha2+1) + newMu2Reuse; %Add to sum of mu^(2)
                reuseMuMeanVariance(alpha1+1,alpha2+1) = reuseMuMeanVariance(alpha1+1,alpha2+1) + newMu2Reuse - newMu1Reuse^2; %Add to sum of mu^(2)-(mu^(1))^2
                reuseMu1MeanNext(alpha1+1,alpha2+1) = reuseMu1MeanNext(alpha1+1,alpha2+1) + newMu1ReuseNext; %Add to sum of mu^(1) using pilot of neighboring cell
                reuseMu1Mean2Next(alpha1+1,alpha2+1) = reuseMu1Mean2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2; %Add to sum of (mu^(1))^2 using pilot of neighboring cell
                
                
                %Worst-case interference
                
                %Compute the mu^(1) and mu^(2) for the next BS, for
                %worst-case interference seen from the base station in the origin
                newMu1Reuse = max((abs(UElocations(:))./abs(UElocations(:)+BSlocation3)).^kappa);
                newMu2Reuse = max((abs(UElocations(:))./abs(UElocations(:)+BSlocation3)).^(2*kappa));
                newMu1ReuseNext = max((abs(UElocations(:))./abs(UElocations(:)+BSlocation3+nextNeighbor)).^kappa);
                
                %Store the results
                reuseMu1Worst(alpha1+1,alpha2+1) = reuseMu1Worst(alpha1+1,alpha2+1) + newMu1Reuse; %Add to sum of mu^(1)
                reuseMu1Worst2(alpha1+1,alpha2+1) = reuseMu1Worst2(alpha1+1,alpha2+1) + newMu1Reuse.^2; %Add to sum of (mu^(1))^2
                reuseMu2Worst(alpha1+1,alpha2+1) = reuseMu2Worst(alpha1+1,alpha2+1) + newMu2Reuse; %Add to sum of mu^(2)
                reuseMuWorstVariance(alpha1+1,alpha2+1) = reuseMuWorstVariance(alpha1+1,alpha2+1) + newMu2Reuse - newMu1Reuse^2; %Add to sum of mu^(2)-(mu^(1))^2
                reuseMu1WorstNext(alpha1+1,alpha2+1) = reuseMu1WorstNext(alpha1+1,alpha2+1) + newMu1ReuseNext; %Add to sum of mu^(1) using pilot of neighboring cell
                reuseMu1Worst2Next(alpha1+1,alpha2+1) = reuseMu1Worst2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2; %Add to sum of (mu^(1))^2 using pilot of neighboring cell
                
                
                %Best-case interference
                
                %Compute the mu^(1) and mu^(2) for the next BS, for
                %best-case interference seen from the base station in the origin
                newMu1Reuse = min((abs(UElocations(:))./abs(UElocations(:)+BSlocation3)).^kappa);
                newMu2Reuse = min((abs(UElocations(:))./abs(UElocations(:)+BSlocation3)).^(2*kappa));
                newMu1ReuseNext = min((abs(UElocations(:))./abs(UElocations(:)+BSlocation3+nextNeighbor)).^kappa);
                
                %Store the results
                reuseMu1Best(alpha1+1,alpha2+1) = reuseMu1Best(alpha1+1,alpha2+1) + newMu1Reuse; %Add to sum of mu^(1)
                reuseMu1Best2(alpha1+1,alpha2+1) = reuseMu1Best2(alpha1+1,alpha2+1) + newMu1Reuse.^2; %Add to sum of (mu^(1))^2
                reuseMu2Best(alpha1+1,alpha2+1) = reuseMu2Best(alpha1+1,alpha2+1) + newMu2Reuse; %Add to sum of mu^(2)
                reuseMuBestVariance(alpha1+1,alpha2+1) = reuseMuBestVariance(alpha1+1,alpha2+1) + newMu2Reuse - newMu1Reuse^2; %Add to sum of mu^(2)-(mu^(1))^2
                reuseMu1BestNext(alpha1+1,alpha2+1) = reuseMu1BestNext(alpha1+1,alpha2+1) + newMu1ReuseNext; %Add to sum of mu^(1) using pilot of neighboring cell
                reuseMu1Best2Next(alpha1+1,alpha2+1) = reuseMu1Best2Next(alpha1+1,alpha2+1) + newMu1ReuseNext.^2; %Add to sum of (mu^(1))^2 using pilot of neighboring cell
                
            end
            
        end
        
    end
    
end
