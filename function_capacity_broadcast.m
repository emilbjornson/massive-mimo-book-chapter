function [sumCapacity,p] = function_capacity_broadcast(H,Pmax)
%This is function computes the downlink sum capacity and is used in the
%article: 
%
%Emil Bj?rnson, Erik G. Larsson, Thomas L. Marzetta, "Massive MIMO: Ten
%Myths and One Critical Question," IEEE Communications Magazine, vol. 54, 
%no. 2, pp. 114-123, February 2016. 
%
%Download article: http://arxiv.org/pdf/1503.06854
%
%This is version 1.01 (Last edited: 2016-04-27)
%
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%INPUT:
%H           = M x K channel matrix
%Pmax        = Maximal transmit power
%
%OUTPUT:
%sumCapacity = Sum capacity
%p           = K-vector with optimized power allocation over the terminals


%Extract system dimensions
K = size(H,1); %Number of terminals
M = size(H,2); %Number of service antennas

%Solve the capacity-achieving power optimization problem using CVX, by
%solving the convex problem stated in Theorem in "Sum Capacity of the
%Vector Gaussian Broadcast Channel and Uplink?Downlink Duality" by Pramod
%Viswanath and David Tse 
cvx_begin
cvx_quiet(true); % this suppresses screen output from the solver
variable p(K);
minimize det_inv(eye(M)+H'*diag(p)*H);
subject to
    sum(p)<=Pmax
    min(p)>=0
cvx_end

%Compute the sum capacity using the optimized power allocation
sumCapacity = real(log2(det(eye(M)+H'*diag(p)*H)));
