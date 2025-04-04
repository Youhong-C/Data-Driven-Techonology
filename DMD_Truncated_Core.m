function [A_tilde, dmd_modes, dmd_amplitudes, dmd_evals, dmd_evecs, U] = DMD_Truncated_Core(X,Y)

%-----------------------------------------------
% Name of file : DMD_Truncated_Core.m
% 
% Created   : 19/02/2025
%
% Purpose   : Core Implementation of the DMD algorithm
%             with Truncated Singular Values
%           
% Author    : Debraj Bhattacharjee
%
% Copyright : Debraj Bhattacharjee, 2025
%------------------------------------------------   

    % Compute economy SVD of X
    [U,S,V]=svd(X,'econ');

    % FInd length to truncate
    h = figure;
    bar(diag(S));

    r = input("Enter truncation rank: ");
    close(h);

    % Truncate the matrices obtained from SVD
    Ur=U(:,1:r);
    Sr=S(1:r,1:r);
    Vr=V(:,1:r);

    % Construct low dimensional DMD matrix A_tilde
    A_tilde=Ur'*Y*Vr/Sr;
    
    % Compute DMD eigenvalues and eigenvectors
    [dmd_evecs,dmd_evals] = eig(A_tilde);
    
    % Compute DMD modes
    dmd_modes = Y*Vr/Sr*dmd_evecs;
    
    % Compute DMD amplitudes
    dmd_amplitudes = pinv(dmd_modes)*X(:,1);
end