function [A_tilde, dmd_modes, dmd_amplitudes, dmd_evals, dmd_evecs, U] = DMD_Truncated(BigX)

%-----------------------------------------------
% Name of file : DMD_Truncated.m
% 
% Created   : 24/01/2025
%
% Purpose   : Implementation of the DMD algorithm
%             with Truncated Singular Values
%           
% Author    : Debraj Bhattacharjee
%
% Copyright : Debraj Bhattacharjee, 2025
%------------------------------------------------

    % Define X
    X = BigX(:,1:end-1);

    % Define Y (or X' as in DMD paper)
    Y = BigX(:,2:end);  

    [A_tilde, dmd_modes, dmd_amplitudes, dmd_evals, dmd_evecs, U] = DMD_Truncated_Core(X,Y);
end