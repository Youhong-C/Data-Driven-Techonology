%-----------------------------------------------
% Name of file : DMD_example.m
% 
% Created   : 04/04/2025
%
% Purpose   : Implementation of the DMD algorithm
%           
% Author    : Youhong Chen
%
% Copyright : Youhong Chen, Debraj Bhattacharjee, 2025
%------------------------------------------------

clc
clear
close

%% Define State-Space
n = 4;
sys = drss(n);
sys.D = 0;

if ~isstable(sys)
    error('Sampled system not stable');
end

nStep = 100;
x0 = randn(n,1);
xCurr = x0;

%% Setup for DMDc - Simulate linear system
BigX = zeros(n,nStep);
Gamma = randn(1,nStep-1);

BigX(:,1) = xCurr;
Z = zeros(nStep - 1, length(xCurr));
Z_orig = Z;

for i = 2:nStep
    u = Gamma(i-1);
    xNext = sys.A*xCurr + sys.B*u*0;

    BigX(:,i) = xNext;
    xCurr = xNext;
    Z(i,:) = xCurr';
end

%% Run DMD
% A_tilde = DMD(BigX);
% [A_tilde, dmd_modes, dmd_amplitudes, dmd_evals, dmd_evecs, U] = DMD(BigX);
[A_tilde, dmd_modes, dmd_amplitudes, dmd_evals, dmd_evecs, U] = DMD_Truncated(BigX);

%% Compare
disp('Eigenvalues of original system');
eig(sys.A)

disp('Eigenvalues of identified system');
eig(A_tilde)

%% Compare responses
xRoute1 = zeros(n, nStep);

xCurr = x0;
xRoute1(:,1) = xCurr;
xRoute2 = xRoute1;

% xRoute1
for i = 2:nStep
    xNext = A_tilde*xCurr;

    xRoute1(:,i) = xNext;
    xCurr = xNext;
    Z_orig(i,:) = xCurr';
end

xCurr = x0;
% xRoute2
for i = 2:nStep
    xNext = (dmd_evecs*dmd_evals/dmd_evecs)*xCurr;

    xRoute2(:,i) = xNext;
    xCurr = xNext;
end
xError = abs(xRoute1) - abs(xRoute2);

%% Verification

Y = BigX(:,2:end);
X = BigX(:,1:end-1);

errorMatrix = Y - A_tilde*X;

%% Find Transformation matrix
T = (Z_orig\Z)';

xReconstruct = zeros(n, nStep);

xCurr = x0;
xReconstruct(:,1) = xCurr;

% xRoute1
for i = 2:nStep
    xNext = A_tilde*xCurr;

    xReconstruct(:,i) = T*xNext;
    xCurr = xNext;    
end

% Plot
stateNumber = 1;

% plot(X(stateNumber,:))
% hold on
% plot(xReconstruct(stateNumber,:))


%% Eigenvector
[V_system,~] = eig(sys.A);

[V_DMD,~] = eig(A_tilde);

V_Transformed = T\V_system;

%% Participation Factor
A_Reconstructed=T*A_tilde*inv(T);

disp('Difference between the original A and reconstructed A');
sys.A-A_Reconstructed;
[V_Reconstructed,~] = eig(A_Reconstructed);

% Oringial PF
evr=V_system;
evl = inv(evr);

pf = abs(evr.*evl.');    % Participation factor matrix

pf_norm1 = pf./max(pf);

% DMD PF with T
evr=dmd_modes;
evl = inv(evr);

pf = abs(evr.*evl.');    % Participation factor matrix

pf_norm2 = pf./max(pf);

%% Reproduction




% Time step and duration
% dt = Tsd; % Time step (seconds)
% T = 1;     % Total time (seconds)
time = 1:1:nStep; % Time vector

% Preallocate for state history
x_project = zeros(n, length(time));
x_project(:, 1) = x0; % Initial state

Phi=U*dmd_evecs;
Lambda=dmd_evals;
% Simulate the system using projected eigenvector decomposition
for k = 2:nStep
    x_project(:, k) = Phi * (Lambda^(k-1)) * (Phi \ x0); % State update
end



% Verify equivalence with exact eigenvector simulation
x_exact = zeros(n, length(time));
x_exact(:, 1) = x0;
for k = 2:nStep
    x_exact(:, k) = dmd_modes * (Lambda^(k-1)) * (dmd_modes \ x0); % State update
end

% Plot results
figure;
hold on;
% Projected data
plot(time, x_project(1, 1:nStep), '--', 'Color', [1, 0, 0, 0.2], 'LineWidth', 5, 'DisplayName', 'x_1(t) (projected)');
plot(time, x_project(2, 1:nStep), '--', 'Color', [0, 0, 1, 0.2], 'LineWidth', 5, 'DisplayName', 'x_2(t) (projected)');

% Exact data
plot(time, x_exact(1, 1:nStep), ':', 'Color', [1, 0, 1, 0.8], 'LineWidth', 3, 'DisplayName', 'x_1(t) (exact)');
plot(time, x_exact(2, 1:nStep), ':', 'Color', [0, 1, 0, 0.8], 'LineWidth', 3, 'DisplayName', 'x_2(t) (exact)');

% Original data
data = Z';
plot(time, data(1, 1:nStep), '-', 'Color', [1, 0, 0, 1], 'LineWidth', 1, 'DisplayName', 'x_1(t) (original)');
plot(time, data(2, 1:nStep), '-', 'Color', [0, 0, 1, 1], 'LineWidth', 1, 'DisplayName', 'x_2(t) (original)');

hold off;
grid on;
xlabel('Time (s)');
ylabel('State Variables');
legend show;
title('Comparison of State Evolution: exact vs dmd');