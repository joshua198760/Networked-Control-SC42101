%% Preliminaries
clear
clc 
ct.A=[-3.7 -7.5; 0 1];ct.B=[0;1];%System matrices
ct.nx=length(ct.A);ct.nu=size(ct.B,2);
if rank(ctrb(ct.A,ct.B))==length(ct.A)%Check controllability
    disp('CT System Controllable')
end
[ct.V,ct.D]=eig(ct.A);

%Closed Loop
ct.clpoles=[-1+2i,-1-2i];%Desired closed loop poles
ct.K=place(ct.A,ct.B,ct.clpoles);%Static Feedback Gain

addpath(genpath("C:\Users\amato\OneDrive\Documenten\TU\MSC S&C\Year 1\Q3 23-24\MPC"))


%% Question 1
q1.h1min = 0.005;q1.h1max = 0.7;%range of h1
numPoints = 150;q1.h1range=linspace(q1.h1min,q1.h1max,numPoints);
plstrategy=['z','h'];%packet loss strategy
q1.sys1stab=zeros(length(q1.h1range),length(plstrategy));

for s=1:length(plstrategy)
    for i=1:length(q1.sys1stab)
        h1 = q1.h1range(i);
        [F1,F0]=cttodt_pl(ct.A,ct.B,h1,plstrategy(s),ct.K); 
        nx=size(F1,1);
        P=sdpvar(nx,nx);Q=1e-3*eye(nx,nx);
        cons=[P>=1e-3*eye(nx),...
              (F0^3)'*P*(F0^3)-P<=-Q,...
              (F0^2*F1)'*P*(F0^2*F1)-P<=-Q];
        obj=0;ops=sdpsettings('verbose',0,'solver','sdpt3');
        sol=optimize(cons,obj,ops);%Solve
         if sol.problem==0
             q1.sys1stab(i,s)=1;%1 means stable
         else
             q1.sys1tab(i,s)=0;%0 means unstable
         end
        clear P
    end
end


% Plotting the results
figure;
hold on;
for s = 1:length(plstrategy)
    % Get stable and unstable points
    stable_points = q1.h1range(q1.sys1stab(:, s) == 1);
    unstable_points = q1.h1range(q1.sys1stab(:, s) == 0);
    
    % Plot stable points
    scatter(stable_points, s * ones(size(stable_points)), 'g', 'filled');
    % Plot unstable points
    scatter(unstable_points, s * ones(size(unstable_points)), 'r', 'filled');
end

% Customize the plot
legend({'Stable', 'Unstable'}, 'Location', 'best');
xlabel('h1');
ylabel('Packet Loss Strategy');
yticks(1:length(plstrategy));
yticklabels({'z', 'h'});
title('Stability of System for Different Packet Loss Strategies and h1 values');
grid on;
hold off;


%% Question 2
% Define the probabilities
p1 = 0.99;p2 = 0.49;p3 = 0.02;

% Create the matrix T
T = [0   0   p1  0   0   0   1-p1;
     0   0   p1  0   0   0   1-p1;
     0   0   0   p2  p2  p3  0;
     p1 1-p1 0   0   0   0   0;
     p1 1-p1 0   0   0   0   0;
     p1 1-p1 0   0   0   0   0;
     0   0   0   p2  p2  p3  0];

%Define the CT system matrices of Sys2
ct.A2=(1/3)*ct.A;ct.B2=ct.B;
ct.K2=place(ct.A2,ct.B2,ct.clpoles);%Static Feedback Gain

%Range of h1
h_values = 0.001:0.005:0.2;
%Store Results
stable_points = [];
unstable_points = [];

for h1 = h_values
    %Sys Matrices (1,3) and (2,7) are equal
    Acl1 = mjls(h1, ct.K, ct.K2, 1);Acl2 = mjls(h1, ct.K, ct.K2, 2);
    Acl4 = mjls(h1, ct.K, ct.K2, 4);Acl5 = mjls(h1, ct.K, ct.K2, 5);
    Acl6 = mjls(h1, ct.K, ct.K2, 6);
    Acl3 = Acl1; Acl7 = Acl2;
    
    nx = length(Acl1);eps = 1e-3; % machine precision
    S = eye(length(T)); % Selection Matrix
    
    %Opti Variables
    P1 = sdpvar(nx, nx); P2 = sdpvar(nx, nx);
    P3 = sdpvar(nx, nx); P4 = sdpvar(nx, nx);
    P5 = sdpvar(nx, nx); P6 = sdpvar(nx, nx);
    P7 = sdpvar(nx, nx);

    %LMIs
    cons = [P1 - ((T(1,:) * S(:,3)) * (Acl3' * P3 * Acl3) ...
           + T(1,:) * S(:,7) * (Acl7' * P7 * Acl7)) >= eps * eye(nx), ...
           
          P2 - ((T(2,:) * S(:,3)) * (Acl3' * P3 * Acl3) ...
          + T(2,:) * S(:,7) * (Acl7' * P7 * Acl7)) >= eps * eye(nx), ...
          
          P3 - ((T(3,:) * S(:,4)) * (Acl4' * P4 * Acl4) ...
          + T(3,:) * S(:,5) * (Acl5' * P5 * Acl5) ...
          + T(3,:) * S(:,6) * (Acl6' * P6 * Acl6)) >= eps * eye(nx), ...
          
          P4 - ((T(4,:) * S(:,1)) * (Acl1' * P1 * Acl1) ...
          + T(3,:) * S(:,2) * (Acl2' * P2 * Acl2)) >= eps * eye(nx), ...
          
          P5 - ((T(5,:) * S(:,1)) * (Acl1' * P1 * Acl1) ...
          + T(5,:) * S(:,2) * (Acl2' * P2 * Acl2)) >= eps * eye(nx), ...
          
          P6 - ((T(6,:) * S(:,1)) * (Acl1' * P1 * Acl1) ...
          + T(6,:) * S(:,2) * (Acl2' * P2 * Acl2)) >= eps * eye(nx), ...
          
          P7 - ((T(7,:) * S(:,4)) * (Acl4' * P4 * Acl4) ...
          + T(7,:) * S(:,5) * (Acl5' * P5 * Acl5) ...
          + T(7,:) * S(:,6) * (Acl6' * P6 * Acl6)) >= eps * eye(nx), ...
          
          P1 >= eps * eye(nx), P2 >= eps * eye(nx), P3 >= eps * eye(nx), ...
          P4 >= eps * eye(nx), P5 >= eps * eye(nx), P6 >= eps * eye(nx), ...
          P7 >= eps * eye(nx)];

    %Optimize
    obj = 0; 
    ops = sdpsettings('verbose', 0, 'solver', 'sdpt3', 'debug', 1);
    sol = optimize(cons, obj, ops); % Solve
    
    %Save sampling time as Stable or Unstable
    if sol.problem == 0
        stable_points = [stable_points, h1];
    else
        unstable_points = [unstable_points, h1];
    end
end

% Plot the results
figure;
hold on;
plot(stable_points, zeros(size(stable_points)), 'go', 'MarkerFaceColor', 'g');
plot(unstable_points, zeros(size(unstable_points)), 'ro', 'MarkerFaceColor', 'r');
xlabel({'$h_1$'},'interpreter','latex');
title({'MSS of $h_1$ (to-hold)'},'interpreter','latex');
legend({'Stable', 'Unconcluded'},'Interpreter','latex');
hold off;

%% Q3 (compare with Q1)
% First Analysis: Polytopic Approach
hrange = linspace(0.005,0.65,30);
stable_points = [];
unstable_points = [];

for h = hrange
    taurange = 0:0.02:h-0.02;
    for tau = taurange
        [Hf, Hg] = Hsecond(h,tau,ct.A,ct.V);
        K = [ct.K 0];
        H1 = Hf.mi_mi - Hg.mi_mi * K;
        H2 = Hf.ma_ma - Hg.ma_ma * K;
        H3 = Hf.comb  - Hg.comb * K;
        nx = length(H1);
        
        P = sdpvar(nx, nx);
        gam=1e-5;
        cons = [P >= 1e-3 * eye(nx, nx), ...
                H1' * P * H1 - (1-gam)*P <= 0, ...
                H2' * P * H2 - (1-gam)*P <= 0, ...
                H3' * P * H3 - (1-gam)*P <= 0];
        
        obj = 0;
        ops = sdpsettings('verbose', 0,'solver','sdpt3',...
            'sdpt3.inftol',1e-5,'sdpt3.gaptol',1e-5,...
            'sdpt3.inftol',1e-5);
        sol = optimize(cons, obj, ops);
        
        if sol.problem == 0 
            stable_points = [stable_points; h, tau];
        else 
            unstable_points = [unstable_points; h, tau];
        end
    end
end

% Finding the points for interpolation
[max_tau_range, max_tau_idx] = max(stable_points(:,2) - 0);
intermediate_point = stable_points(max_tau_idx, :);

% Finding the point with the highest h that is stable
[~, max_h_idx] = max(stable_points(:, 1));
highest_h_point = stable_points(max_h_idx, :);

% Define the x and y coordinates for the dividing line
x = [0, intermediate_point(1), highest_h_point(1)];
y = [0, intermediate_point(2), highest_h_point(2)];

% Second Analysis: Eigenvalue Approach
exp2.hmin = 0.005;
exp2.hmax = 0.65;
exp2.hrange = linspace(exp2.hmin, exp2.hmax, 30); % range of inter-sampling times
stable_points_eigen = [];
unstable_points_eigen = [];

for i = 1:length(exp2.hrange)
    taumax = exp2.hrange(i) - 0.02;
    for tau = 0:0.02:taumax
        [F, G] = ctsstodtssdel(ct.A, ct.B, exp2.hrange(i), tau);
        cloopmat = F - G * [ct.K zeros(1, ct.nu)]; % DT closed-loop matrix
        specrad = max(abs(eig(cloopmat)));
        if specrad < 1
            stable_points_eigen = [stable_points_eigen; exp2.hrange(i), tau];
        else
            unstable_points_eigen = [unstable_points_eigen; exp2.hrange(i), tau];
        end
    end
end

%% Plot Comparison
% Finding the points for interpolation
[max_tau_range_eigen, max_tau_idx_eigen] = max(stable_points_eigen(:,2) - 0);
intermediate_point_eigen = stable_points_eigen(max_tau_idx_eigen, :);

% Finding the point with the highest h that is stable
[~, max_h_idx_eigen] = max(stable_points_eigen(:, 1));
highest_h_point_eigen = stable_points_eigen(max_h_idx_eigen, :);

% Define the x and y coordinates for the dividing line
x_eigen = [0, intermediate_point_eigen(1), highest_h_point_eigen(1)];
y_eigen = [0, intermediate_point_eigen(2), highest_h_point_eigen(2)];

figure;

% Subplot 1: Polytopic Approach
subplot(2, 1, 1);
hold on;
title('Polytopic Approach');
fill([0, x, hrange(end)], [0, y, 0], 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([x(2), x(3), hrange(end), hrange(end)], [y(2), 0, 0, hrange(end) - 0.02], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(x, y, 'k', 'LineWidth', 1.5);
ylabel('$\tau$','interpreter','latex');
legend({'Stable Region', 'Unstable Region'}, 'Location', 'Best');
grid on;
hold off;

% Subplot 2: Eigenvalue Approach
subplot(2, 1, 2);
hold on;
title('Eigenvalue Approach');
fill([0, x_eigen, exp2.hrange(end)], [0, y_eigen, 0], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([x_eigen(2), x_eigen(3), exp2.hrange(end), exp2.hrange(end)], [y_eigen(2), 0, 0, exp2.hrange(end) - 0.02], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(x_eigen, y_eigen, 'k', 'LineWidth', 1.5);

xlabel('$h$','Interpreter','latex');
ylabel('$\tau$','interpreter','latex');
legend({'Stable Region', 'Unstable Region'}, 'Location', 'Best');
grid on;
hold off;


%% Visualization of Approximation
% Define parameters
h = 0.2;
taurange = linspace(0, h-0.01, 100);
F13 = [];F23 = [];
G11 = [];G21 = [];

for tau = taurange
    [bigF, bigG] = BigFG(h, tau, ct.A, ct.V);
    F13 = [F13; bigF(1, 3)];
    F23 = [F23; bigF(2, 3)];
    G11 = [G11; bigG(1, 1)];
    G21 = [G21; bigG(2, 1)];
end

%Obtain Vertices
taumax=taurange(end);
[vertF,vertG]=BigFGflex(h,taumax,ct.A,ct.V);

% Plot Uncertain matrix entries in F
figure;
plot(F13, F23, 'b');
hold on;

% Add vertices of polytopic approximation
extreme_F = [vertF.mi_mi(1,3), vertF.mi_mi(2,3);
             vertF.ma_ma(1,3), vertF.ma_ma(2,3);
             vertF.ma_ma(1,3), vertF.mi_mi(2,3)]; 
plot(extreme_F(:,1), extreme_F(:,2), 'ro', 'MarkerFaceColor', 'r');
fill(extreme_F(:,1), extreme_F(:,2), 'r', 'FaceAlpha', 0.1);
hold off;
title('Uncertain matrix entries in $F$', 'Interpreter', 'latex');
xlabel('$F_{13}$', 'Interpreter', 'latex');
ylabel('$F_{23}$', 'Interpreter', 'latex');
legend({'$F_{13}-F_{23}$','$\mathcal{H}_F$','$\overline{\mathcal{F}}$'},'Interpreter','latex')

% Plot Uncertain matrix entries in G
figure;
plot(G11, G21, 'b');
hold on;
% Add vertices of polytopic approximation
extreme_G = [vertG.mi_mi(1,1), vertG.mi_mi(2,1);
             vertG.ma_ma(1,1), vertG.ma_ma(2,1);
             vertG.ma_ma(1,1), vertG.mi_mi(2,1)];  
plot(extreme_G(:,1), extreme_G(:,2), 'ro', 'MarkerFaceColor', 'r');
fill(extreme_G(:,1), extreme_G(:,2), 'r', 'FaceAlpha', 0.1);
hold off;
title('Uncertain matrix entries in $G$', 'Interpreter', 'latex');
xlabel('$G_{11}$', 'Interpreter', 'latex');
ylabel('$G_{21}$', 'Interpreter', 'latex');
legend({'$G_{11}-G_{21}$','$\mathcal{H}_G$','$\overline{\mathcal{G}}$'},'Interpreter','latex')

%% Visualization Refined Polytope
% Define parameters
h = 0.2;
taurange = linspace(0, h-0.01, 100);
F13 = [];F23 = [];
G11 = [];G21 = [];

for tau = taurange
    [bigF, bigG] = BigFG(h, tau, ct.A, ct.V);
    F13 = [F13; bigF(1, 3)];
    F23 = [F23; bigF(2, 3)];
    G11 = [G11; bigG(1, 1)];
    G21 = [G21; bigG(2, 1)];
end

%Obtain Vertices
taumax=taurange(end);
[vertF,vertG]=BigFGflex(h,taumax,ct.A,ct.V,...
                        0.2,0.45,0.05,0.2);

% Plot Uncertain matrix entries in F
figure;
plot(F13, F23, 'b');
hold on;

% Add vertices of polytopic approximation
extreme_F = [vertF.mi_mi(1,3), vertF.mi_mi(2,3);
             vertF.ma_ma(1,3), vertF.ma_ma(2,3);
             vertF.comb(1,3), vertF.comb(2,3);
             vertF.comb2(1,3), vertF.comb2(2,3);]; 
plot(extreme_F(:,1), extreme_F(:,2), 'ro', 'MarkerFaceColor', 'r');
fill(extreme_F(:,1), extreme_F(:,2), 'r', 'FaceAlpha', 0.1);
hold off;
title('Uncertain matrix entries in $F$', 'Interpreter', 'latex');
xlabel('$F_{13}$', 'Interpreter', 'latex');
ylabel('$F_{23}$', 'Interpreter', 'latex');
legend({'$F_{13}-F_{23}$','$\mathcal{H}_F$','$\overline{\mathcal{F}}$'},'Interpreter','latex')

% Plot Uncertain matrix entries in G
figure;
plot(G11, G21, 'b');
hold on;
% Add vertices of polytopic approximation
extreme_G = [vertG.mi_mi(1,1), vertG.mi_mi(2,1);
             vertG.ma_ma(1,1), vertG.ma_ma(2,1);
             vertG.comb(1,1), vertG.comb(2,1);
             vertG.comb2(1,1), vertG.comb2(2,1);];  
plot(extreme_G(:,1), extreme_G(:,2), 'ro', 'MarkerFaceColor', 'r');
fill(extreme_G(:,1), extreme_G(:,2), 'r', 'FaceAlpha', 0.1);
hold off;
title('Uncertain matrix entries in $G$', 'Interpreter', 'latex');
xlabel('$G_{11}$', 'Interpreter', 'latex');
ylabel('$G_{21}$', 'Interpreter', 'latex');
legend({'$G_{11}-G_{21}$','$\mathcal{H}_G$','$\overline{\mathcal{G}}$'},'Interpreter','latex')

hrange = linspace(0.005,0.65,15);
stable_points = [];
unstable_points = [];

%% Analysis Refined Polytope
for h = hrange
    taurange = 0:0.02:h-0.02;
    for tau = taurange
        [Hf, Hg] = Hrefined(h,tau,ct.A,ct.V,0.2,0.45,0.05,0.2);
        K = [ct.K 0];
        H1 = Hf.mi_mi - Hg.mi_mi * K;
        H2 = Hf.ma_ma - Hg.ma_ma * K;
        H3 = Hf.comb  - Hg.comb * K;
        H4 = Hf.comb2  - Hg.comb2 * K;
        nx = length(H1);
        
        P = sdpvar(nx, nx);
        gam=1e-7;
        cons = [P >= 1e-3 * eye(nx, nx), ...
                H1' * P * H1 - (1-gam)*P <= 0, ...
                H2' * P * H2 - (1-gam)*P <= 0, ...
                H3' * P * H3 - (1-gam)*P <= 0,...
                H4' * P * H4 - (1-gam)*P <= 0];
        
        obj = 0;
        ops = sdpsettings('verbose', 0,'solver','sdpt3',...
            'sdpt3.inftol',1e-5,'sdpt3.gaptol',1e-5,...
            'sdpt3.inftol',1e-5);
        sol = optimize(cons, obj, ops);
        
        if sol.problem == 0 
            stable_points = [stable_points; h, tau];
        else 
            unstable_points = [unstable_points; h, tau];
        end
    end
end
% Finding the points for interpolation
[max_tau_range, max_tau_idx] = max(stable_points(:,2) - 0);
intermediate_point = stable_points(max_tau_idx, :);

% Finding the point with the highest h that is stable
[~, max_h_idx] = max(stable_points(:, 1));
highest_h_point = stable_points(max_h_idx, :);

% Define the x and y coordinates for the dividing line
x = [0, intermediate_point(1), highest_h_point(1)];
y = [0, intermediate_point(2), highest_h_point(2)];

% Finding the points for interpolation LOAD THE EIGENVALUE ANALYSIS ABOVE
% FIRST
[max_tau_range_eigen, max_tau_idx_eigen] = max(stable_points_eigen(:,2) - 0);
intermediate_point_eigen = stable_points_eigen(max_tau_idx_eigen, :);

% Finding the point with the highest h that is stable
[~, max_h_idx_eigen] = max(stable_points_eigen(:, 1));
highest_h_point_eigen = stable_points_eigen(max_h_idx_eigen, :);

% Define the x and y coordinates for the dividing line
x_eigen = [0, intermediate_point_eigen(1), highest_h_point_eigen(1)];
y_eigen = [0, intermediate_point_eigen(2), highest_h_point_eigen(2)];

% Create figure with two subplots
figure;

% Subplot 1: Polytopic Approach
subplot(2, 1, 1);
hold on;
title('Polytopic Approach');
fill([0, x, hrange(end)], [0, y, 0], 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([x(2), x(3), hrange(end), hrange(end)], [y(2), 0, 0, hrange(end) - 0.02], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(x, y, 'k', 'LineWidth', 1.5);
ylabel('$\tau$','interpreter','latex');
legend({'Stable Region', 'Unstable Region'}, 'Location', 'Best');
grid on;
hold off;

% Subplot 2: Eigenvalue Approach
subplot(2, 1, 2);
hold on;
title('Eigenvalue Approach');
fill([0, x_eigen, exp2.hrange(end)], [0, y_eigen, 0], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([x_eigen(2), x_eigen(3), exp2.hrange(end), exp2.hrange(end)], [y_eigen(2), 0, 0, exp2.hrange(end) - 0.02], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(x_eigen, y_eigen, 'k', 'LineWidth', 1.5);
xlabel('$h$','Interpreter','latex');
ylabel('$\tau$','interpreter','latex');
legend({'Stable Region', 'Unstable Region'}, 'Location', 'Best');
grid on;
hold off;
%% Question 4
% Define the state-space matrices
A = [ct.A]; B = [ct.B];
K = ct.K; Q = eye(2);
P = lyap((A - B * K), Q);
% Initial conditions
initial_conditions = {[3; 0]};
sigma_values = [0.1, 0.5, 0.8];
t0 = 0; tend = 15;  % Simulation time span

figure;

for i = 1:length(sigma_values)
    sigma = sigma_values(i);

    subplot(3, 1, i);
    hold on;

    for j = 1:length(initial_conditions)
        x0 = initial_conditions{j};
        sk = 0;  % Initial sampling time
        xi_sk = x0;
        u = -K * xi_sk;  % Initial Input
        tsim = [t0, tend];

        % Simulation
        T = []; Xi = [];
        eT = []; eS = [];

        while tsim(1) < tsim(2)
            % Simulate until trigger
            options = odeset('Events', @(t, xi) event_func(t, xi, xi_sk, sigma));
            [Tp, Xp, TE, YE, IE] = ode45(@(t, xi) update_sysdyn(t, xi, u), tsim, x0, options);

            % Store times and states
            T = [T; Tp];
            Xi = [Xi; Xp];

            if isempty(TE) % no event triggered
                break;
            end

            % Update initial condition and time
            tsim(1) = TE;
            x0 = YE';

            % Update the last sampling time
            sk = TE;
            xi_sk = YE';
            u = -K * xi_sk;

            % Store event time and state
            eT = [eT; TE];
            eS = [eS; YE];
        end

        % Plot results
        plot(T, Xi);
        plot(eT, eS, 'r*', 'MarkerSize', 5);
    end

    xlabel('Time');
    ylabel('Amplitude');
    title(['Sigma = ', num2str(sigma)]);
    legend('$x_1$', '$x_2$', 'Trigger', 'interpreter', 'latex');
    grid on;
    hold off;
end

function dxi = update_sysdyn(t, xi, u)
    A = [-3.7 -7.5; 0 1];
    B = [0; 1];
    dxi = A * xi + B * u;
end

function [phi_val, isterminal, direction] = event_func(t, xi, xk, sigma)
    A = [-3.7 -7.5; 0 1];
    B = [0; 1];
    poles = [-1+2i, -1-2i];
    K = place(A, B, poles);
    Q = eye(2);
    P = lyap((A-B*K)', Q);
    %Check triggering condition
    phi_mat = [A'*P + P*A + sigma*Q, -P*B*K; -(B*K)'*P, zeros(size(K, 2))];
    z = [xi; xk];
    phi_val = z' * phi_mat * z;
    isterminal = 1;  % Stop the integration
    direction = 0;  % Detect all zero crossings
end