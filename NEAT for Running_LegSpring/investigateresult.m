% Friedl De Groote
% May 26, 2015

close all
clear all

global time_fb p_COM v_COM a_COM e_con inputdata perturbation

setinputdata

% Load result
load neatsave

% Find individual with maximal fitness
fitness = 0;
index_max = 0;
for i = 1:length(population)
    population(i).fitness
    if population(i).fitness > fitness
        fitness = population(i).fitness;
        index_max = i;
    end
end

color_palette= [255 0 0; 255 150 150; 0 114 54; 153 255 102; 0 0 255; 44 176 207; 0 0 0; 155 155 155; 255 153 0; 255 185 0; 153 51 102; 171 218 77]*(1/255);
color=color_palette([6 1 2 10 11 4 3 5 12 7 8 9 1 10 4 6 11 3 2 5 12],:); %S0052 functional vars (EMGGRFandCOM file)

perturbation = 2;   % Backward perturbation: 1.
                    % Forward perturbation: 2.

figure(1)
figure(2)
for k = 1:10
    %% Forward simulation of posture after perturbation
    % model state x = [a1 a2 a3 lm1 lm2 lm3 theta theta_dot]
    % Initial condition
    % Initial ankle angle from uniform distribution between theta_0_min and
    % theta_0_max
%     q_0_min = -0.1;
%     q_0_max = 0.1;
    q_0_min = -0.05;
    q_0_max = 0.05;
    q_0 = rand(1)*(q_0_max - q_0_min) - 0.5 * (q_0_max - q_0_min)
    % => static posture: theta_dot_0, theta_ddot_0 = 0
    qdot_0 = 0;
    qddot_0 = 0;
    % Compute corresponding inverse dynamic torque
    theta_ref = 0.75;
    aPLAT = interp1(inputdata.atime, inputdata.aPLAT, 0);
    Tlimit = - exp(20*(q_0-theta_ref)) + exp(20*(-q_0-theta_ref));
    TID = qddot_0 * inputdata.m - inputdata.m * 9.81 * inputdata.h * sin(q_0) - inputdata.m * aPLAT * inputdata.h * cos(q_0) - Tlimit;
    % Solve initial muscle force sharing problem
    Tmus = zeros(3,1);
    % Co-contraction from uniform distribution between c_min and c_max (in
    % percentage)
    c_max = 0.2;
    c = rand(1)*(c_max - 0.05)+0.05;
    % Determine plantar flexor and dorsi flexor torque
    if TID < 0 % SOL, GAS
        T_SG = (1+c) * TID;
        Tmus(3) = -c * TID; % TA
    else
        Tmus(3) = (1+c) * TID;
        T_SG = -c * TID;
    end
    % Divide plantar flexor torque over SOL and GAS
    F_0_GAS = inputdata.params(1,1);
    F_0_SOL = inputdata.params(1,2);
    F_0_SG = F_0_GAS + F_0_SOL;
    % Specify relative contribution of SOL and GAS
    delta_min = -0.2;
    delta_max = 0.2;
    delta = rand(1)*(delta_max - delta_min) - 0.5 * (delta_max - delta_min);
    Tmus(1) = (F_0_GAS/F_0_SG - delta) * T_SG;
    Tmus(2) = (F_0_SOL/F_0_SG + delta) * T_SG;
    % Muscles isometric at onset perturbation
    % => Pick lM0 such that vM0 = 0
%     load muscle_splines
    muscle_ind = [7 8 9]; % GAS SOL TA
    
    lM0 = zeros(inputdata.M,1);
    a0 = zeros(inputdata.M,1);
    for i = 1:inputdata.M
        ma = ppval(inputdata.q_ankle(muscle_ind(i)).dM_ankle, q_0);
        lMT0 = ppval(q_ankle(muscle_ind(i)).dM_Length, q_0);
        fse = Tmus(i)/(ma*inputdata.params(1,i));
        lTtilda = log(5*(fse + 0.25))/35 + 0.995;
        lT = lTtilda * inputdata.params(3,i);
        lM0(i) = (lMT0 - lT)./inputdata.params(2,i);
        a = fsolve(@(a) ForceEquilibriumFTstate(a,fse,0,lMT0,0,inputdata.params(:,i), inputdata.Fvparam, inputdata.Fpparam, inputdata.Faparam), 0.02);
        % ForceEquilibrium(a0(i),x,0,lMT0,inputdata.params(:,i), inputdata.Fvparam, inputdata.Fpparam, inputdata.Faparam)
        a(a < 0.01) = 0.01;
        a0(i) = a;
    end
    
    
    % Initial condition
    x0 = [a0; lM0; q_0; qdot_0]; % upright position
    
    t0 = 0;
    tf = 1;
    
    % xdot = OneSegmentThreeMuscles(t0, x0, inputdata)
    
    % Set integrator absolute and relative integrator error tolerances
    AbsTol = 1e-5;
    RelTol = 1e-3;
    options = odeset('AbsTol',AbsTol,'RelTol',RelTol);
    
    p_COM_0 = inputdata.h * sin(x0(7));
    v_COM_0 = 0; % assume for the moment that the initial velocity is zero
    a_COM_0 = 0; % assume for the moment that the initial acceleration is zero
    
    time_fb = [t0-tau_fb; t0];
    p_COM = [p_COM_0; p_COM_0];
    v_COM = [v_COM_0; v_COM_0];
    a_COM = [a_COM_0; a_COM_0];
    e_con = [zeros(1,3);zeros(1,3)];
    
    tspan = [t0: (tf-t0)/100 :tf];
    
    [t,x] = ode15s(@OneSegmentThreeMusclesController,tspan,x0, options, population(index_max));
    
    xCOM = inputdata.h * sin(x(:,7));
    
    % Test
    figure(1)
    subplot(3,1,1)
    plot(time_fb(2:end), p_COM(2:end), 'Color', color(k,:), 'LineWidth', 2); hold on;
    title('COM position')
    xlabel('t [s]')
    ylabel('p [m]')
    
    subplot(3,1,2)
    plot(time_fb(2:end), v_COM(2:end), 'Color', color(k,:), 'LineWidth', 2); hold on;
    title('COM velocity')
    xlabel('t [s]')
    ylabel('v [m/s]')
    
    subplot(3,1,3)
    plot(time_fb(2:end), a_COM(2:end), 'Color', color(k,:), 'LineWidth', 2); hold on;
    title('COM acceleration')
    xlabel('t [s]')
    ylabel('a [m/s]')
    
    figure(2)
    subplot(3,1,1)
    plot(tspan, x(:, 1), 'Color', color(k,:), 'LineWidth', 2); hold on;
    title('GAS')
    xlabel('t [s]')
    ylabel('a []')
    
    subplot(3,1,2)
    plot(tspan, x(:, 2), 'Color', color(k,:), 'LineWidth', 2); hold on;
    title('SOL')
    xlabel('t [s]')
    ylabel('a []')
    
    subplot(3,1,3)
    plot(tspan, x(:, 3), 'Color', color(k,:), 'LineWidth', 2); hold on;
    title('TA')
    xlabel('t [s]')
    ylabel('a []')
    
    e = interp1(time_fb, e_con, tspan);
    %%
    
end
