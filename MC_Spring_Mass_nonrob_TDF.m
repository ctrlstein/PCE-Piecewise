clear all; clc; close all;


%%%%%%%%%%%%%%%%%%%%
%%% Outer System %%%
%%%%%%%%%%%%%%%%%%%%

% samples_vec = [1000 5000 10000 50000 100000];
samples_vec = [10000];
omgn_m = pi;
omgn_d_1 = 0.25*pi;
omgn_d_2 = 0.5*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nonrobust input shaper %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TDF_vec = [0.5 0.5 1];
u = 1;

%%%%%%%%%%%%%%%%%%%%%
%%% MC Simulation %%%
%%%%%%%%%%%%%%%%%%%%%

for size_i = 1:length(samples_vec)
    size_i
    omgn_vec_1 = omgn_m + omgn_d_1*(2*rand(samples_vec(size_i), 1) - 1);
    omgn_vec_2 = omgn_m + omgn_d_2*(2*rand(samples_vec(size_i), 1) - 1);

    tic
    res_temp = fct_cost_TDF_MC(omgn_vec_1, omgn_vec_2, samples_vec(size_i), TDF_vec, u);
    time_vec(:, size_i) = [toc samples_vec(size_i)];
    fval_opt_vec(:, size_i) = [res_temp{1} samples_vec(size_i)];
    fval_opt_var_vec(:, size_i) = [res_temp{5} samples_vec(size_i)];
    
end


%%%%%%%%%%%%%%%%
%%% Analysis %%%
%%%%%%%%%%%%%%%%

figure;
bar(fval_opt_vec(1, :))
xticks(linspace(1, length(samples_vec), length(samples_vec)))
xticklabels(string(fval_opt_vec(2, :)))
xlabel('#Samples');
ylabel('E[V_R] in J');
title('Residual Energy (Mean) per Sample Size');

figure;
bar(fval_opt_var_vec(1, :))
xticks(linspace(1, length(samples_vec), length(samples_vec)))
xticklabels(string(fval_opt_var_vec(2, :)))
xlabel('#Samples');
ylabel('Var(V_R) in J^2');
title('Residual Energy (Variance) per Sample Size');

figure;
bar(time_vec(1, :))
xticks(linspace(1, length(samples_vec), length(samples_vec)))
xticklabels(string(time_vec(2, :)))
xlabel('#Samples');
ylabel('Seconds');
title('Time needed per Sample Size');


%%%%%%%%%%%%%%%%%
%%% Functions %%%
%%%%%%%%%%%%%%%%%

function dxdt = fct_inside_real_system(t,x,omgn,x_TDF,u)
    u_mod = (x_TDF(1) + x_TDF(2)*heaviside(t-x_TDF(3))) * u;
    dxdt(1,1) = x(2);
    dxdt(2,1) = -omgn^2*x(1) + omgn^2*u_mod;
end

function J = fct_cost_TDF_MC(fct_omgn_vec_1, fct_omgn_vec_2, N, TDF_MC_vec, u)
% Function that computes solves spring-mass differential equation and
% returns x, x_dot as well as expected residual energy and its variance at
% t_2 for given pair (omega_n, omega_m)
    t1_vec = linspace(0,10,101);
    t2_vec = linspace(10,20,101);
    x0_MC_1 = [0, 0];
    
    t = [t1_vec t2_vec];
    y_1_t_vals_mat = zeros(length(t), N);
    y_2_t_vals_mat = zeros(length(t), N);
    
    options = odeset('RelTol', 1e-12, 'AbsTol',1e-12);

    parfor sample_i = 1:N
        sample_i
        [TT1, YY_MC_1] = ode45(@(t,x) fct_inside_real_system(t,x,fct_omgn_vec_1(sample_i),TDF_MC_vec,u), t1_vec, x0_MC_1, options);
        x0_MC_2 = YY_MC_1(end, :);
        [TT2, YY_MC_2] = ode45(@(t,x) fct_inside_real_system(t,x,fct_omgn_vec_2(sample_i),TDF_MC_vec,u), t2_vec, x0_MC_2, options);
    
        y_1_t_vals_mat(:, sample_i) = [YY_MC_1(:, 1); YY_MC_2(:, 1)];
        y_2_t_vals_mat(:, sample_i) = [YY_MC_1(:, 2); YY_MC_2(:, 2)];
    
        V_res_TDF_MC_vec(sample_i) = 0.5*(YY_MC_2(end,2))^2 + fct_omgn_vec_2(sample_i)^2 * 0.5 * (1 - YY_MC_2(end,1))^2;
    
    end
    
    ve = mean(V_res_TDF_MC_vec(:));
    vv = var(V_res_TDF_MC_vec(:));
    
    J = {ve, t, y_1_t_vals_mat, y_2_t_vals_mat, vv};
end
