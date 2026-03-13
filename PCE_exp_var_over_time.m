clear all; clc; close all;
%%%%%%%%%%%%%%%%%%%%
%%% Outer System %%%
%%%%%%%%%%%%%%%%%%%%

omgn_m = pi;
omgn_d_1 = 0.25*pi;
omgn_d_2 = 0.5*pi;
pdf_1 = 0.5;
pdf_2 = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Robust Input shaper %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb_ris = [0 0 0 0 0];
ub_ris = [1 1 1 2 5];
x0_ris = [0.3 0.3 0.3 1 2];
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
[x_opt, fval_opt] = fmincon(@(x)fct_cost_TDF_rob(x,omgn_m), x0_ris, [], [], [], [], lb_ris, ub_ris, @(x)fct_noncon_TDF_rob(x,omgn_m), options);

A0_rob = x_opt(1); A1_rob = x_opt(2); A2_rob = x_opt(3); T1_rob = x_opt(4); T2_rob = x_opt(5);

TDF_rob_vec = [0.25, 0.5, 0.25, 1, 2];
TDF_vec = TDF_rob_vec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Robust input shaper: GSA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise weights and nodes for gauss legendre
N = 100;
gauss_legendre_coll = fct_gauss_legendre_wn(N);

degrees_vec = [10 20 30];

means_a_mat = zeros(3, 202);
means_b_mat = zeros(3, 202);
var_a_mat = zeros(3, 202);
var_b_mat = zeros(3, 202);

 for deg_i=1:length(degrees_vec)
    p = degrees_vec(deg_i)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PCE: Basis Functions %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    poly_fct_eval_mat = zeros(p+1, N);
    syms zeta zeta_2 u real

    % Pre-evaluate omega_2^2 at Gauss-Legendre nodes (increases performance)
    omega_2_fct = matlabFunction((omgn_m+omgn_d_2*zeta)^2, 'Vars', zeta);
    omega_2_eval_vec = omega_2_fct(gauss_legendre_coll{1}).';
    
    poly_x_mask = [];
    counter = 0;
    
    for n_1 = 0:p
        % Create Legendre basis
        poly_temp_1 = 1/(2^n_1*factorial(n_1)) * diff((zeta^2-1)^n_1,n_1);
        poly_vec(n_1+1) = poly_temp_1;

        % Pre-evaluate Legendre polynomials at Gauss-Legendre nodes (increases performance)
        poly_fct_temp = matlabFunction(poly_temp_1, 'Vars', zeta);
        poly_fct_eval_mat(n_1+1, :) = poly_fct_temp(gauss_legendre_coll{1}).';

        % Safe basis index combinations for x, x_dot and x_dot_dot for second interval
        for n_2 = 0:p-n_1
            counter = counter + 1;
            poly_x_mask = [poly_x_mask; [n_1 n_2]]; % Mask for x_dot is the same as for x
            
        end
    end

    % PCE coefficients
    x_dd_vec = sym('x_dd', [p+1 1]);
    x_d_vec = sym('x_d', [p+1 1]);
    x_vec = sym('x', [p+1 1]);
    
    k = (p+1)*(p+2)/2;
    x_dd_2_vec = sym('x_dd', [k 1]);
    x_d_2_vec = sym('x_d', [k 1]);
    x_2_vec = sym('x', [k 1]);

    LHS_expanded = poly_vec*x_dd_vec;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PCE: Instrusive Approach %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    RHS_1_expanded = -(omgn_m + omgn_d_1*zeta)^2*poly_vec*x_vec + (omgn_m + omgn_d_1*zeta)^2*u;

    % Galerkin Projecion
    counter = 0;
    PCE_A2_mat = zeros(k, k);

    for ind_1 = 1:p+1
        % Galerkin for first interval
        LHS_1_temp = LHS_expanded*poly_vec(ind_1); % F(Zeta)*Basis Function(ind)
        LHS_1_proj(ind_1) = int(LHS_1_temp,zeta,-1,1)*pdf_1;
        LHS_1_proj_coeff_temp = coeffs(LHS_1_proj(ind_1)); % extract coefficients from left-hand-side
        LHS_1_proj(ind_1) = LHS_1_proj(ind_1)/LHS_1_proj_coeff_temp;
        RHS_1_temp = RHS_1_expanded*poly_vec(ind_1);
        RHS_1_proj(ind_1) = int(RHS_1_temp,zeta,-1,1)*pdf_1/LHS_1_proj_coeff_temp;
        
        % Create expressions analytically to avoid int()
        for ind_2 = 1:p+2-ind_1
            counter = counter + 1;
            if ind_2 == 1
                PCE_b2_vec(counter) = (2*(ind_1-1)+1)/2 *dot(gauss_legendre_coll{2}, poly_fct_eval_mat(ind_1, :).*omega_2_eval_vec);
            else
                PCE_b2_vec(counter) = 0;
            end

            for i = 1:k
                multi_i = poly_x_mask(i, :);
                if multi_i(2) == ind_2-1
                    A2_temp_vec(i) = (2*(ind_1-1)+1)/2*dot(gauss_legendre_coll{2}, poly_fct_eval_mat(multi_i(1)+1, :).*poly_fct_eval_mat(ind_1, :).*(-omega_2_eval_vec));
                else
                    A2_temp_vec(i) = 0;
                end
            end
            PCE_A2_mat(counter, :) = A2_temp_vec;
        end
    end

    % Equations to Matrix
    syms u real
    [PCE_A1_mat, PCE_b1_vec] = equationsToMatrix(RHS_1_proj, x_vec);
    PCE_A1_mat = double(PCE_A1_mat); % otherwise matrix is symbolic in ode45
    [PCE_b1_vec, temp] = equationsToMatrix(-PCE_b1_vec, u); % extract u (note eqn2mat bring b to RHS)
    PCE_b1_vec = double(PCE_b1_vec);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Robust input shaper: GSA %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u = 1;
    results_coll = fct_cost_TDF_GSA(p, PCE_A1_mat, PCE_b1_vec, PCE_A2_mat, PCE_b2_vec, T2_rob, TDF_vec, u, k);
    time_vec(deg_i) = toc;
    
    t_vals_vec = results_coll{1};
    y1_vals_mat = results_coll{2};
    y2_vals_mat = results_coll{3};
    
    a_mean_vec = [y1_vals_mat(:, 1).' y2_vals_mat(:, 1).'];
    b_mean_vec = [y1_vals_mat(:, p+2).' y2_vals_mat(:, k+1).'];
    
    means_a_mat(deg_i, :) = a_mean_vec;
    means_b_mat(deg_i, :) = b_mean_vec;
    
    a_var_vec = fct_get_var(y1_vals_mat(:, 1:p+1), y2_vals_mat(:, 1:k), poly_x_mask);
    b_var_vec = fct_get_var(y1_vals_mat(:, p+2:end), y2_vals_mat(:, k+1:end), poly_x_mask);

    var_a_mat(deg_i, :) = a_var_vec;
    var_b_mat(deg_i, :) = b_var_vec;

end

save('means_PCE_mat.mat', 'means_a_mat', 'means_b_mat');
save('vars_PCE_mat.mat', 'var_a_mat', 'var_b_mat');


%%%%%%%%%%%%%%%%
%%% Analysis %%%
%%%%%%%%%%%%%%%%

figure;
bar(time_vec(:))
hold on;
yline(59.6135, 'r--', 'Time needed MC', 'LabelHorizontalAlignment','left')
xlabel('Degree of Polynom')
xticks(linspace(1, length(degrees_vec), length(degrees_vec)))
xticklabels(string(degrees_vec(:)))
ylabel('Seconds')
title('Time needed per max. Degree')
ax = gca;                       
ax.FontSize = 30;               
ax.XLabel.FontSize = 30;        
ax.YLabel.FontSize = 30;        
ax.Title.FontSize  = 30; 
hold off;

%%%%%%%%%%%%%%%%%
%%% Functions %%%
%%%%%%%%%%%%%%%%%

function J = fct_cost_TDF_rob(x,omgn)
    J = x(5); % time T2
end

function [c, ceq] = fct_noncon_TDF_rob(x,omgn)
    A0 = x(1); A1 = x(2); A2 = x(3); T1 = x(4); T2 = x(5);
    c = x(4) - x(5);
    ceq(1) = A0 + A1*cos(T1*omgn) + A2*cos(T2*omgn);
    ceq(2) = -A1*sin(T1*omgn) - A2*sin(T2*omgn);
    ceq(3) = -A1*T1*cos(T1*omgn) - A2*T2*cos(T2*omgn);
    ceq(4) = A1*T1*sin(T1*omgn) + A2*T2*sin(T2*omgn);
    ceq(5) = A0 + A1 + A2 - 1;
end

function set_up = fct_gauss_legendre_wn(num_points)
% Function to construct weights (w) and nodes (x) for gauss legendre
% integration
    k     = (1:num_points-1)';
    beta  = k ./ sqrt(4*k.^2 - 1);
    J     = diag(beta,1) + diag(beta,-1);
    [V,D] = eig(J);
    x     = diag(D);
    [x,ix]= sort(x);
    V     = V(:,ix);
    w     = 2 * (V(1,:).^2)';
    set_up = {x, w};
end

function J = fct_cost_TDF_GSA(p, PCE_A1_mat, PCE_b1_vec, PCE_A2_mat, PCE_b2_vec, T2_rob, x_TDF, u, k_1)
% This function solves the spring-mass differential eq. and returns the
% PCE coefficients
    options = odeset('RelTol', 1e-12, 'AbsTol',1e-12);
    t1_vec = linspace(0,10,101);
    t2_vec = linspace(10,20,101);
    a1_init_vec = zeros(2*(p+1),1);

    [TT_1,YY_1] = ode45(@(t,y) fct_inside_PCE(t,y,p,PCE_A1_mat,PCE_b1_vec,x_TDF,T2_rob,u), t1_vec, a1_init_vec, options);

    init_2_vec = [YY_1(end, 1:p+1) zeros(1, (p^2+p)/2) YY_1(end, p+2:end) zeros(1, (p^2+p)/2)];
    
    [TT_2,YY_2] = ode45(@(t,y) fct_inside_PCE(t,y,k_1-1,PCE_A2_mat,PCE_b2_vec.',x_TDF,T2_rob,u), t2_vec, init_2_vec, options);

    J = {[t1_vec t2_vec], YY_1, YY_2};
end

function dxdt = fct_inside_PCE(t,x,p,PCE_A_mat,PCE_b_vec,x_TDF, T2_rob,u)
% Function that applies TDF to PCE approximation
    u_mod = (x_TDF(1) + x_TDF(2)*heaviside(t-x_TDF(4)) + x_TDF(3)*heaviside(t-T2_rob)) * u;
    dxdt(1:p+1,1) = x(p+2:2*(p+1));
    dxdt(p+2:2*(p+1),1) = PCE_A_mat*x(1:p+1) + PCE_b_vec*u_mod;
end

function var_vec = fct_get_var(states_1_mat, states_2_mat, mask)
% Function to compute variance of state/velocity at t_2 using PCE
% coefficients
    for time_i = 1:size(states_1_mat, 1)
        var_vec_1(time_i) = sum(arrayfun(@(val, idx) val^2*1/(2*(idx-1) + 1),  states_1_mat(time_i, 2:end), 2:size(states_1_mat, 2)));
        var_vec_2(time_i) = sum(arrayfun(@(val, idx) val^2*1/(2*mask(idx, 1) + 1)*1/(2*mask(idx, 2) + 1),  states_2_mat(time_i, 2:end), 2:size(states_2_mat, 2)));
    end

    var_vec = [var_vec_1 var_vec_2];
end
