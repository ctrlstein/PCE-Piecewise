clear all; clc; close all;


%%%%%%%%%%%%%%
%%% System %%%
%%%%%%%%%%%%%%

omgn_m = pi;
omgn_d_1 = 0.25*pi;
omgn_d_2 = 0.5*pi;
pdf_1 = 1/2;
pdf_2 = 1/2;

% Initialise weights and nodes for gauss legendre integration
N = 100;
gauss_legendre_coll = fct_gauss_legendre_wn(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Robust input shaper %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


lb = [0 0 0 0 0];
ub = [1 1 1 2 5];
x0 = [0.3 0.3 0.3 1 2];
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
[x_opt, fval_opt] = fmincon(@(x)fct_cost_TDF_rob(x,omgn_m), x0, [], [], [], [], lb, ub, @(x)fct_noncon_TDF_rob(x,omgn_m), options);
A0_rob = x_opt(1); A1_rob = x_opt(2); A2_rob = x_opt(3); T1_rob = x_opt(4); T2_rob = x_opt(5);

x_opt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE: Basis Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 30; % polynomial degree
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Robust input shaper: GSA var %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = 1;
lb = [0 0 0 0];
ub = [1 1 1 T2_rob]; % set final time same as robust TDF
x0 = [0.25 0.5 0.25 1];
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');

[x_opt, var_opt] = fmincon(@(x)fct_cost_TDF_GSA_var(p, PCE_A1_mat, PCE_A2_mat, PCE_b1_vec, PCE_b2_vec, T2_rob, x, u, k, poly_x_mask, poly_fct_eval_mat, omega_2_eval_vec, gauss_legendre_coll{2}), ...
    x0, [], [], [], [], lb, ub, @(x)fct_noncon_TDF_GSA(x), options);
A0_GSA = x_opt(1); A1_GSA = x_opt(2); A2_GSA = x_opt(3); T1_GSA = x_opt(4); T2_GSA = T2_rob;
'GSA TDF optimized for variance of residual energy'
[A0_GSA, A1_GSA, A2_GSA, T1_GSA, T2_GSA]
var_opt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test res. energy GSA var %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TDF_rob_vec = [A0_rob, A1_rob, A2_rob, T1_rob, T2_rob];
TDF_GSA_vec = [A0_GSA, A1_GSA, A2_GSA, T1_GSA, T2_GSA];

t1_vec = linspace(0,10,1001);
t2_vec = linspace(10,20,1001);

omgn_d_1_vec = linspace(omgn_m-omgn_d_1,omgn_m+omgn_d_1,21);
omgn_d_2_vec = linspace(omgn_m-omgn_d_2,omgn_m+omgn_d_2,21);

for ind_1 = 1:length(omgn_d_1_vec)
    omgn_1 = omgn_d_1_vec(ind_1);
    x0_1 = [0, 0];
    [TT, YY_rob_1] = ode45(@(t,x) fct_inside_real_system(t,x,omgn_1,TDF_rob_vec,u), t1_vec, x0_1);
    [TT, YY_GSA_1] = ode45(@(t,x) fct_inside_real_system(t,x,omgn_1,TDF_GSA_vec,u), t1_vec, x0_1);

    for ind_2 = 1:length(omgn_d_2_vec)
        omgn_2 = omgn_d_2_vec(ind_2);
        x0_TDF_rob_2 = YY_rob_1(end,:);
        x0_TDF_GSA_2 = YY_GSA_1(end,:);

        [TT, YY_rob_2] = ode45(@(t,x) fct_inside_real_system(t,x,omgn_2,TDF_rob_vec,u), t2_vec, x0_TDF_rob_2);
        [TT, YY_GSA_2] = ode45(@(t,x) fct_inside_real_system(t,x,omgn_2,TDF_GSA_vec,u), t2_vec, x0_TDF_GSA_2);

        V_res_TDF_rob_mat(ind_1,ind_2) = 0.5*(YY_rob_2(end,2))^2 + omgn_2^2 * 0.5*(1 - YY_rob_2(end,1))^2;
        V_res_TDF_GSA_mat(ind_1,ind_2) = 0.5*(YY_GSA_2(end,2))^2 + omgn_2^2 * 0.5*(1 - YY_GSA_2(end,1))^2;
    end
end

figure;
surf(omgn_d_1_vec, omgn_d_2_vec, V_res_TDF_rob_mat, 'Facecolor', 'r', 'EdgeColor','none')
xlabel('omega_1')
ylabel('omega_2')
zlabel('Residual Energy')
hold on;
surf(omgn_d_1_vec, omgn_d_2_vec, V_res_TDF_GSA_mat, 'FaceColor','b', 'EdgeColor','none')
grid on; 
t = 'omega_1~U(0.75*pi,1.25*pi), omega_2~U(0.5*pi,1.5*pi), deg=30';
title(t)
ax = gca;                       
ax.FontSize = 30;               
ax.XLabel.FontSize = 30;        
ax.YLabel.FontSize = 30;        
ax.Title.FontSize  = 30; 
view(0, 90);
hold off;

figure;
lu_bounds = heatmap_f(V_res_TDF_rob_mat, V_res_TDF_GSA_mat, 'omega_1', 'omega_2', round(omgn_d_1_vec, 2), round(omgn_d_2_vec, 2), t, 0, []);

set(gcf,'renderer','painters');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Robust input shaper: GSA exp_val %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x_opt, expected_val_opt] = fmincon(@(x)fct_cost_TDF_GSA_exp_val(p, PCE_A1_mat, PCE_A2_mat, PCE_b1_vec, PCE_b2_vec, T2_rob, x, u, k, poly_x_mask, poly_fct_eval_mat, omega_2_eval_vec, gauss_legendre_coll{2}), ...
    x0, [], [], [], [], lb, ub, @(x)fct_noncon_TDF_GSA(x), options);
A0_GSA = x_opt(1); A1_GSA = x_opt(2); A2_GSA = x_opt(3); T1_GSA = x_opt(4); T2_GSA = T2_rob;
'GSA TDF optimized for expected residual energy'
[A0_GSA, A1_GSA, A2_GSA, T1_GSA, T2_GSA]
expected_val_opt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test res. energy GSA exp_val %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TDF_rob_vec = [A0_rob, A1_rob, A2_rob, T1_rob, T2_rob];
TDF_GSA_vec = [A0_GSA, A1_GSA, A2_GSA, T1_GSA, T2_GSA];

t1_vec = linspace(0,10,1001);
t2_vec = linspace(10,20,1001);

omgn_d_1_vec = linspace(omgn_m-omgn_d_1,omgn_m+omgn_d_1,21);
omgn_d_2_vec = linspace(omgn_m-omgn_d_2,omgn_m+omgn_d_2,21);

for ind_1 = 1:length(omgn_d_1_vec)
    omgn_1 = omgn_d_1_vec(ind_1);
    x0_1 = [0, 0];
    [TT, YY_rob_1] = ode45(@(t,x) fct_inside_real_system(t,x,omgn_1,TDF_rob_vec,u), t1_vec, x0_1);
    [TT, YY_GSA_1] = ode45(@(t,x) fct_inside_real_system(t,x,omgn_1,TDF_GSA_vec,u), t1_vec, x0_1);

    for ind_2 = 1:length(omgn_d_2_vec)
        omgn_2 = omgn_d_2_vec(ind_2);
        x0_TDF_rob_2 = YY_rob_1(end,:);
        x0_TDF_GSA_2 = YY_GSA_1(end,:);

        [TT, YY_rob_2] = ode45(@(t,x) fct_inside_real_system(t,x,omgn_2,TDF_rob_vec,u), t2_vec, x0_TDF_rob_2);
        [TT, YY_GSA_2] = ode45(@(t,x) fct_inside_real_system(t,x,omgn_2,TDF_GSA_vec,u), t2_vec, x0_TDF_GSA_2);

        V_res_TDF_rob_mat(ind_1,ind_2) = 0.5*(YY_rob_2(end,2))^2 + omgn_2^2 * 0.5*(1 - YY_rob_2(end,1))^2;
        V_res_TDF_GSA_mat(ind_1,ind_2) = 0.5*(YY_GSA_2(end,2))^2 + omgn_2^2 * 0.5*(1 - YY_GSA_2(end,1))^2; 
    end
end

figure;
surf(omgn_d_1_vec, omgn_d_2_vec, V_res_TDF_rob_mat, 'Facecolor', 'r', 'EdgeColor','none')
xlabel('omega_1')
ylabel('omega_2')
zlabel('Residual Energy')
hold on;
surf(omgn_d_1_vec, omgn_d_2_vec, V_res_TDF_GSA_mat, 'FaceColor','b', 'EdgeColor','none')
grid on; 
t = 'omega_1~U(0.75*pi,1.25*pi), omega_2~U(0.5*pi,1.5*pi), deg=30';
title(t)
ax = gca;                       
ax.FontSize = 30;               
ax.XLabel.FontSize = 30;        
ax.YLabel.FontSize = 30;        
ax.Title.FontSize  = 30; 
view(0, 90);
hold off;

figure;
heatmap_f(V_res_TDF_rob_mat, V_res_TDF_GSA_mat, 'omega_1', 'omega_2', round(omgn_d_1_vec, 2), round(omgn_d_2_vec, 2), t, 1, lu_bounds)

set(gcf,'renderer','painters');


%%%%%%%%%%%%%%%%%
%%% Functions %%%
%%%%%%%%%%%%%%%%%

function bounds = heatmap_f(tdf, gsa, xlab, ylab, xt, yt, tit, load_cmap, bounds)
    heatmap_mat = gsa - tdf;
    h = heatmap(xt, yt, heatmap_mat);
    
    if isempty(bounds)
        upper = max(heatmap_mat(:));
        lower = min(heatmap_mat(:));
    else
        upper = bounds(2);
        lower = bounds(1);
    end
    rel = -lower/(upper - lower);
    n = 512;
    n1 = round(n * rel);
    n2 = n - n1;
    
    green = [0 1 0];
    yellow = [1 1 0];
    red = [1 0 0];

    cmap1 = [linspace(green(1), yellow(1), n1)', ...
         linspace(green(2), yellow(2), n1)', ...
         linspace(green(3), yellow(3), n1)'];

    cmap2 = [linspace(yellow(1), red(1), n2)', ...
         linspace(yellow(2), red(2), n2)', ...
         linspace(yellow(3), red(3), n2)'];
    
    % Use same colormap for both GSA TDFs to improve comparability
    if load_cmap == 0
        cmap = [cmap1; cmap2];
        save('cmap.mat', 'cmap');
    else
        load('cmap.mat', 'cmap');
    end
    
    colormap(cmap);
    h.ColorLimits = [lower upper];
    h.FontSize = 30;  
    h.XLabel = xlab;
    h.YLabel = ylab;
    h.Title = tit;
    h.GridVisible = 'off';
    h.CellLabelColor = 'none';

    bounds = [lower, upper];
   
end

function exp_res_energy = fct_cost_TDF_GSA_exp_val(p, PCE_A1_mat, PCE_A2_mat, PCE_b1_vec, PCE_b2_vec, T2_rob, x_TDF, u, k_1, mask_pol, poly_eval_mat, omgn_2_eval_vec, gl_weights)
% This function solves the spring-mass differential eq. and returns the
% expected residual energy at t_2
    t1_vec = linspace(0,10,101);
    t2_vec = linspace(10,20,101);
    a1_init_vec = zeros(2*(p+1),1);
    [TT_1,YY_1] = ode45(@(t,y) fct_inside_PCE(t,y,p,PCE_A1_mat,PCE_b1_vec,x_TDF,T2_rob,u), t1_vec, a1_init_vec);
    
    init_2_vec = [YY_1(end, 1:p+1) zeros(1, (p^2+p)/2) YY_1(end, p+2:end) zeros(1, (p^2+p)/2)];
    
    [TT_2,YY_2] = ode45(@(t,y) fct_inside_PCE(t,y,k_1-1,PCE_A2_mat,PCE_b2_vec.',x_TDF,T2_rob,u), t2_vec, init_2_vec);
    
    
    % calculate the expected value of the res. energy
    exp_res_energy = get_exp_res_energy(YY_2(end, :), mask_pol, poly_eval_mat, omgn_2_eval_vec, gl_weights);
    
    % var_res_en = get_var_res_energy(YY_2(end, :), mask_pol, poly_eval_mat, omgn_2_eval_vec, gl_weights,exp_res_energy) % can be used to print variance
end

function var_res_en  = fct_cost_TDF_GSA_var(p, PCE_A1_mat, PCE_A2_mat, PCE_b1_vec, PCE_b2_vec, T2_rob, x_TDF, u, k_1, mask_pol, poly_eval_mat, omgn_2_eval_vec, gl_weights)
% This function solves the spring-mass differential eq. and returns the
% variance of the residual energy at t_2
    t1_vec = linspace(0,10,101);
    t2_vec = linspace(10,20,101);
    a1_init_vec = zeros(2*(p+1),1);
    [TT_1,YY_1] = ode45(@(t,y) fct_inside_PCE(t,y,p,PCE_A1_mat,PCE_b1_vec,x_TDF,T2_rob,u), t1_vec, a1_init_vec);
    
    init_2_vec = [YY_1(end, 1:p+1) zeros(1, (p^2+p)/2) YY_1(end, p+2:end) zeros(1, (p^2+p)/2)];
    
    [TT_2,YY_2] = ode45(@(t,y) fct_inside_PCE(t,y,k_1-1,PCE_A2_mat,PCE_b2_vec.',x_TDF,T2_rob,u), t2_vec, init_2_vec);
    
    
    % calculate the variance of the res. energy
    exp_en = get_exp_res_energy(YY_2(end, :), mask_pol, poly_eval_mat, omgn_2_eval_vec, gl_weights);
    var_res_en = get_var_res_energy(YY_2(end, :), mask_pol, poly_eval_mat, omgn_2_eval_vec, gl_weights,exp_en);
    % exp_en % can be used to print expected values
end


function [c, ceq] = fct_noncon_TDF_GSA(x_TDF)
c = [];
ceq = x_TDF(1) + x_TDF(2) + x_TDF(3) - 1;
end


function dxdt = fct_inside_PCE(t,x,p,PCE_A_mat,PCE_b_vec,x_TDF, T2_rob,u)
% Function that applies TDF to PCE approximation
u_mod = (x_TDF(1) + x_TDF(2)*heaviside(t-x_TDF(4)) + x_TDF(3)*heaviside(t-T2_rob)) * u;
dxdt(1:p+1,1) = x(p+2:2*(p+1));
dxdt(p+2:2*(p+1),1) = PCE_A_mat*x(1:p+1) + PCE_b_vec*u_mod;
end

function dxdt = fct_inside_real_system(t,x,omgn,TDF_vec,u)
% Function that applies TDF to spring-mass system
u_mod = (TDF_vec(1) + TDF_vec(2)*heaviside(t-TDF_vec(4)) + TDF_vec(3)*heaviside(t-TDF_vec(5))) * u;
dxdt(1,1) = x(2);
dxdt(2,1) = -omgn^2*x(1) + omgn^2*u_mod;
end


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

function e_res = get_exp_res_energy(coeff_vec, mask_mat, mask_eval_mat, omgn_2_eval, weights_gl)
% Function that computes the expected residual energy at t_2 using gauss
% legendre integration
    K = length(mask_mat);
    N = length(weights_gl);
    e_res = 0;

    for i = 1:N
        fi_vec = zeros(1,N);
        for j = 1:N
            pol_x = 0;
            pol_x_dot = 0;
            for k_2 = 1:K
                multi_k = mask_mat(k_2, :);
                pol_x = pol_x + coeff_vec(k_2)*mask_eval_mat(multi_k(1)+1, j)*mask_eval_mat(multi_k(2)+1, i);
                pol_x_dot = pol_x_dot + coeff_vec(K+k_2)*mask_eval_mat(multi_k(1)+1, j)*mask_eval_mat(multi_k(2)+1, i);
            end
            
            fi_vec(j) = pol_x_dot^2 + omgn_2_eval(j)*(1-pol_x)^2;
        end
        e_res = e_res + weights_gl(i) * dot(weights_gl, fi_vec);
    end
    
    e_res = 1/8 * e_res;
end

function var_e_res = get_var_res_energy(coeff_vec, mask_mat, mask_eval_mat, omgn_2_eval, weights_gl, e_res)
% Function that computes the varaince of the residual energy at t_2 using gauss
% legendre integration
    K = length(mask_mat);
    N = length(weights_gl);
    e_res_2 = 0;

    for i = 1:N
        for j = 1:N
            pol_x = 0;
            pol_x_dot = 0;
            for k_2 = 1:K
                multi_k = mask_mat(k_2, :);
                pol_x = pol_x + coeff_vec(k_2)*mask_eval_mat(multi_k(1)+1, j)*mask_eval_mat(multi_k(2)+1, i);
                pol_x_dot = pol_x_dot + coeff_vec(K+k_2)*mask_eval_mat(multi_k(1)+1, j)*mask_eval_mat(multi_k(2)+1, i);
            end
            
            fi_vec(j) = (pol_x_dot^2 + omgn_2_eval(j)*(1-pol_x)^2)^2;
        end
        e_res_2 = e_res_2 + weights_gl(i) * dot(weights_gl, fi_vec);
    end

    e_res_2 = 1/16 * e_res_2;

    var_e_res = e_res_2 - e_res^2;
end
