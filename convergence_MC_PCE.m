load('means_MC_vec.mat')
load('vars_MC_vec.mat')

load('means_PCE_mat.mat')
load('vars_PCE_mat.mat')

colors = [ 
    0.0   0.4   0.8
    0.2   0.6   1.0
    0.4   0.75  1.0
];

labels = { ...
    'PCE deg = 10', ...
    'PCE deg = 20', ...
    'PCE deg = 30', ...
    'Monte Carlo'};


% Expected State
h_1 = gobjects(1,4);
figure;
hold on;
for i= 1:3
    h_1(i) = plot(means_a_mat(i, ([1:100, 102:end])), 'Color', colors(3-i+1, :));
end
h_1(4) = plot(y_1_mean_vec([1:100, 102:end]), 'r');

xlabel('t');
ylabel('E[X]');
title('E[X] Monte Carlo vs PCE');
legend(h_1, labels, 'Location', 'best');
ax = gca;                       
ax.FontSize = 30;               
ax.XLabel.FontSize = 30;        
ax.YLabel.FontSize = 30;        
ax.Title.FontSize  = 30; 
hold off;

% Expected velocity

h_2 = gobjects(1,4);
figure;
hold on;
for i= 1:3
    h_2(i) = plot(means_b_mat(i, ([1:100, 102:end])), 'Color', colors(3-i+1, :));
end
h_2(4) = plot(y_2_mean_vec([1:100, 102:end]), 'r');

xlabel('t');
ylabel('E[X(dot)]');
title('E[X(dot)] Monte Carlo vs PCE');
legend(h_2, labels, 'Location', 'best');
ax = gca;                       
ax.FontSize = 30;               
ax.XLabel.FontSize = 30;        
ax.YLabel.FontSize = 30;        
ax.Title.FontSize  = 30; 
hold off;


% Variance State

h_3 = gobjects(1,4);
figure;
hold on;
for i= 1:3
    h_3(i) = plot(var_a_mat(i, ([1:100, 102:end])), 'Color', colors(3-i+1, :));
end
h_3(4) = plot(y_1_var_vec([1:100, 102:end]), 'r');

xlabel('t');
ylabel('Var(X)');
title('Var(X) Monte Carlo vs PCE');
legend(h_3, labels, 'Location', 'best');
ax = gca;                       
ax.FontSize = 30;               
ax.XLabel.FontSize = 30;        
ax.YLabel.FontSize = 30;        
ax.Title.FontSize  = 30; 
hold off;

% Variance velocity

h_4 = gobjects(1,4);
figure;
hold on;
for i= 1:3
        h_4(i) = plot(var_b_mat(i, ([1:100, 102:end])), 'Color', colors(3-i+1, :));
end
h_4(4) = plot(y_2_var_vec([1:100, 102:end]), 'r');

xlabel('t');
ylabel('Var(X(dot))');
title('Var(X(dot)) Monte Carlo vs PCE');
legend(h_4, labels, 'Location', 'best');
ax = gca;                       
ax.FontSize = 30;               
ax.XLabel.FontSize = 30;        
ax.YLabel.FontSize = 30;        
ax.Title.FontSize  = 30; 
hold off;
