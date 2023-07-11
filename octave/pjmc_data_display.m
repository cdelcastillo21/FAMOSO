% Carlos del-Castillo-Negrete
% April, 2015
% FAMOSO - a FAst MOnte carlo SOlver - pjmc_data_display.m
% Octave/matlab function to plot one by one the recorded distrubution function estimations of 
% a .pjmc_data simulation data file corresponding to a PMC simulation.
function pjmc_data_display(file_name)

[delta_t, pj_iter, pj_ratio, pj_aff, N, szs, lower_bnds, upper_bnds, t_i, cpu_t_i, f_yx] = load_pjmc_data(file_name)

% Define the final distributuion function over the given grid
x = linspace(lower_bnds(1),upper_bnds(1),szs(1));
y = linspace(lower_bnds(2),upper_bnds(2),szs(2));
[X,Y] = meshgrid(y,x);

figure(1)
set(gcf,'defaulttextfontsize',16,'defaultaxesfontsize',20)
set(gcf,'defaultlinelinewidth',2.0)
set(gca,'box','on')
for i=1:N
    contour(X,Y, f_yx(:,:,i),10,'linewidth',2);
    title( sprintf ('Time = %f', t_i(i)) )
    xlabel('V_{||}');
    ylabel('V_{\perp}');
    input("");
end
