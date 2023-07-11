
clear all

file_name1 = '../../data/pjmc_static_baseline.pjmc_data'
file_name2 = '../../data/scans/ranks/pjmc_static_r05by05_np1e5_r50_aff1_thresh99.pjmc_data'
[delta_t, pj_iter, pj_ratio, pj_aff, N_1, szs, lower_bnds, upper_bnds, t_i_1, cpu_t_i, f_yx_1] = load_pjmc_data(file_name1);
[delta_t, pj_iter, pj_ratio, pj_aff, N, szs, lower_bnds, upper_bnds, t_i, cpu_t_i, f_yx] = load_pjmc_data(file_name2);

% Define the final distributuion function over the given grid
v_perp = linspace(lower_bnds(1),upper_bnds(1),szs(1));
v_par = linspace(lower_bnds(2),upper_bnds(2),szs(2));
[V_PAR_GRID,V_PERP_GRID] = meshgrid(v_par,v_perp);
F_inf = (2/sqrt(pi))*V_PERP_GRID.*exp(-(V_PERP_GRID.^2+V_PAR_GRID.^2));

n_data = floor(1/pj_aff);
e_thresh = 0.99;

% Calculate error for each distributuion function
for i=1:N
    t_i(i);
    err_1(i) = sum(sum((F_inf-f_yx(:,:,i)).^2));
end 
% Calculate error for each distributuion function
for i=1:N_1
    err_2(i) = sum(sum((F_inf-f_yx_1(:,:,i)).^2));
end 


figure(1)
set(gcf,'defaulttextfontsize',18,'defaultaxesfontsize',22)
set(gcf,'defaultlinelinewidth',2.0)
set(gca,'box','on')
set(gca,'linewidth',2.5)
set(gca,'fontsize',25)

plot(t_i(1:30),err_1(1:30),'r-o', 'Markersize', 10, 'linewidth', 3.5); hold on;
%plot(t_i_1,err_2(:),'b');
%title(['Convergence of PMC - PJ_{ratio} = ', num2str(pj_ratio), ...
%        ', (r_1,r_2) = (5,5),  E_{th} = ', num2str(e_thresh)]);
xlabel('Time');
ylabel('||F_i - F_{\infty}||^2');
%legend('Projective Method','Standard Integration');
input('');
