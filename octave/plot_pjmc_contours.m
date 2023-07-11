

%file_name1 = '../../data/pjmc_baseline.pjmc_data'
file_name2 = '../../data/pjmc_static_parallel_r10by10_np1e6_r75_aff1_th99.pjmc_data';
file_name2 = '../../data/scans/ranks/pjmc_static_r05by05_np1e5_r50_aff1_thresh99.pjmc_data'
%[delta_t, pj_iter, pj_ratio, pj_aff, N_1, szs, lower_bnds, upper_bnds, t_i_1, cpu_t_i, f_yx_1] = load_pjmc_data(file_name1);
[delta_t, pj_iter, pj_ratio, pj_aff, N, szs, lower_bnds, upper_bnds, t_i, cpu_t_i, f_yx] = load_pjmc_data(file_name2);

% Define the final distributuion function over the given grid
x = linspace(lower_bnds(1),upper_bnds(1),szs(1));
y = linspace(lower_bnds(2),upper_bnds(2),szs(2));
[X,Y] = meshgrid(y,x);
F_inf = (2/sqrt(pi))*Y.*exp(-(Y.^2+X.^2));

for i=1:3
    init_idx = (n_data_per_iter+1)*(i-1) + 1
    final_mic_idx = (n_data_per_iter+1)*(i-1) + n_data_per_iter 
    proj_idx = (n_data_per_iter+1)*(i-1) + n_data_per_iter + 1 
    subplot(3,3,1+(i-1)*3)
    contour(X,Y,f_yx(:,:,init_idx),10,'linewidth',2);
    title(['t = ', num2str(t_i(init_idx))]);
    xlabel('V_{||}');
    ylabel('V_{\perp}');
    subplot(3,3,2+(i-1)*3);
    contour(X,Y,f_yx(:,:,final_mic_idx),10,'linewidth',2);
    title(['t = ', num2str(t_i(final_mic_idx))]);
    xlabel('V_{||}');
    ylabel('V_{\perp}');
    subplot(3,3,3+(i-1)*3);
    contour(X,Y,f_yx(:,:,proj_idx),10,'linewidth',2);
    title(['t = ', num2str(t_i(proj_idx))]);
    xlabel('V_{||}');
    ylabel('V_{\perp}');
end 








