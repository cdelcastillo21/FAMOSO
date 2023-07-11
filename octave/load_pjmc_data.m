% Carlos del-Castillo-Negrete
% April, 2015
% FAMOSO - a FAst MOnte carlo SOlver - load_pjmc_data.m
% Function to load a .pjmc_data file that contains simulation data as ouptuted by the
% pjmc_data_output() routine in the projective_montecarlo module. 
function [delta_t, pj_iter, pj_ratio, pj_aff, N, szs, lower_bnds, upper_bnds, t_i, cpu_t_i, f_yx] = load_pjmc_data(file_name)

%ORDER OF OUTPUT IN PJMC_DATA_OUTPUT:
%write(1,*) pjmc_prob%delta_t, pjmc_prob%num_iter, pjmc_prob%pj_ratio, pjmc_prob%pj_aff, &
%    &N, max_gd%szs, max_gd%lower_bnds, max_gd%upper_bnds, pjmc_prob%mc_data%t_i, pjmc_prob%mc_data%cpu_t_i, A 

raw_data = load(file_name);
delta_t = raw_data(1);
pj_iter = raw_data(2);
pj_ratio = raw_data(3);
pj_aff = raw_data(4);
N = raw_data(5);
szs = raw_data(6:7);
lower_bnds = raw_data(8:9);
upper_bnds = raw_data(10:11);
t_i = raw_data(12:(12+N-1));
cpu_t_i = raw_data((12+N):(12+2*N-1));
f_yx = reshape(raw_data((12+2*N):length(raw_data)),szs(1),szs(2),N);

