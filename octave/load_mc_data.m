% Carlos del-Castillo-Negrete
% April, 2015
% FAMOSO - a FAst MOnte carlo SOlver - load_mc_data.m
% Function to load a .mc_data file that contains simulation data as ouptuted by the
% mc_data_output() routin in the montecarlo module. 
function [delta_t, N, nx, ny, lower_bnds, upper_bnds, cpu_t_i, t_i, f_yx] = load_mc_data(file_name)

%ORDER OF OUTPUT IN MC_DATA_OUTPUT:
%write(1,*) mc_prob%delta_t, max_gd%szs, max_gd%lower_bnds, max_gd%upper_bnds, &
%   &N, mc_prob%mc_data%cpu_t_i, mc_prob%mc_data%t_i, A 

raw_data = load(file_name);
delta_t = raw_data(1);
nx = raw_data(2);
ny = raw_data(3);
lower_bnds = raw_data(4:5);
upper_bnds = raw_data(6:7);
N = raw_data(8);
cpu_t_i = raw_data(9:(9+N-1));
t_i = raw_data((9+N):(9+2*N-1));
f_yx = reshape(raw_data((9+2*N):length(raw_data)),nx,ny,N);

