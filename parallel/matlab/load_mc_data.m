% Carlos del-Castillo-Negrete
% April, 2019
% FAMOSO - a FAst MOnte carlo SOlver - load_mc_data.m
% Function to load a .mc_data file that contains simulation data as ouptuted by the
% mc_data_output() routin in the montecarlo module. 

function [mc_data] = load_mc_data(file_name)

%ORDER OF OUTPUT IN MC_DATA_OUTPUT:
%write(1,*) mc_prob%delta_t, max_gd%szs, max_gd%lower_bnds, max_gd%upper_bnds, &
%   &N, mc_prob%mc_data%cpu_t_i, mc_prob%mc_data%t_i, A 

raw_data =   textread(file_name);
mc_data.nt = raw_data(1,1);
mc_data.np = raw_data(2,1);
mc_data.delta_t = raw_data(3,1);
mc_data.szs = raw_data(4,:);
mc_data.lower_bnds = raw_data(5,1:2);
mc_data.upper_bnds = raw_data(6,1:2);
mc_data.N = raw_data(7,1);
mc_data.t_i = raw_data(8,1:mc_data.N);
mc_data.cpu_init = raw_data(9,1:2);
cpu_t_i_flat = raw_data(10,1:4*(mc_data.N-1));
f_flat = raw_data(11:size(raw_data,1),1);

mc_data.cpu_t_i = reshape(cpu_t_i_flat,mc_data.N-1,4);
mc_data.cpu_t_i(:,:) = mc_data.cpu_t_i(:,:) - mc_data.cpu_init(1);
mc_data.cpu_init(:) = mc_data.cpu_init(:)-mc_data.cpu_init(1);
mc_data.f_yx = reshape(f_flat,mc_data.szs(1),mc_data.szs(2),mc_data.N);



