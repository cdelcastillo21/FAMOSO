% Carlos del-Castillo-Negrete
% April, 2019
% create_figs - Script used to read mc data and produce, weak/strong
% scaling graphs and compute time bar graphs.

close all


5 plot_nscan('../data/10M-N-Scan',1)

function plot_nscan(dr, fi)
    dinfo = dir(dr);
    files = {dinfo.name}
    start_idx = 2;
    nf = size(files,2);
    cpu_ts = zeros(nf-start_idx,4);
    np = zeros(nf-start_idx,1);
    nt = zeros(nf-start_idx,1);
    tt = zeros(nf-start_idx,1);
    N = zeros(nf-start_idx,1);

    for i=1:(size(files,2)-start_idx)
        d = load_mc_data(strcat(dr,"/",files{i+start_idx}))
        N_iters = size(d.cpu_t_i,1);
        
        N(i) = d.N;
        nt(i) = d.nt;
        np(i) = d.np;
        cpu_ts(i,1) = d.cpu_init(2);
        for j=1:N_iters
            cpu_ts(i,2) = cpu_ts(i,2) + d.cpu_t_i(j,2)-d.cpu_t_i(j,1); % step 
            cpu_ts(i,3) = cpu_ts(i,3) + d.cpu_t_i(j,3)-d.cpu_t_i(j,2); % fit_grid
            cpu_ts(i,4) = cpu_ts(i,4) + d.cpu_t_i(j,4)-d.cpu_t_i(j,3); % xy_to_f
        end
        tt(i) = d.cpu_t_i(N_iters,4);
    end
      
    figure(fi)
    percent_time = cpu_ts(:,1:4)./tt;
    percent_time(:,1) = percent_time(:,2);
    percent_time(:,2) = 1-percent_time(:,1);
    
    bar(N,tt,'stacked')
    xlabel('Number of Snapshots'); ylabel('Total Runtime');
    title(strcat('NP = ', num2str(np(1)), ', NT = ', num2str(nt(1))))
    
    figure(fi+1)
    percent_time = cpu_ts(:,1:4)./tt;
    percent_time(:,1) = percent_time(:,2);
    percent_time(:,2) = 1-percent_time(:,1);
    title(strcat('NP = ', num2str(np), ', NT = ', num2str(nt)))
    bar(N,percent_time(:,1:2),'stacked')
    legend(['compute';'comm   '],'location','northwest');
    xlabel('Number of Snapshots'); ylabel('Percent Time Spent');
    title(strcat('NP = ', num2str(np(1)), ', NT = ', num2str(nt(1))))
    
    % Find min nt and index
    [m,m_i] = min(nt); 
    
end

function plot_rt(dr, fi)
    dinfo = dir(dr);
    files = {dinfo.name}
    start_idx = 2;
    nf = size(files,2);
    cpu_ts = zeros(nf-start_idx,4);
    np = zeros(nf-start_idx,1);
    nt = zeros(nf-start_idx,1);
    tt = zeros(nf-start_idx,1);

    for i=1:(size(files,2)-start_idx)
        d = load_mc_data(strcat(dr,"/",files{i+start_idx}))
        N_iters = size(d.cpu_t_i,1);
        
        nt(i) = d.nt;
        np(i) = d.np;
        cpu_ts(i,1) = d.cpu_init(2);
        for j=1:N_iters
            cpu_ts(i,2) = cpu_ts(i,2) + d.cpu_t_i(j,2)-d.cpu_t_i(j,1); % step 
            cpu_ts(i,3) = cpu_ts(i,3) + d.cpu_t_i(j,3)-d.cpu_t_i(j,2); % fit_grid
            cpu_ts(i,4) = cpu_ts(i,4) + d.cpu_t_i(j,4)-d.cpu_t_i(j,3); % xy_to_f
        end
        tt(i) = d.cpu_t_i(N_iters,4);
    end
      
    figure(fi)
    percent_time = cpu_ts(:,1:4)./tt;
    percent_time(:,1) = percent_time(:,2);
    percent_time(:,2) = 1-percent_time(:,1);
    title(strcat('NP = ', num2str(np), ', NT = ', num2str(nt)))
    bar(nt,tt,'stacked')
    xlabel('Number of MPI Tasks'); ylabel('Total Runtime');
    
    figure(fi+1)
    percent_time = cpu_ts(:,1:4)./tt;
    percent_time(:,1) = percent_time(:,2);
    percent_time(:,2) = 1-percent_time(:,1);
    title(strcat('NP = ', num2str(np), ', NT = ', num2str(nt)))
    bar(nt,percent_time(:,1:2),'stacked')
    legend(['compute';'comm   '],'location','northwest');
    xlabel('Number of MPI Tasks'); ylabel('Percent Time Spent');
    
    % Find min nt and index
    [m,m_i] = min(nt); 
    
    % Weak Scaling
    if 1
        % Scaled speed up
        tt = tt./np;
        t_1 = tt(m_i)/m;
        tt = t_1./tt;
        figure(fi+2)
        scatter(nt,tt,'or')
        hold on;
        
        p = polyfit(nt,tt,1)
        xs = 1:1:max(nt);
        figure(fi+2)
        plot(xs,(1-p(1))+p(1)*xs,'--k');
        xlabel('Number of Tasks','FontSize',20);
        ylabel('Scaled Speed-Up','FontSize',20);
        title(strcat('S = ', num2str(1-p(1))),'FontSize',26);
    end
    
    % Strong Scaling
    if 1
        t_1 = tt(m_i)/m;
        y = nt-(t_1./tt);
        x = (t_1./tt).*(nt-1);
        p = polyfit(x,y,1)
        xs = 1:1:max(nt);
        figure(fi+2)
        tt = t_1./tt;
        scatter(nt,tt,'or')
        hold on;
        plot(xs,1./(p(1)+(1-p(1))./xs),'--k');
        xlabel('Number of Tasks','FontSize',20);
        ylabel('Speed-Up','FontSize',20);
        title(strcat('S = ', num2str(p(1))),'FontSize',26);
    end
    
    
    plot_mc(strcat(dr,"/",files{1+start_idx}),fi+3)
end

function plot_mc(file, fi)
    [mc_data] = load_mc_data(file);

    delta_t = mc_data.delta_t;
    szs = mc_data.szs;
    lower_bnds = mc_data.lower_bnds;
    upper_bnds = mc_data.upper_bnds;
    N = mc_data.N;
    t_i = mc_data.t_i;
    cpu_init = mc_data.cpu_init;
    cpu_t_i = mc_data.cpu_t_i;
    f_yx = mc_data.f_yx;

    % Define the final distributuion function over the given grid
    v_perp = linspace(lower_bnds(1),upper_bnds(1),szs(1));
    v_par = linspace(lower_bnds(2),upper_bnds(2),szs(2));
    [V_PAR_GRID,V_PERP_GRID] = meshgrid(v_par,v_perp);
    F_inf = (2/sqrt(pi))*V_PERP_GRID.*exp(-(V_PERP_GRID.^2+V_PAR_GRID.^2));


    err = zeros(1,N);

    % Calculate error for each distributuion function
    for i=1:N
        t_i(i);
        err(i) = sum(sum((F_inf-f_yx(:,:,i)).^2))/norm(F_inf);
    end
    
    
    figure(fi)
    set(gcf,'defaulttextfontsize',18,'defaultaxesfontsize',22)
    set(gcf,'defaultlinelinewidth',2.0)
    set(gca,'box','on')
    set(gca,'linewidth',2.5)
    set(gca,'fontsize',25)

    plot((1:1:N),err,'r-o', 'Markersize', 10, 'linewidth', 3.5); hold on;
    xlabel('N');
    ylabel('||F_i - F_{\infty}||^2');

    idxs = N:-ceil(N/4):1;
    
    figure(fi+1)
    for i=1:size(idxs,2)
        subplot(2,2,i)
        contourf(V_PAR_GRID,V_PERP_GRID,f_yx(:,:,i))
    end
end
