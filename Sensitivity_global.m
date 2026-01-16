%% =========================================================
%  main_pipeline_analysis.m
%  Academic pipeline model
%  Option A + OAT sensitivity + PRCC
% =========================================================
clear; clc; close all;

%% ===========================
%  Force LaTeX everywhere
% ===========================
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ===========================
%  Load observed degree data
% ===========================
T = readtable('NCESdata.csv');

T.BS  = str2double(erase(T.BS, ','));
T.MS  = str2double(erase(T.MS, ','));
T.PHD = str2double(erase(T.PHD, ','));

year_data = T.Year';
BS_data   = T.BS';
MS_data   = T.MS';
PhD_data  = T.PHD';

%% ===========================
%  Baseline parameters
% ===========================
p.gU   = 0.14;
p.gG   = 0.17;
p.aG   = 0.08;
p.aP   = 0.25;
p.aF   = 0.03;

p.rMS  = 0.798;
p.rPhD = 0.202;

p.pGF     = 0.08;
p.pGP     = 0.45;
p.pPF_max = 0.18;

p.KF      = 4000;
p.alphaF  = 1.0;

p.P0 = 50;
p.F0 = 200;

%% ===========================
%  Baseline simulation
% ===========================
sim_base = simulate_pipeline_optionA(year_data, BS_data, MS_data, PhD_data, p);

%% ===========================
%  One at a time sensitivity
% ===========================
sens_params   = {'pGF','pGP','pPF_max','aP','aF','KF','alphaF'};
scale_factors = [0.5 0.75 1.0 1.25 1.5];

param_latex = struct( ...
    'pGF',      '$p^{GF}$', ...
    'pGP',      '$p^{GP}$', ...
    'pPF_max',  '$p^{PF}_{\max}$', ...
    'aP',       '$a^P$', ...
    'aF',       '$a^F$', ...
    'KF',       '$K_F$', ...
    'alphaF',   '$\alpha^F$' );

results_OAT = sensitivity_OAT(year_data, BS_data, MS_data, PhD_data, ...
                              p, sens_params, scale_factors);

for s = 1:numel(sens_params)

    param = sens_params{s};
    sf    = scale_factors;

    F_end = results_OAT.(param).F_end;
    F_max = results_OAT.(param).F_max;
    P_end = results_OAT.(param).P_end;
    P_max = results_OAT.(param).P_max;

    ptex = param_latex.(param);

    figure;

    subplot(2,2,1)
    plot(sf, F_end, 'o-','LineWidth',1.5)
    xlabel([ptex ])
    ylabel('$F_{\mathrm{final}}$')
    %title(['$F_{\mathrm{final}}$ vs ' ptex])

    subplot(2,2,2)
    plot(sf, F_max, 'o-','LineWidth',1.5)
    xlabel([ptex ])
    ylabel('$F_{\max}$')
    %title(['$F_{\max}$ vs ' ptex])

    subplot(2,2,3)
    plot(sf, P_end, 'o-','LineWidth',1.5)
    xlabel([ptex ])
    ylabel('$P_{\mathrm{final}}$')
    %title(['$P_{\mathrm{final}}$ vs ' ptex])

    subplot(2,2,4)
    plot(sf, P_max, 'o-','LineWidth',1.5)
    xlabel([ptex ])
    ylabel('$P_{\max}$')

end

%% ===========================
%  PRCC analysis
% ===========================
param_names = sens_params;
N_samples   = 500;

prcc_results = prcc_analysis(year_data, BS_data, MS_data, PhD_data, ...
                             p, param_names, N_samples);

figure
bar(prcc_results(:,1))

xt = {'$p^{GF}$','$p^{GP}$','$p^{PF}_{\max}$','$a^P$','$a^F$','$K_F$','$\alpha^F$'};
set(gca,'XTick',1:numel(xt),'XTickLabel',xt,'XTickLabelRotation',45)
set(gca,'TickLabelInterpreter','latex')

ylabel('$\mathrm{PRCC\ with\ }F_{\mathrm{final}}$')
%title('$\mathrm{PRCC:\ sensitivity\ of\ final\ faculty\ size}$')

%% =========================================================
%  FUNCTIONS
% =========================================================
function sim = simulate_pipeline_optionA(year_data, BS_data, MS_data, PhD_data, p)

years = year_data(:)';
TT    = numel(years);

U = BS_data(:)' / p.gU;
G = (MS_data(:)' + PhD_data(:)') / p.gG;

P = zeros(1,TT);
F = zeros(1,TT);

P(1) = p.P0;
F(1) = p.F0;

for k = 1:(TT-1)
    pPF = p.pPF_max / (1 + p.alphaF * F(k)/p.KF);

    P(k+1) = (1 - pPF - p.aP) * P(k) + p.pGP * p.gG * G(k);
    F(k+1) = (1 - p.aF) * F(k) + p.pGF * p.gG * G(k) + pPF * P(k);
end

sim.years   = years;
sim.P       = P;
sim.F       = F;
end

function results = sensitivity_OAT(year_data, BS_data, MS_data, PhD_data, p_base, sens_params, scale_factors)

for s = 1:numel(sens_params)
    param = sens_params{s};

    for i = 1:numel(scale_factors)
        p = p_base;
        p.(param) = p_base.(param) * scale_factors(i);

        sim = simulate_pipeline_optionA(year_data, BS_data, MS_data, PhD_data, p);

        results.(param).F_end(i) = sim.F(end);
        results.(param).F_max(i) = max(sim.F);
        results.(param).P_end(i) = sim.P(end);
        results.(param).P_max(i) = max(sim.P);
    end
end
end

function prcc_matrix = prcc_analysis(year_data, BS_data, MS_data, PhD_data, p_base, param_names, N)

nPar = numel(param_names);
Y    = zeros(N,4);
X    = rand(N,nPar);

for j = 1:nPar
    base = p_base.(param_names{j});
    X(:,j) = 0.5*base + (1.5*base - 0.5*base)*X(:,j);
end

for k = 1:N
    p = p_base;
    for j = 1:nPar
        p.(param_names{j}) = X(k,j);
    end

    sim = simulate_pipeline_optionA(year_data, BS_data, MS_data, PhD_data, p);
    Y(k,:) = [sim.F(end), max(sim.F), sim.P(end), max(sim.P)];
end

for m = 1:4
    prcc_matrix(:,m) = compute_prcc(X,Y(:,m));
end
end

function r = compute_prcc(X,y)
Xr = tiedrank(X);
yr = tiedrank(y);

p = size(X,2);
r = zeros(p,1);

for i = 1:p
    idx = setdiff(1:p,i);
    rx  = Xr(:,i) - Xr(:,idx)*(Xr(:,idx)\Xr(:,i));
    ry  = yr      - Xr(:,idx)*(Xr(:,idx)\yr);
    r(i) = corr(rx,ry);
end
end
