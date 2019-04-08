% This script is used to generate the deviation results of 100 randomized
% results
clear all;
load('Thai.mat')
load('chia_network_new.mat')


%% Without variation
pa = 3;
b_real = Thai_abundance_chia_full_norm(:,pa);

numLayer_max = 4; % define maximum number of trophic layers
f_byproduct = 0.9;  % define "leakage" f
% i_selfish calculates all microbe ids which
% do not produce any byproduct. You can set their f to 0.
% b2m is saved in 'microbiome_PageRank.mat' which is a 
% matrix determines the nutrient splitting among microbes
[m2b, b2m] = Ain_out(b_real, i_all_filt, j_all_filt, v_all_filt);
out_degree = reshape(sum(full(b2m),1), [2244,1]);
i_selfish = i_b_all(find(out_degree(i_b_all) == 0));
f = f_byproduct .* ones(2244,1);
f(i_selfish) = 0.0;

% A_in, A_out matrix [c2b_real, bp_out] and [mu_matrix, ctot_matrix] is calculated
[m2b, b2m] = Ain_out(b_real, i_all_filt, j_all_filt, v_all_filt);
[m2b_total, m2m_total, m2m_layer] = mu(f, m2b, b2m, numLayer_max);

% This following part is used to remove all nutrients in the nutrient
% intake not used by the patient:
x_full_init = zeros(2244,1);
x_full_init(i_nutrient_intake) = 1;
i_nutrient_intake_used = find((sum(m2b,2) .* x_full_init));
i_intake = i_nutrient_intake_used;

lb = zeros(size(i_intake,1),1)'; % lb is lower bounds of variables
ub = ones(size(i_intake,1),1)' * 100; % ub is upper bounds of variables

x0 = zeros(length(i_intake),1);  % x0 is initial conditions of variables
x0(:) = 0.1;
x0 = x0 /sum(x0);

funct = @(x)pred_error_func(x,f, b_real, m2b_total, m2b, b2m,i_intake);
options = optimoptions(@lsqnonlin,'Algorithm','Trust-region-reflective','TolFun',1e-4,'TolX',1e-4);
x = lsqnonlin(funct,x0,lb, ub, options);

% predicted metagenomics vs real metagenomics
x_full = zeros(2244,1);
x_full(i_intake) = x;
predicted_b_a = m2b_total * x_full;

% predicted metabolome vs real metabolome
c_layer = zeros(2244,numLayer_max);
c_layer_not_used = zeros(2244,numLayer_max);
c_layer_not_used_now = zeros(2244,1);

for ii = 1:numLayer_max
    % m2m_layer transform the nutrient intake into nutrient at each layer
    c_layer(:,ii) = m2m_layer(:,:,ii) * x_full;
    i_layer_not_used = find(sign(c_layer(:,ii)) - sign(sum(m2b,2) .* c_layer(:,ii)));
    c_layer_not_used(i_layer_not_used,ii) = c_layer(i_layer_not_used,ii); 
    c_layer_not_used_now = c_layer_not_used_now + c_layer_not_used(:,ii);
    c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used_now;
end

pred_metabolome_without_variation = c_layer(:,4);
pred_metagenome_without_variation = predicted_b_a;


%% 100 random variation of uptake rates and byproduct production
pred_metabolome_all = zeros(2244, 100);
pred_metagenome_all = zeros(2244, 100);

for numRandom = 1:100
    disp(numRandom);
    %pa = 3;
    b_real = Thai_abundance_chia_full_norm(:,pa);

    %numLayer_max = 4; % define maximum number of trophic layers
    %f_byproduct = 0.9;  % define "leakage" f
    % i_selfish calculates all microbe ids which
    % do not produce any byproduct. You can set their f to 0.
    % b2m is saved in 'microbiome_PageRank.mat' which is a 
    % matrix determines the nutrient splitting among microbes
    [m2b, b2m] = Ain_out_randomization(b_real, i_all_filt, j_all_filt, v_all_filt);
    out_degree = reshape(sum(full(b2m),1), [2244,1]);
    i_selfish = i_b_all(find(out_degree(i_b_all) == 0));
    f = f_byproduct .* ones(2244,1);
    f(i_selfish) = 0.0;

    % A_in, A_out matrix [c2b_real, bp_out] and [mu_matrix, ctot_matrix] is calculated
    [m2b, b2m] = Ain_out_randomization(b_real, i_all_filt, j_all_filt, v_all_filt);
    [m2b_total, m2m_total, m2m_layer] = mu(f, m2b, b2m, numLayer_max);

    % This following part is used to remove all nutrients in the nutrient
    % intake not used by the patient:
    x_full_init = zeros(2244,1);
    x_full_init(i_nutrient_intake) = 1;
    i_nutrient_intake_used = find((sum(m2b,2) .* x_full_init));
    i_intake = i_nutrient_intake_used;


    lb = zeros(size(i_intake,1),1)'; % lb is lower bounds of variables
    ub = ones(size(i_intake,1),1)' * 100; % ub is upper bounds of variables

    x0 = zeros(length(i_intake),1);  % x0 is initial conditions of variables
    x0(:) = 0.1;
    x0 = x0 /sum(x0);

    funct = @(x)pred_error_func(x,f, b_real, m2b_total, m2b, b2m,i_intake);
    options = optimoptions(@lsqnonlin,'Algorithm','Trust-region-reflective','TolFun',1e-4,'TolX',1e-4);
    x = lsqnonlin(funct,x0,lb, ub, options);

    % predicted metagenomics vs real metagenomics
    x_full = zeros(2244,1);
    x_full(i_intake) = x;
    predicted_b_a = m2b_total * x_full;

    % predicted metabolome vs real metabolome
    c_layer = zeros(2244,numLayer_max);
    c_layer_not_used = zeros(2244,numLayer_max);
    c_layer_not_used_now = zeros(2244,1);

    for ii = 1:numLayer_max
        % m2m_layer transform the nutrient intake into nutrient at each layer
        c_layer(:,ii) = m2m_layer(:,:,ii) * x_full;
        i_layer_not_used = find(sign(c_layer(:,ii)) - sign(sum(m2b,2) .* c_layer(:,ii)));
        c_layer_not_used(i_layer_not_used,ii) = c_layer(i_layer_not_used,ii); 
        c_layer_not_used_now = c_layer_not_used_now + c_layer_not_used(:,ii);
        c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used_now;
    end
    
    pred_metabolome_all(:,numRandom) = c_layer(:,4);
    pred_metagenome_all(:,numRandom) = predicted_b_a;
end

save('100_randomized_results.mat','pred_metabolome_all','pred_metabolome_without_variation','pred_metagenome_all','pred_metagenome_without_variation')