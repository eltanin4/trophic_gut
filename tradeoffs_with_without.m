%% Non-linear optimization for one individual 3

%% with tradeoffs
pa = 3;
b_real = Thai_abundance_chia_full_norm(:,pa);

numLayer_max = 4; % define maximum number of trophic layers
f_byproduct = 0.9;  % define "leakage" f
% i_selfish calculates all microbe ids which
% do not produce any byproduct. You can set their f to 0.
% b2m is saved in 'microbiome_PageRank.mat' which is a 
% matrix determines the nutrient splitting among microbes
out_degree = reshape(sum(full(b2m),1), [2244,1]);
i_selfish = i_b_all(find(out_degree(i_b_all) == 0));
f = f_byproduct .* ones(2244,1);
f(i_selfish) = 0.0;

% A_in, A_out matrix [c2b_real, bp_out] and [mu_matrix, ctot_matrix] is calculated
[m2b, b2m] = Ain_out(b_real, i_all_filt, j_all_filt, v_all_filt);
[m2b_total, m2m_total, m2b_layer] = mu(f, m2b, b2m, numLayer_max);

% i_nutrient_intake containins labels of diets we want to optimize over
i_list  = unique(j_all_filt(find((v_all_filt==7) .* (i_all_filt~=2170)))); 
i_nutrient_intake = [i_list([3:(end-2),end]); i_sugar_diet([3,7])];

% This following part is used to remove all nutrients in the nutrient
% intake not used by the patient:
x_full_init = zeros(2244,1);
x_full_init(i_nutrient_intake) = 1;
i_nutrient_intake_used = find((sum(m2b,2) .* x_full_init));
i_intake_diet = i_nutrient_intake_used;


lb = zeros(size(i_fit_diet,1),1)'; % lb is lower bounds of variables
ub = ones(size(i_fit_diet,1),1)' * 100; % ub is upper bounds of variables

x0 = zeros(length(i_fit_diet),1);  % x0 is initial conditions of variables
x0(:) = 0.1;
x0 = x0 /sum(x0);

funct = @(x)pred_error_func(x,f, b_real, m2b_total, m2b, b2m,i_fit_diet);
%options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','iter','TolFun',1e-4,'TolX',1e-4);
options = optimoptions(@lsqnonlin,'Algorithm','Trust-region-reflective','Display','iter','TolFun',1e-4,'TolX',1e-4);
x = lsqnonlin(funct,x0,lb, ub, options);

% predicted metagenomics vs real metagenomics
predicted_error = pred_error_func(x,f, b_real, m2b_total, m2b, b2m,i_fit_diet);
%predicted_error(find(b_real))
x_full = zeros(2244,1);
x_full(i_fit_diet) = x;
predicted_b_a = m2b_total * x_full;
corr(predicted_b_a, b_real)
%[a,b]=corr(predicted_b_a, b_real,'Type','Spearman')
figure;
plot(predicted_b_a + 1e-6, b_real + 1e-6, 'ko')
hold on
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1e-6 1 1e-6 1])
xlabel('Predicted abundance from optimized diet','FontSize',15,'Fontweight','Bold')
ylabel('real abundance','FontSize',15,'Fontweight','Bold')

% predicted metabolome vs real metabolome
x_full = zeros(2244,1);
x_full(i_fit_diet) = x;
%[a,b] = sort(x_full(i_fit_diet),'descend');  
%x_full(i_fit_diet(b(5:end))) = 0;   % keep only top 4 metabolites
%predicted_b_a_layer = zeros(2244,2244);
c_layer = zeros(2244,numLayer_max);
c_layer_not_used = zeros(2244,numLayer_max);
c_layer_not_used_now = zeros(2244,1);
for ii = 1:numLayer_max
    %predicted_b_a_layer = predicted_b_a_layer + mu_layer(:,:,ii);
    %pred_error_layer = 
    c_layer(:,ii) = m2b_layer(:,:,ii) * x_full;
    i_layer_not_used = find(sign(c_layer(:,ii)) - sign(sum(m2b,2) .* c_layer(:,ii)));
    %disp(i_layer_not_used);
    c_layer_not_used(i_layer_not_used,ii) = c_layer(i_layer_not_used,ii); 
    c_layer_not_used_now = c_layer_not_used_now + c_layer_not_used(:,ii);
    c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used_now;
    %c_layer_not_used(i_layer_not_used,ii) = x_full(i_layer_not_used); 
    %c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used(:,ii);
end
metabolome = Thai_metabolome_chia_full_norm(:,pa);

%corr1 = zeros(numLayer_max,1);
corr2 = zeros(numLayer_max,1);

%predicted_b_a_layer = (1-f) * c2b_real' * predicted_b_a_layer * x_full;
%i_common_all = zeros(2244,uu_max);
for ii=1:numLayer_max
    i_common = find(sign(c_layer(:,ii)) .* sign(metabolome));
    %i_common_all(i_common,ii) = 1;
    if length(i_common) ~= 0
        %[corr1(ii),pval1(ii)] = corr(c_layer(i_common,ii), metabolome(i_common),'type','Spearman');
        [corr2(ii),pval2(ii)] = corr(c_layer(i_common,ii), metabolome(i_common),'type','Pearson');
    end
end
%figure;
%plot(corr1,'ko')
figure;
plot(corr2,'ko')

[a,b] = sort(c_layer(i_common,ii),'descend'); names_all(i_common(b))

predicted_b_a_with_tradeoffs = predicted_b_a;

%% without tradeoffs
pa = 3;
b_real = Thai_abundance_chia_full_norm(:,pa);

numLayer_max = 4; % define maximum number of trophic layers
f_byproduct = 0.9;  % define "leakage" f
% i_selfish calculates all microbe ids which
% do not produce any byproduct. You can set their f to 0.
% b2m is saved in 'microbiome_PageRank.mat' which is a 
% matrix determines the nutrient splitting among microbes
out_degree = reshape(sum(full(b2m),1), [2244,1]);
i_selfish = i_b_all(find(out_degree(i_b_all) == 0));
f = f_byproduct .* ones(2244,1);
f(i_selfish) = 0.0;

% A_in, A_out matrix [c2b_real, bp_out] and [mu_matrix, ctot_matrix] is calculated
[m2b, b2m] = Ain_out_without_tradeoffs(b_real, i_all_filt, j_all_filt, v_all_filt);
[m2b_total, m2m_total, m2b_layer] = mu(f, m2b, b2m, numLayer_max);

% i_nutrient_intake containins labels of diets we want to optimize over
i_list  = unique(j_all_filt(find((v_all_filt==7) .* (i_all_filt~=2170)))); 
i_nutrient_intake = [i_list([3:(end-2),end]); i_sugar_diet([3,7])];

% This following part is used to remove all nutrients in the nutrient
% intake not used by the patient:
x_full_init = zeros(2244,1);
x_full_init(i_nutrient_intake) = 1;
i_nutrient_intake_used = find((sum(m2b,2) .* x_full_init));
i_intake_diet = i_nutrient_intake_used;


lb = zeros(size(i_fit_diet,1),1)'; % lb is lower bounds of variables
ub = ones(size(i_fit_diet,1),1)' * 100; % ub is upper bounds of variables

x0 = zeros(length(i_fit_diet),1);  % x0 is initial conditions of variables
x0(:) = 0.1;
x0 = x0 /sum(x0);

funct = @(x)pred_error_func(x,f, b_real, m2b_total, m2b, b2m,i_fit_diet);
%options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','iter','TolFun',1e-4,'TolX',1e-4);
options = optimoptions(@lsqnonlin,'Algorithm','Trust-region-reflective','Display','iter','TolFun',1e-4,'TolX',1e-4);
x = lsqnonlin(funct,x0,lb, ub, options);

% predicted metagenomics vs real metagenomics
predicted_error = pred_error_func(x,f, b_real, m2b_total, m2b, b2m,i_fit_diet);
%predicted_error(find(b_real))
x_full = zeros(2244,1);
x_full(i_fit_diet) = x;
predicted_b_a = m2b_total * x_full;
corr(predicted_b_a, b_real)
%[a,b]=corr(predicted_b_a, b_real,'Type','Spearman')
figure;
plot(predicted_b_a + 1e-6, b_real + 1e-6, 'ko')
hold on
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1e-6 1 1e-6 1])
xlabel('Predicted abundance from optimized diet','FontSize',15,'Fontweight','Bold')
ylabel('real abundance','FontSize',15,'Fontweight','Bold')

% predicted metabolome vs real metabolome
x_full = zeros(2244,1);
x_full(i_fit_diet) = x;
%[a,b] = sort(x_full(i_fit_diet),'descend');  
%x_full(i_fit_diet(b(5:end))) = 0;   % keep only top 4 metabolites
%predicted_b_a_layer = zeros(2244,2244);
c_layer = zeros(2244,numLayer_max);
c_layer_not_used = zeros(2244,numLayer_max);
c_layer_not_used_now = zeros(2244,1);
for ii = 1:numLayer_max
    %predicted_b_a_layer = predicted_b_a_layer + mu_layer(:,:,ii);
    %pred_error_layer = 
    c_layer(:,ii) = m2b_layer(:,:,ii) * x_full;
    i_layer_not_used = find(sign(c_layer(:,ii)) - sign(sum(m2b,2) .* c_layer(:,ii)));
    %disp(i_layer_not_used);
    c_layer_not_used(i_layer_not_used,ii) = c_layer(i_layer_not_used,ii); 
    c_layer_not_used_now = c_layer_not_used_now + c_layer_not_used(:,ii);
    c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used_now;
    %c_layer_not_used(i_layer_not_used,ii) = x_full(i_layer_not_used); 
    %c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used(:,ii);
end
metabolome = Thai_metabolome_chia_full_norm(:,pa);

%corr1 = zeros(numLayer_max,1);
corr2 = zeros(numLayer_max,1);

%predicted_b_a_layer = (1-f) * c2b_real' * predicted_b_a_layer * x_full;
%i_common_all = zeros(2244,uu_max);
for ii=1:numLayer_max
    i_common = find(sign(c_layer(:,ii)) .* sign(metabolome));
    %i_common_all(i_common,ii) = 1;
    if length(i_common) ~= 0
        %[corr1(ii),pval1(ii)] = corr(c_layer(i_common,ii), metabolome(i_common),'type','Spearman');
        [corr2(ii),pval2(ii)] = corr(c_layer(i_common,ii), metabolome(i_common),'type','Pearson');
    end
end
%figure;
%plot(corr1,'ko')
figure;
plot(corr2,'ko')

[a,b] = sort(c_layer(i_common,ii),'descend'); names_all(i_common(b))

predicted_b_a_without_tradeoffs = predicted_b_a;

%%

in_degree = sum(sign(m2b),1);
figure;
semilogy(in_degree, predicted_b_a_with_tradeoffs,'ko')
hold on;
semilogy(in_degree, predicted_b_a_without_tradeoffs,'bs')
hold on;
semilogy(in_degree, b_real,'rd')
legend({'predicted data with trade-offs','predicted data without trade-offs','real data'})

corr(in_degree', predicted_b_a_with_tradeoffs,'type','Pearson')
corr(in_degree', predicted_b_a_without_tradeoffs,'type','Pearson')
corr(in_degree', b_real,'type','Pearson')