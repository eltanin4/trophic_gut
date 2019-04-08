%% Non-linear optimization for one individual (one child from Thai data)

clear all;

load('chia_network_new.mat')
load('Thai.mat')

pa = 20;   % specify which child data to be used
b_real = Thai_abundance_chia_full_norm(:,pa);  % load the real metagenomic 
                                               % data of the child pa.
numLayer_max = 4; % define maximum number of trophic layers
f_byproduct = 0.9;  % define "leakage" f

% i_selfish calculates all microbe ids which
% do not produce any byproduct. You can set their f to 0.
f = f_byproduct .* ones(2244,1);
f(i_selfish) = 0.0;


% b2m is a matrix determines the nutrient splitting among microbes
% m2b is a matrix determines the byproducts generation
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
options = optimoptions(@lsqnonlin,'Algorithm','Trust-region-reflective','Display','iter','TolFun',1e-4,'TolX',1e-4);
x = lsqnonlin(funct,x0,lb, ub, options);

% predicted metagenomics vs real metagenomics
x_full = zeros(2244,1);
x_full(i_intake) = x;
% m2b_total transform the nutrient intake into predicted bacterial abundance
predicted_b_a = m2b_total * x_full;    
corr(predicted_b_a, b_real)
figure;
plot(predicted_b_a + 1e-6, b_real + 1e-6, 'ko','Markersize',15)
hold on
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1e-6 1 1e-6 1])
xlabel('Predicted abundance from optimized diet','FontSize',15,'Fontweight','Bold')
ylabel('real abundance','FontSize',15,'Fontweight','Bold')

% predicted metabolome vs real metabolome
x_full = zeros(2244,1);
x_full(i_intake) = x;
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

% predicted metabolome vs real metabolome
metabolome = Thai_metabolome_chia_full_norm(:,pa);
i_common = find(sign(c_layer(:,4)) .* sign(metabolome));
figure;
plot(c_layer(i_common,4) / sum(c_layer(i_common,4)) + 1e-6, metabolome(i_common) / sum(metabolome(i_common)) + 1e-6, 'ko','Markersize',15)
hold on
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1e-6 1 1e-6 1])
xlabel('Predicted metabolome','FontSize',15,'Fontweight','Bold')
ylabel('real metabolome','FontSize',15,'Fontweight','Bold')

[a,b] = sort(c_layer(i_common,ii),'descend'); 
disp(['There are ',num2str(length(i_common)),' shared metabolites in metabolome:'])
disp(names_all(i_common(b)))
