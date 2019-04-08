clear all;
load('100_randomized_results.mat')
load('Thai.mat')

pa = 3;  % 100_randomized_results.mat is run for individual #3

%% Plot two simulated results
%{
figure; 
loglog(pred_metabolome_all(:,1) / sum(pred_metabolome_all(:,1)),pred_metabolome_all(:,2) / sum(pred_metabolome_all(:,2)),'ko','Markersize',15)
hold on;
plot([1e-6, 1],[1e-6, 1],'k-')
%}


%% Plot the local sensitivity analysis of metabolome

norm = repmat(sum(pred_metabolome_all,1),2244,1);size(norm)
pred_metabolome_all = pred_metabolome_all ./ norm;

metabolome = Thai_metabolome_chia_full_norm(:,pa);
i_common = find(sign(mean(pred_metabolome_all,2)) .* sign(metabolome));

%{
%mean_pred_metabolome = mean(pred_metabolome_all(i_common,:),2);
%error_pred_metabolome = std(pred_metabolome_all(i_common,:),0,2);
figure;
%loglog(mean_pred_metabolome, metabolome(i_common), 'ko')
h=errorbar(metabolome(i_common),mean_pred_metabolome, error_pred_metabolome,'ko');
hold on;
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'YScale','log')
set(gca,'XScale','log')
ylim manual
h.LData = h.YData - max(1e-6,h.YData-h.LData);
%}

mean_pred_metabolome = mean(log10(pred_metabolome_all(i_common,:)+1e-6),2);
error_pred_metabolome = std(log10(pred_metabolome_all(i_common,:)+1e-6),0,2);

figure;
%loglog(mean_pred_metabolome, metabolome(i_common), 'ko')
h=errorbar(log10(metabolome(i_common)+1e-6),mean_pred_metabolome, error_pred_metabolome,'ko','Markersize',15);
hold on;
plot([-6, 0],[-6, 0],'k-')
hold on;
plot(log10(metabolome(i_common)+1e-6),log10(pred_metabolome_without_variation(i_common) / sum(pred_metabolome_without_variation(i_common))+1e-6),'rs','Markersize',15)
xlabel('Real metabolome')
ylabel('100 randomized predicted metabolome')
saveas(gcf,'./saved_Figures/SupplFig2A.svg')

% overpredicted microbes
% names_all(i_common(find(log10(metabolome(i_common)+1e-6)-mean_pred_metabolome < 0)))

%% Plot the local sensitivity analysis of metagenome

norm = repmat(sum(pred_metagenome_all,1),2244,1);size(norm)
pred_metagenome_all = pred_metagenome_all ./ norm;

b_real = Thai_abundance_chia_full_norm(:,pa);

%{
mean_pred_metagenome = mean(pred_metagenome_all,2);
error_pred_metagenome = std(pred_metagenome_all,0,2);

figure;
h=errorbar(b_real,mean_pred_metagenome, error_pred_metagenome,'ko');
hold on;
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'YScale','log')
set(gca,'XScale','log')
ylim manual
h.LData = h.YData - max(1e-6,h.YData-h.LData);
%}

mean_pred_metagenome = mean(log10(pred_metagenome_all+1e-6),2);
error_pred_metagenome = std(log10(pred_metagenome_all+1e-6),0,2);

figure;
h=errorbar(log10(b_real+1e-6),mean_pred_metagenome, error_pred_metagenome,'ko','Markersize',10);
hold on;
plot([-6, 0],[-6, 0],'k-')
hold on;
plot(log10(b_real+1e-6),log10(pred_metagenome_without_variation+1e-6),'rs','Markersize',10)
xlabel('Real metagenome')
ylabel('100 randomized predicted metagenome')
saveas(gcf,'./saved_Figures/SupplFig2B.svg')
