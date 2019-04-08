%load('digested_poly_diet.mat')
%load('phase_diagram.mat')
%load('lambda_random.mat')
clear all;
load('Calibration_of_the_model_new.mat')
load('Thai.mat')

%% Read calculate the the correlation

a = sum(Thai_metabolome_chia_full_norm,2);   % summed metabolites in real metabolome
b = sum(reshape(c_layer_all(:,5,3,:),2244,41),2);  % summed metabolites in predicted metabolome

i_all_metabolites = find(a .* b);   % find the intersection between real metabolome and predicted metabolome

% Recalculate the Pearson correlation by using the common shared set of
% metabolites in the metabolome:
for kk = 1:length(f_list)
    for ll = 1:length(numLayer_max_list)
        f_byproduct = f_list(kk);
        numLayer_max = numLayer_max_list(ll);
        for pa=1:41
            metabolome = Thai_metabolome_chia_full_norm(:,pa);
            [corrP,pvalP] = corr(c_layer_all(i_all_metabolites,kk,ll,pa), metabolome(i_all_metabolites),'type','Pearson');
            corrP_all(kk,ll,pa) = corrP;
        end
    end
    
end

for sp=1:41; 
%    dev_best(sp)=min(min(dev(:,:,sp))); 
    aux=corrP_all(:,:,sp);
    [i,j,v]=find(aux);
    [a,b]=max(v);
    ff_best(sp)=i(b);
    numLayer_best(sp)=j(b);
    dev_best(sp)=v(b);
    %i_common_best(sp)=i_common_array(i(b),j(b),sp);
end;


% Mean Pearson correlation coefficient between real metabolome and predicted metabolome
% visualization
aux=sum(corrP_all,3)/41; 
figure; 
imagesc(f_list,numLayer_max_list,aux');
set(gca, 'XTick', f_list, 'XTickLabel', f_list, 'YTick', [2:2:10], 'YTickLabel', numLayer_max_list,'fontsize',15,'fontweight','Bold')
xlabel('leakage fraction f')
ylabel('number of layers')
title('Average Pearson correlation coefficient')
colorbar();
shading flat;
saveas(gcf,'./saved_Figures/Fig2A.svg')

figure; 
histogram(f_list(ff_best),[0.0, f_list+0.1])
title('Histogram of best f for each individual')
set(gca, 'XTick', f_list)
saveas(gcf,'./saved_Figures/SupplFig1A.svg')

figure; 
histogram(numLayer_max_list(numLayer_best))
title('Histogram of number of layers for each individual')
saveas(gcf,'./saved_Figures/SupplFig1B.svg')



%% Scatter plots of bacterial metagenomics and metabolome of 
% one individual #3 (f = 0.9, number of layers = 4)
pa=20;
b_real = Thai_abundance_chia_full_norm(:,pa);

c_layer = c_layer_all(:,5,3,pa);
predicted_b_a = pred_abun_all(:,5,3,pa);

metabolome = Thai_metabolome_chia_full_norm(:,pa);
i_common = i_all_metabolites;
%i_common = find(sign(c_layer) .* sign(metabolome));

corr(predicted_b_a, b_real)
%[a,b]=corr(predicted_b_a, b_real,'Type','Spearman')
set(0,'defaultLineLineWidth',2);
figure;
plot(predicted_b_a + 1e-5, b_real + 1e-5, 'ko','Markersize', 15)
hold on
plot([1e-5, 1],[1e-5, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'Fontsize',16)
axis([1e-5 1 1e-5 1])
%xlabel('Predicted abundance','FontSize',15,'Fontweight','Bold')
%ylabel('real abundance','FontSize',15,'Fontweight','Bold')
xlabel('Predicted abundance','FontSize',20)
ylabel('real abundance','FontSize',20)
title('f = 0.9, number of layers = 4','FontSize',20,'Fontweight','Bold')
saveas(gcf,'./saved_Figures/Fig2B.svg')

figure;
plot(c_layer(i_common) / sum(c_layer(i_common)) + 1e-5, metabolome(i_common) / sum(metabolome(i_common)) + 1e-5, 'ko','Markersize', 15)
%plot(c_layer(:,numLayer_max) + 1e-6, metabolome + 1e-6, 'ko')
hold on
plot([1e-5, 1],[1e-5, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'Fontsize',16)
axis([1e-5 1 1e-5 1])
%xlabel('Predicted metabolome','FontSize',15,'Fontweight','Bold')
%ylabel('Real metabolome','FontSize',15,'Fontweight','Bold')
xlabel('Predicted metabolome','FontSize',20)
ylabel('Real metabolome','FontSize',20)
title('f = 0.9, number of layers = 4','FontSize',20,'Fontweight','Bold')
saveas(gcf,'./saved_Figures/Fig2C.svg')

% [cor, p] = corr(c_layer(i_common) / sum(c_layer(i_common)) + 1e-5, metabolome(i_common) / sum(metabolome(i_common)) + 1e-5)
% [cor, p] = corr(predicted_b_a, b_real)

%% mean metabolome and metagenomics - Scatter plots of the case (f = 0.9, number of layers = 4)
%pa=3;
b_real = sum(Thai_abundance_chia_full_norm,2) ./ max(sum(Thai_abundance_chia_full_norm>0,2), 1);

c_layer = sum(c_layer_all(:,5,3,:),4) ./ max(sum(c_layer_all(:,5,3,:)>0,4), 1);
predicted_b_a = sum(pred_abun_all(:,5,3,:),4) ./ max(sum(pred_abun_all(:,5,3,:)>0,4), 1);

metabolome = sum(Thai_metabolome_chia_full_norm,2) ./ max(sum(Thai_metabolome_chia_full_norm>0,2), 1);

i_common = find(sign(c_layer) .* sign(metabolome));

corr(predicted_b_a, b_real)
%[a,b]=corr(predicted_b_a, b_real,'Type','Spearman')
set(0,'defaultLineLineWidth',2);
figure;
plot(predicted_b_a + 1e-6, b_real + 1e-6, 'ko','Markersize', 15)
hold on
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1e-6 1 1e-6 1])
xlabel('Predicted abundance','FontSize',15,'Fontweight','Bold')
ylabel('real abundance','FontSize',15,'Fontweight','Bold')
title('f = 0.9, number of layers = 4','FontSize',15,'Fontweight','Bold')

figure;
plot(c_layer(i_common) / sum(c_layer(i_common)) + 1e-5, metabolome(i_common) / sum(metabolome(i_common)) + 1e-5, 'ko','Markersize', 15)
%plot(c_layer(:,numLayer_max) + 1e-6, metabolome + 1e-6, 'ko')
hold on
plot([1e-5, 1],[1e-5, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'Fontsize',18)
axis([1e-5 1 1e-5 1])
%xlabel('Predicted metabolome','FontSize',15,'Fontweight','Bold')
%ylabel('Real metabolome','FontSize',15,'Fontweight','Bold')
xlabel('Predicted mean metabolome','FontSize',25)
ylabel('Real mean metabolome','FontSize',25)
title('f = 0.9, number of layers = 4','FontSize',25,'Fontweight','Bold')


%% Scatter plots of mean metabolome of lambda = 1 vs mean metabolome of varied lambda (f = 0.9, number of layers = 4)
%{
%load('digested_poly_diet.mat')
c_layer = sum(c_layer_all(:,5,3,:),4) ./ max(sum(c_layer_all(:,5,3,:)>0,4), 1);
load('lambda_random.mat')
c_layer2 = sum(c_layer_all(:,5,3,:),4) ./ max(sum(c_layer_all(:,5,3,:)>0,4), 1);
figure;
plot(c_layer / sum(c_layer) + 1e-6, c_layer2 / sum(c_layer2) + 1e-6, 'ko','Markersize', 15)
%plot(c_layer + 1e-6, metabolome + 1e-6, 'ko')
hold on
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1e-6 1 1e-6 1])
xlabel('mean Predicted metabolome of \lambda = 1','FontSize',15,'Fontweight','Bold')
ylabel('mean Predicted metabolome of random \lambda','FontSize',15,'Fontweight','Bold')
title('f = 0.9, number of layers = 4','FontSize',15,'Fontweight','Bold')
%}