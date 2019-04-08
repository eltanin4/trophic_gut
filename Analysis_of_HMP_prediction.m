clear all;
load('HMP_prediction.mat')
%load('HMP_prediction_old.mat')
load('chia_network_new.mat')

%% calculating the HMP data from Veronika's new mappings
% read Veronika's new mapping of microbes in 380 patients to Chia namespace
[a,b]=xlsread('./abundance_table_matched.xlsx');
abundance_table_full=a(:,3:end);
i1=find(a(:,1)); whos i1;
bugs_2_microbes_full=a(:,1);
for m=1:380; 
    a1=sparse(bugs_2_microbes_full(i1),ones(size(i1)), abundance_table_full(i1,m), 2244,1); 
    abundance_chia_full(1:2244,m)=a1; 
end;

%% Layer by layer activity analysis (metabolic and bacterial)
% metabolites produced at each layer
c1=reshape(sum(c_lbyl,1),4,380);
cdiets = reshape(intake_all(:,1,1,:),2244,380);
c1_norm=c1./repmat(sum(cdiets,1),4,1);

% metabolites leftover at each layer.
c2=reshape(sum(c_left,1),4,380);
c2_norm=c2./repmat(sum(cdiets,1),4,1);

b1= reshape(sum(b_lbyl,1),4,380);
b1_norm=b1(1:end-1,:)./repmat(sum(cdiets,1),3,1);

%figure; bar(sum(c1_norm,2)/380)
% hold on
% bar(sum(c2_norm,2)/380,'r')
% bar(sum(b1_norm,2)/380,'g')d

c1p = sum(c1_norm,2)/380;
c2p=sum(c2_norm,2)/380;
b1p=sum(b1_norm,2)/380;
b1p=vertcat([0],b1p);
c1pp=c1p-c2p;
%
c2p(4)=c2p(4)+c1pp(4)
c1pp(4)=0;
%
figure; 
bar(horzcat(c1pp,c2p,b1p), 'stacked'); legend('byproducts','leftover','bacteria');
saveas(gcf,'./saved_Figures/Fig3A.svg')
figure; 
%bar(horzcat(c1pp./sum(c1pp),c2p./sum(c2p),b1p./sum(b1p)), 'grouped'); legend('byproducts','leftover','bacteria');
bar(horzcat(c2p./sum(c2p),b1p./sum(b1p)), 'grouped'); legend('leftover','bacteria');
saveas(gcf,'./saved_Figures/Fig3B.svg')

%% Calculate layer by layer conversion and diversity for each layer
% calculate the layer by layer conversion and fixation of biomass
b_lbyl1 = b_lbyl(:,1:end-1,:);
c_lbyl1 = c_lbyl(:,2:end,:);
c_left1 = c_left(:,2:end,:);

c1=reshape(sum(c_lbyl,1),4,380);
%cdiets = reshape(intake_all(:,1,1,:),2244,380);

norm_factor = repmat(reshape(sum(cdiets,1), 1, 1, 380), [2244,3,1]);

b_norm = b_lbyl1./norm_factor;

c_left_norm = c_left1./norm_factor;
c_lbyl_norm = c_lbyl1./norm_factor;
c_left_norm(:,end,:) = c_lbyl_norm(:,end,:);

% alpha and gamma diversity for microbes
for pa=1:380; 
    b2p=b_norm(:,:,pa)./repmat(sum(b_norm(:,:,pa),1),2244,1,1); 
    alpha_D(1:3,pa)=1./sum(b2p.^2,1); 
    alpha_D_shannon(1:3,pa)=exp(-sum(b2p.*log(max(b2p,1e-10)),1)); 
%    
    [i,j,v]=find(b2p);
    b2p_genus=full(sparse(bact2genus(i), j, v, 241, 3));
    alpha_D_g(1:3,pa) = 1./sum(b2p_genus.^2,1);
    alpha_D_g_shannon(1:3,pa)=exp(-sum(b2p_genus.*log(max(b2p_genus,1e-10)),1));   
end;

b_norm_total = sum(b_norm,3);
b2p = b_norm_total ./ repmat(sum(b_norm_total,1),2244,1);
gamma_D(1:3) = 1./sum(b2p.^2,1);
gamma_D_shannon(1:3)=exp(-sum(b2p.*log(max(b2p,1e-10)),1)); 
%
[i,j,v]=find(b2p);
b2p_genus=full(sparse(bact2genus(i), j, v, 241, 3));
gamma_D_g(1:3) = 1./sum(b2p_genus.^2,1)';
gamma_D_g_shannon(1:3)=exp(-sum(b2p_genus.*log(max(b2p_genus,1e-10)),1))'; 
%
beta_D = gamma_D' ./ mean(alpha_D,2);
beta_D2 = gamma_D' .* mean(1./alpha_D,2);

beta_D_shannon = gamma_D_shannon' ./ mean(alpha_D_shannon,2);
beta_D2_shannon = gamma_D_shannon' .* mean(1./alpha_D_shannon,2);

beta_D_g = gamma_D_g' ./ mean(alpha_D_g,2)
beta_D2_g = gamma_D_g' .* mean(1./alpha_D_g,2)

beta_D_g_shannon = gamma_D_g_shannon'./mean(alpha_D_g_shannon,2)
beta_D2_g_shannon = gamma_D_g_shannon'.* mean(1./alpha_D_g_shannon,2)

for pa=1:380; 
    b2p=sum(b_norm(:,:,pa),2) ./ sum(sum(b_norm(:,:,pa),2));
    alpha_D_1l(pa)=1./sum(b2p.^2,1); 
end;
alpha_D2_mean=1./mean(1./alpha_D,2)
alpha_D2_mean_g=1./mean(1./alpha_D_g,2)
%figure; 
%bar(2:4,[alpha_D2_mean,beta_D2,gamma_D'./4])
% [AX,H1,H2] =plotyy(2:4,[alpha_D2_mean,beta_D2], 2:4,gamma_D', 'bar', 'bar');
% set(H1,'FaceColor','r') % a
% set(H2,'FaceColor','b') % b
% hLgnd=legend([H1(1) H1(2) H2(1)],'alpha','beta','gamma');


% alpha, beta and gamma diversity for diets
for pa=1:380; 
    b2p=c_lbyl(:,1,pa)./repmat(sum(c_lbyl(:,1,pa),1),2244,1); 
    alpha_D_diet(pa)=1./sum(b2p.^2,1); 
    alpha_D_shannon_diet(pa)=exp(-sum(b2p.*log(max(b2p,1e-10)),1)); 
end;

diet_mean = sum(c_lbyl(:,1,:),3);
diet_mean = diet_mean ./ sum(diet_mean);
gamma_D_diet = 1./sum(diet_mean.^2,1);
gamma_D_shannon_diet=exp(-sum(diet_mean.*log(max(diet_mean,1e-10)),1));

alpha_D2_diet_mean = 1 / mean(1./alpha_D_diet,2);
beta_D_diet = gamma_D_diet' ./ mean(alpha_D_diet,2);
beta_D2_diet = gamma_D_diet' .* mean(1./alpha_D_diet,2);

% alpha and gamma diversity for metabolites using all byproducts 
% leftovers and next layer inputs
for pa=1:380; 
    b2p=c_lbyl_norm(:,:,pa)./repmat(sum(c_lbyl_norm(:,:,pa),1),2244,1,1); 
    alpha_D_met(1:3,pa)=1./sum(b2p.^2,1); 
    alpha_D_shannon_met(1:3,pa)=exp(-sum(b2p.*log(max(b2p,1e-10)),1)); 
end;
alpha_D2_met_mean = 1 ./ mean(1./alpha_D_met,2);

b_norm_total = sum(c_lbyl_norm,3);
b2p = b_norm_total ./ repmat(sum(b_norm_total,1),2244,1);
gamma_D_met(1:3) = 1./sum(b2p.^2,1);
gamma_D_shannon_met(1:3)=exp(-sum(b2p.*log(max(b2p,1e-10)),1)); 
%
beta_D_met = gamma_D_met' ./ mean(alpha_D_met,2);
beta_D2_met = gamma_D_met' .* mean(1./alpha_D_met,2);
%beta_D2_met = gamma_D_met' ./ alpha_D2_met_mean;

beta_D_shannon_met = gamma_D_shannon_met' ./ mean(alpha_D_shannon_met,2);
beta_D2_shannon_met = gamma_D_shannon_met' .* mean(1./alpha_D_shannon_met,2);

for pa=1:380; 
    b2p=sum(c_lbyl_norm(:,:,pa),2) ./ sum(sum(c_lbyl_norm(:,:,pa),2));
    alpha_D_1l_met(pa)=1./sum(b2p.^2,1); 
end;


% Plot all diversity calculated above
alpha_all = vertcat([alpha_D2_diet_mean],alpha_D2_met_mean);
beta_all = vertcat([beta_D2_diet],beta_D2_met);
gamma_all = vertcat([gamma_D_diet],gamma_D_met');
beta_bact_all = vertcat([0],beta_D2_g);

%figure;
%bar(horzcat(alpha_all, beta_all, gamma_all, beta_bact_all))
%legend('\alpha diversity of nutrients','\beta diversity of nutrients','\gamma diversity of nutrients','\beta diversity of bacteria')

%figure;
%bar(horzcat(alpha_all, beta_all, gamma_all, beta_bact_all)')
%legend('\alpha diversity of nutrients','\beta diversity of nutrients','\gamma diversity of nutrients','\beta diversity of bacteria')

figure;
subplot(3,3,3)
bar(alpha_all,'g')
ylim([0, ceil(max(vertcat(alpha_D2_mean, alpha_all)))])
title('\alpha diversity of nutrients')

subplot(3,3,6)
bar(beta_all,'g')
ylim([0, ceil(max(vertcat(beta_D2, beta_all)))])
title('\beta diversity of nutrients')

subplot(3,3,9)
bar(gamma_all,'g')
ylim([0, ceil(max(vertcat(gamma_D', gamma_all)))])
title('\gamma diversity of nutrients')

subplot(3,3,1)
bar(alpha_D2_mean,'r')
ylim([0, ceil(max(vertcat(alpha_D2_mean, alpha_all)))])
title('\alpha diversity of bacteria')

subplot(3,3,4)
bar(beta_D2,'r')
ylim([0, ceil(max(vertcat(beta_D2, beta_all)))])
title('\beta diversity of bacteria')

subplot(3,3,7)
bar(gamma_D,'r')
ylim([0, ceil(max(vertcat(gamma_D', gamma_all)))])
title('\gamma diversity of bacteria')

subplot(3,3,2)
bar(alpha_D2_mean_g)
ylim([0, ceil(max(vertcat(alpha_D2_mean, alpha_all)))])
title('\alpha diversity of bacterial genus')

subplot(3,3,5)
bar(beta_D2_g)
ylim([0, ceil(max(vertcat(beta_D2, beta_all)))])
title('\beta diversity of bacterial genus')

subplot(3,3,8)
bar(gamma_D_g)
ylim([0, ceil(max(vertcat(gamma_D', gamma_all)))])
title('\gamma diversity of bacterial genus')

saveas(gcf,'./saved_Figures/Fig4.svg')