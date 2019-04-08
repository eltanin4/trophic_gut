% 11/9/2017 FULL NETWORK FROM CHIA ET AL SI
a=dlmread('microbe_to_compound_network_new_full.txt');
%% NEXT LINE IS WRONG AS IT CLUMPS PAIRS WITH THE SAME a(:,1) AND a(:,2) TOGETHER
% all2all_net=sparse(a(:,1),a(:,2),a(:,3),2244,2244);
% BETTER USE:
i_all=a(:,1);
j_all=a(:,2);
v_all=a(:,3);
%% filter currency metabolites with little nutritional value
i_remove=[2034  2039  2054  2055  2057  2113  2124  2168 2169  2176  2177  2244];
%     {'Bicarbonate (HCO3-, H2CO3, Carbonic acid, Carbonate)'}
%     {'Ca2+ (Calcium)'                                      }
%     {'Cl- (Chloride)'                                      }
%     {'CO2'                                                 }
%     {'Cu2+ (Copper)'                                       }
%     {'I- (Iodide)'                                         }
%     {'K+ (Potassium)'                                      }
%     {'Mg2+ (Magnesium)'                                    }
%     {'Mn2+ (Manganese)'                                    }
%     {'N2'                                                  }
%     {'Na+ (Sodium)'                                        }
%     {'Zn2+ (Zinc)'                                         }
a_keep=ones(2244,1);
a_keep(i_remove)=0;
i1=find(a_keep(i));
i_all_filt=i_all(i1);
j_all_filt=j_all(i1);
v_all_filt=v_all(i1);
% Description of interaction types by Veronika
% 2	 The metabolite is directly consumed (imported) by the microbe�
% 3	 The metabolite is directly produced (exported) by the microbe�
% 5	 The metabolite is both directly consumed (imported) and produced (exported) by the microbe�
%
% 6	 The macromolecule is subject to extracellular degradation by the microbe (this type of links I used for our subnetworks)
% 7	 The metabolite in the second column (where we usually have microbes) is derived from the extracellular degradation of the macromolecule in the first column.�
% 8	 The metabolite is indirectly produced (exported) by the microbe via extracellular macromolecule degradation 
%    (because we may have microbes which degrade polysaccharides, but not use some of by-products).
%%
% Advise on which compounds to remove from Nick Chia SI
% Step 1. Remove metabolites that are in large abundances from unspecified 
% intestinal sources (e.g., sloughed epithelial cells and secreted 
% proteins), because these metabolites would not be limiting factors in 
% metabolic interactions between microbes and the host. These metabolites
% include metal ions, CO2, N2, and bicarbonate.
% compounds_to_remove = [2034  2039  2054  2055  2057  2113  2124  2168
% 2169  2176  2177  2244];
[i,j,v]=find(all2all_net);
%%
analyze links of type 6, 7, 8
i1=find(v==6); whos i1; s1=sparse(i(i1),j(i1),ones(size(i1)),2244,2244);
i1=find(v==7); whos i1; s2=sparse(i(i1),j(i1),ones(size(i1)),2244,2244);
s3=s2'*s1; sum(sum(s3))
max(max(s3))
sum(sum(s3>0))
i1=find(v==8); whos i1; s4=sparse(i(i1),j(i1),ones(size(i1)),2244,2244);
sum(sum(s4.*(s3>0)))
c2m_ec=sign(s3);
% we found that matrix s3 has all 717 elements of type 8.
% Hence we will ignore type 8. 
%
%%
% We use the data for 381 patients compiled by Veronika
% from HMP, MetaHit study for Chinese healthy subjects
% all of the data come via WGS
%%
%import and map abundance data
[a,b]=xlsread('microbe_abundance_healthy_wgs.xlsx');
patients_new=b(1,:);
abundance_table_new=a(2:end,:);
bugs_new=b(2:end,1);
%%
% species prevalence and mean abundance before mapping to Chia et al 
species_prevalence_new=sum(abundance_table_new>0,2);
[a,b]=hist(species_prevalence_new./381, 0:0.1:1);
figure; semilogy(b,a,'ko-')
i2=find(species_prevalence_new); 
species_mean_abundance=zeros(720,1);
for m=1:length(i2); 
    i1=find(abundance_table_new(i2(m),:)); 
    species_mean_abundance(i2(m))=mean(abundance_table_new(i2(m),i1)); 
end;
figure; 
loglog(species_mean_abundance, species_prevalence_new,'ko')
[x y]=corr(species_mean_abundance(i2), species_prevalence_new(i2), 'type','Spearman')
i3=find(species_mean_abundance>0.1);
hold on;
[x1b,x2b,x2b_std,x2stat]=log_bin(species_mean_abundance(i2), species_prevalence_new(i2),2);
plot(x1b,x2b,'rd-')
%%
% now we map bugs to "microbes" in Chia et al network
bugs_2_microbes=zeros(size(bugs_new));
i1=find(a_nonempty_names_all);
l_i1=length(i1);
m_i1=max(y(z((l_i1+1):end)));
%sum(y(z((l_i1+1):end))<=814)
i1((l_i1+1):m_i1)=0;
bugs_2_microbes=i1(y(z((l_i1+1):end)));
%sum(bugs_2_microbes>0)
%names_all(bugs_2_microbes(find(bugs_2_microbes>0)))
%%
i1=find(bugs_2_microbes);
expl=sum(abundance_table_new(i1,:),1);
max(expl)
min(expl)
figure; hist(expl)
sum(expl>50)
sum(expl>50)./length(expl)
%%
% select the best patient to study
n_b_new=sum(abundance_table_new>0); 
i1=find(bugs_2_microbes);
n_b_new_mapped=sum(abundance_table_new(i1,:)>0); whos n*
mean(n_b_new)
mean(n_b_new_mapped)
figure; hist(n_b_new_mapped)
[a,b]=max(n_b_new_mapped)
expl(b)
[a1,b1]=max(expl)
n_b_new_mapped(b1)
%%
%patient=203; % is weird as it is 
% dominated by most abundant species has 24 mapped bugs explaining 99.3% of the abundance
patient=27; % 65 mapped species (45% total) and ~60% abundance explained
b_a=zeros(2244,1);
b_a(bugs_2_microbes(i_bugs_mapped))=abundance_table_new(i_bugs_mapped,patient);
b_a=b_a./sum(b_a);
i_b=find(b_a);c_
%%
%sum(c2m_in_norm,2)
%
% to deal with externally processed metabolites (such as fiber)
% for which we have information on which byproducts come from which inputs
% we may consider multiple "replicas" of organism: one per each input
% TO BE IMPLEMENTED LATER
[i,j,v]=find(all2all_net);
%%
%analyze links of type 6, 7, 8
i1=find(v==6); whos i1; s1=sparse(i(i1),j(i1),ones(size(i1)),2244,2244);
i1=find(v==7); whos i1; s2=sparse(i(i1),j(i1),ones(size(i1)),2244,2244);
s3=s2'*s1; sum(sum(s3))
max(max(s3))
sum(sum(s3>0))
i1=find(v==8); whos i1; s4=sparse(i(i1),j(i1),ones(size(i1)),2244,2244);
sum(sum(s4.*(s3>0)))
c2m_ec=sign(s3);
%%
% Description of interaction types by Veronika
% 2	 The metabolite is directly consumed (imported) by the microbe
% 3	 The metabolite is directly produced (exported) by the microbe
% 5	 The metabolite is both directly consumed (imported) and produced (exported) by the microbe
%
% 6	 The macromolecule is subject to extracellular degradation by the microbe (this type of links I used for our subnetworks)
% 7	 The metabolite in the second column (where we usually have microbes) is derived from the extracellular degradation of the macromolecule in the first column.�
% 8	 The metabolite is indirectly produced (exported) by the microbe via extracellular macromolecule degradation 
%    (because we may have microbes which degrade polysaccharides, but not use some of by-products).
%
b_a_real=b_a;
%%
b_a=rand(2244,1).*sign(b_a);
b_a=b_a./sum(b_a);
%%
[i,j,v]=find(all2all_net);
i1=find((v==2)+(v==5)+(v==6)); length(i1)
c2m_in=sparse(i(i1),j(i1),ones(size(i1)),2244,2244); 
c2m_in_norm=sparse(i(i1),j(i1),b_a(j(i1)),2244,2244);
%
i1=find((v==3)+(v==5)+(v==8)); length(i1)
c2m_out=sparse(i(i1),j(i1),ones(size(i1)),2244,2244);
c_out=sum(c2m_in_norm,2);
sum(c_out>0)
sum(sum(c2m_in_norm>0))
sum(sum(c2m_in>0))
[i2,j2,v2]=find(c2m_in_norm);
c2m_in_norm=sparse(i2,j2,v2./c_out(i2),2244,2244);
%
m_out=full(sum(c2m_out, 1)');
[i2,j2,v2]=find(c2m_out);
c2m_out_norm=sparse(i2,j2,sign(b_a(j2)).*v2./m_out(j2),2244,2244);
%%
%iterative loop to implement PageRank
% s3 is the nutrients-to-nutrient matrix taking through one step of
% byproduct generation
s3=c2m_in_norm*c2m_out_norm';
i1=find(sum(s3+s3'))'; whos i1
i2=find(a_polysach_diet(i1)>0);
ad=zeros(2244,1);
ad(i1(i2))=1;
sum(ad)
sum(ad.*a_polysach_diet)
sum(ad(i1))
%%
f=0.5;
a2o=ad; 
for m=1:10; 
    a2=(a2o'*s3)'; 
    a2=(1-f).*ad+f.*a2; 
%    disp(sum(a2)); 
%    disp(max(abs((a2o-a2)./a2))); 
    a2o=a2; 
end;
% [x,y]=sort(a2,'descend');
% sum(a2>0)
% x(1:79)'
% names_all(y(1:79))
b_pred=(a2'*c2m_in_norm)';
b_pred=b_pred./sum(b_pred);
hold on; loglog(b_a,b_pred,'o');
%max(abs((b_pred-b_a)./max(b_a, b_pred)))
b_a=b_pred;
%
%i_a=find(b_a_real); 
%[x y]=corr(b_a_real(i_a), b_pred(i_a), 'type', 'Spearman')
%%
% to iterate bacterial abundances use
%% b_a=b_pred;
%%
% clear the network from metabolites that cannot be used as food
% "currency metabolites"
a_keep=ones(2244,1);
a_keep(i_remove)=0;
[i,j,v]=find(all2all_net);
i1=find(a_keep(i)); whos i i1
all2all_net_filt=sparse(i(i1),j(i1),v(i1),2244,2244);
%%
% bacterial abundances in Chia namespace
i1=find(bugs_2_microbes>0); whos i1
species_chia_prevalence=zeros(2244,1);
species_chia_prevalence(bugs_2_microbes(i1))=species_prevalence_new(i1);
species_chia_mean_abundance=zeros(2244,1);
species_chia_mean_abundance(bugs_2_microbes(i1))=species_mean_abundance(i1);
species_chia_mean_abundance_all_samples=zeros(2244,1);
species_chia_mean_abundance_all_samples(bugs_2_microbes(i1))=mean(abundance_table_new(i1,:)')';
abundance_chia_table(bugs_2_microbes(i1),:)=abundance_table_new(i1,:);
%%
i2=find(b_pred>1e-6);
figure; loglog(species_chia_mean_abundance(i2), b_pred(i2), 'ko');
[a,b]=corr(species_chia_mean_abundance(i2), b_pred(i2),'Type', 'Spearman')
%%
% we varied diets by selecting 10 random nutrients in 
% equal abundance and calculated correlations between
% abundances with different cutoffs
p_pred_uniform=b_pred_uniform./sum(b_pred_uniform);
co_array=[1e-7,1e-6,1e-5,1e-4,1e-3,1e-2];
i_c_all=find(a_nonempty_names_all(2001:end))+2000;
tic;
clear corr_array corr_pvalue_array diet_array
for uuu=1:200;
    r1=randperm(244);
    diet_array(uuu,1:10)=i_c_all(r1(1:10));
    c_diet=zeros(2244,1);
    c_diet(i_c_all(r1(1:10)))=ones(10,1);
%    
    %c_diet=a_polysach_diet;
    %c_diet=a_select_diet;
%
    mac_arthur_loop;
%    
    p_pred=b_pred./sum(b_pred);
    for bbb=1:length(co_array);
        co=co_array(bbb);
        i5=find((species_chia_mean_abundance>co).*(p_pred>co));
        if isempty(i5)~=1;
            [a,b]=corr(species_chia_mean_abundance(i5),p_pred(i5),'type','Spearman');
        else;
            a=0; b=1;
        end;
%         i5=find((p_pred_uniform>co).*(p_pred>co));
%         [a,b]=corr(p_pred_uniform(i5),p_pred(i5),'type','Spearman');
        corr_array(uuu,bbb)=a;
        corr_pvalue_array(uuu,bbb)=b;
        %    disp([co,a,b]);
    end;
%     [a,b]=corr(p_pred_uniform,p_pred,'type','Pearson');
%     corr_array(uuu,bbb+1)=a;
%     corr_pvalue_array(uuu,bbb+1)=b;
end;
toc;
%%
c_diet=a_select_diet+10.*a_polysach_diet;
% add mucin 2170 to diet
c_diet(2170)=10;
c_diet=c_diet./sum(c_diet);
mac_arthur_loop;
%%
% we varied diets by selecting 10 random nutrients in 
% equal abundance and calculated correlations between
% abundances with different cutoffs
c_diet=a_select_diet+10.*a_polysach_diet;
% add mucin 2170 to diet
c_diet(2170)=10;
c_diet=c_diet./sum(c_diet);
%
co_array=[1e-7,1e-6,1e-5,1e-4,1e-3,1e-2];
tic;
clear corr_array corr_pvalue_array lambda_array
for uuu=1:200;
    mac_arthur_loop;
    lambda_array(uu,1:length(lambda_in_vec))=lambda_in_vec;
%    
    p_pred=b_pred./sum(b_pred);
    for bbb=1:length(co_array);
        co=co_array(bbb);
%         i5=find((species_chia_mean_abundance>co).*(p_pred>co));
%         if isempty(i5)~=1;
%         [a,b]=corr(species_chia_mean_abundance(i5),p_pred(i5),'type','Spearman');
%         else; 
%             a=0; b=1;
%         end;  
        i5=find((p_pred_tradeoff>co).*(p_pred>co));
        [a,b]=corr(p_pred_tradeoff(i5),p_pred(i5),'type','Spearman');
        corr_array(uuu,bbb)=a;
        corr_pvalue_array(uuu,bbb)=b;
        %    disp([co,a,b]);
    end;
%     [a,b]=corr(p_pred_tradeoff,p_pred,'type','Pearson');
%     corr_array(uuu,bbb+1)=a;
%     corr_pvalue_array(uuu,bbb+1)=b;
end;
toc;
%%
% clear the network from metabolites that cannot be used as food
% "currency metabolites"
a_keep=ones(2244,1);
a_keep(i_remove)=0;
[i,j,v]=find(all2all_net);
i1=find(a_keep(i)); whos i i1
all2all_net_filt=sparse(i(i1),j(i1),v(i1),2244,2244);
%%
% Metropolis prediction of all patients diet
%clear c_diet_all_patients b_pred_all_patients cc_all_patients cc_pvalue_all_patients
tic;
%for sp=1:100;
for sp=27:27;
    b_a_real=abundance_chia_full(:,sp)./sum(abundance_chia_full(:,sp));
    i_b=find(b_a_real);
    n_b=length(i_b);
    % Letting all diet components free.
    % i_diet_all=find(a_nonempty_names_all(2001:end))+2000;
    % a_diet_all=zeros(2244,1);
    % a_diet_all(i_diet_all)=1;
    % n_diet_all=length(i_diet_all);
    % c_diet=rand(2244,1).*a_diet_all;
    % c_diet=c_diet./sum(c_diet);
    % diet consisting of 10 polysaccharides + mucin at ~90%s,
    % ~12 of sugars at ~%10
    c_diet=a_select_diet+10.*a_polysach_diet;
%    c_diet=10.*a_polysach_diet;
    c_diet_original = c_diet;
    % add mucin 2170 to diet
    c_diet(2170)=10;
    c_diet=c_diet./sum(c_diet);
%
    i_diet=find(c_diet);
    n_diet=length(i_diet);
    c_diet=rand(2244,1).*sign(c_diet);
    c_diet=c_diet./sum(c_diet);
%   
    b_a=sign(b_a_real);
    b_a=rand(2244,1).*b_a;
    b_a=b_a./sum(b_a);
    mac_arthur_loop;
%    a=sum((log10(b_a_real(i_b))-log10(max(b_pred(i_b),1e-6))).^2)./n_b;
%    a=sum((b_a_real(i_b)-b_pred(i_b)).^2./(b_a_real(i_b)+b_pred(i_b)).^2)./n_b;
    [a,b]=corr(b_a_real(i_b),b_pred(i_b),'type','Spearman');
    corr_curr=a;
    pvalue_curr=b;
    i_s=find(b_a_real>1e-6);
    lagrange=4;
    for t=1:900;
        c_diet_old=c_diet;
        i_c=floor(n_diet.*rand)+1;
        %    c_diet(i_diet(i_c))=c_diet(i_diet(i_c)).*(0.1+2.*rand);
        c_diet(i_diet)=c_diet(i_diet).*(0.01+8.*rand(n_diet,1));
        c_diet=c_diet./sum(c_diet);
        b_a=sign(b_a_real);
        b_a=rand(2244,1).*b_a;
        b_a=b_a./sum(b_a);
        mac_arthur_loop;
 %       a=sum((b_a_real(i_b)-b_pred(i_b)).^2./(b_a_real(i_b)+b_pred(i_b)).^2)./n_b;
 %       a=sum((log10(b_a_real(i_b))-log10(max(b_pred(i_b),1e-6))).^2)./n_b;
        i_s=find(b_a_real>1e-6);
        [a,b]=corr(log10(b_a_real(i_s)),log10(max(b_pred(i_s),1e-6)),'type','Pearson');
        a=a+lagrange.*length(i_s)./n_b;
        if (a-corr_curr)./0.0001>log(rand);
%         if a<corr_curr;
            corr_curr=a;
            corr_true=corr_curr-lagrange.*length(i_s)./n_b;
            pvalue_curr=b;
            b_pred_opt=b_pred;
        else;
            c_diet_tried=c_diet;
            c_diet=c_diet_old;
        end;
        disp([t, length(i_s), a, corr_curr, corr_true, log10(pvalue_curr)]);
%        disp([t, a, corr_curr]);
    end;   
%     c_diet_all_patients(1:2244,sp)=c_diet;
%     c_tot_all_patients(1:2244,sp)=c_tot;
%     b_pred_all_patients(1:2244,sp)=b_pred_opt;
%     cc_all_patients(sp)=corr_curr;
%     cc_pvalue_all_patients(sp)=pvalue_curr;
    disp([sp, length(i_b), corr_curr, log10(pvalue_curr)]);
%     disp([sp, n_b, corr_curr]);
end;
toc;fi
%%
i_diet=find(c_diet); whos i_diet
cd=c_diet_all_patients(i_diet,:);
% clustering plot
CGobj = clustergram(cd, 'Standardize','Column', 'RowLabels',names_all(i_diet),'ColumnLabels',patients_new);
% PCA analysis
[pca_coeff, pca_score, pca_latent, pca_tsquared, pca_explained] = pca(cd');
figure; plot( pca_score(:,1),  pca_score(:,2),'ko');
%%
figure; plot(n_chia_species, cc_all_patients,'ko')
figure; loglog(sort(c_diet_all_patients,'descend'),'-')
diets_tsne = tsne(cd'); % only works onSergei's laptop
figure; plot(diets_tsne(:,1),diets_tsne(:,2),'ko')
patients_chinese=(1:157);
patients_american=(158:295);
patients_european=(296:380);
hold on;
plot(diets_tsne(patients_chinese,1),diets_tsne(patients_chinese,2),'rd')
plot(diets_tsne(patients_american,1),diets_tsne(patients_american,2),'bs')
plot(diets_tsne(patients_european,1),diets_tsne(patients_european,2),'bs')
%
figure; 
plot(pca_score(patients_chinese,1),  pca_score(patients_chinese,2),'rd')
hold on;
plot(pca_score(patients_american,1),  pca_score(patients_american,2),'bs')
plot(pca_score(patients_european,1),  pca_score(patients_european,2),'go')
i1=find((pca_score(:,1)<-0.3).*(pca_score(:,2)<-0.3)); whos i1
% i1 is a cluster of 6 Chinese and one american patient whose predicted
% diets are rich in Hyaluronan
%%
for patient=27:27;
% Working with layers, searching for optimal f and understanding how deep 
% a typical gut microbiome is for a given patient according to the 
% Chia network of utilization and production.
% IMPORTANT YOU NEED TO SET: f_byproduct=0.1
%c_diet=a_select_diet+10.*a_polysach_diet;
c_diet=10.*a_polysach_diet;
% add mucin 2170 to diet
c_diet(2170)=10;
%c_diet(2170)=1;
c_diet=c_diet./sum(c_diet);
b_a_real=abundance_chia_table(:,patient)./sum(abundance_chia_table(:,patient));
i_b=find(b_a_real);
b_a=b_a_real;
% b_a=sign(b_a_real);
% b_a=rand(2244,1).*b_a;
b_a=b_a./sum(b_a);
mac_arthur_loop;
i_b_present=find(b_pred);
i_c_present=find(c_pred);
%
% Plotting now
%figure; semilogy(sum(c_layer,1),'ko-')
% a3=sum(c_layer,1);
% sum(a3.*0.5.^(0:19))
% figure; pcolor(log10(c_layer(i_c_all,:))); shading flat; 
% caxis([-6,0]); colorbar;
%figure; semilogy(c_layer(i_c_present,:)','s-')
%legend(names_all(nonzeros(i_c_present)));
%
% figure; pcolor(log10(b_layer(i_b_present,:))); shading flat; 
% caxis([-6,0]); colorbar;
% figure; semilogy(b_layer(i_b_present,:)','o-')
% legend(names_all(i_b_present))
%%
% assigning layers to finally present bacteria
%[b_layer_max, b_assigned_layer]=max(b_layer');
[b_layer_max, b_assigned_layer]=max(sign(b_layer'));
b_assigned_layer= b_assigned_layer'.*sign(b_pred);
max(b_assigned_layer); % typical max is 3.
for ii=1:max(b_assigned_layer);
    disp(['layer ' num2str(ii)]);
    disp([names_all(find(b_assigned_layer==ii))]);
end;

%[c_layer_max, c_assigned_layer]=max(c_layer');
[c_layer_max, c_assigned_layer]=max(sign(c_layer'));
c_assigned_layer= c_assigned_layer'.*sign(c_tot);
% max(c_assigned_layer); % typical max is 3.
% for ii=1:max(c_assigned_layer);
%     disp(['layer ' num2str(ii)]);
%     disp([names_all(find(c_assigned_layer==ii))]);
% end;

% figure; semilogy(assigned_layer(i_b_present), b_pred(i_b_present), 'go') % resulting graph with b_pred
                                               % shows no clear trend of decreasing
% figure; semilogy(b_assigned_layer(i_b_present), b_a_real(i_b_present), 'ro')
[x1b,x2b,x2b_std,x2stat]=log_bin(b_assigned_layer(i_b_present),b_a_real(i_b_present),5);
% hold on; plot(x1b, x2b, 'rd-');
% figure; plot(x1b,log(x2b),'ko');
a=fit(x1b(1:2), x2b(1:2),'exp1');
f_fitted(patient)=exp(a.b);
disp([patient, exp(a.b)]);
end;

%%
% analysis of predicted abundances
hold on; semilogy(b_assigned_layer(i_b_present), b_pred_all_patients(i_b_present, 27), 'bs')
% According to this, f=0.9 fits better to the observed data for patient 27.
%%
% read Veronika's new mapping of microbes in 380 patients to Chia namespace
[a,b]=xlsread('abundance_table_matched.xlsx');
abundance_table_full=a(:,3:end);
i1=find(a(:,1)); whos i1;
bugs_2_microbes_full=a(:,1);
for m=1:380; 
    a1=sparse(bugs_2_microbes_full(i1),ones(size(i1)), abundance_table_full(i1,m), 2244,1); 
    abundance_chia_full(1:2244,m)=a1; 
end;
sum(sum(abundance_table_full(i1,:)))
sum(sum(abundance_chia_full))
abundance_chia_full=abundance_chia_full./100;
%%
% using SVD to find the least square solution
sp=27;     
b_a_real=abundance_chia_full(:,sp)./sum(abundance_chia_full(:,sp));
b_a=b_a_real;
i_b=find(b_a);
%%
% use for Thai data in Thai.mat
sp=3;
b_a_real=Thai_abundance_chia_full(:,sp)./sum(Thai_abundance_chia_full(:,sp));
b_a=b_a_real;
i_b=find(b_a);
% make sure mac_arthur loop does not iterate b_pred
mac_arthur_loop;
[x y]=sort(c_tot,'descend');
log10(x(1:92))
sum(x>1e-3)
% For patient 27 
% there are 83 c_tot which are larger than 1e-3
% then big drop to ~1e-6
i_diet_svd=y(1:83);
%
[U,S,V] = svd(full(c2b(i_diet_svd,i_b))');
whos U S V
sum(diag(S>1e-15))
%there are 50 non-zero eigenvalues
r1=1:sum(diag(S>1e-15));
d1=U(:,r1)*S(r1,r1)*V(:,r1)';
max(max(d1-full(c2b(i_diet_svd,i_b))'))
min(min(d1-full(c2b(i_diet_svd,i_b))'))
%d1 is equal to full(c2b(i_diet_svd,i_b))'
[i j v]=find(S(r1,r1));
Sinv=full(sparse(i,j,1./v, length(r1), length(r1)));
d_inv=U(:,r1)*Sinv*V(:,r1)';
c1=(b_a(i_b)'*d_inv)';
c_met=zeros(2244,1);
c_met(i_diet_svd)=c1;
b_pred=(c_met'*c2b)';
figure; loglog(b_a,b_pred,'kx')
hold on;
b_pred1=(max(c_met,0)'*c2b)';
loglog(b_a,b_pred1,'rd')
% loglog(b_a,b_pred1,'rd')
% figure; semilogy(sort(c1,'descend'),'ko-')
sum(c1>1e-2)
c_met_th=c_met.*(c_met>(3e-2));
b_pred2=(c_met_th'*c2b)';
loglog(b_a,b_pred2,'bd')
%%
% Thai data analysis
Thai_k_in=zeros(2244,41);
Thai_k_out=zeros(2244,41);
Thai_k_out_w=zeros(2244,41);
for sp=1:41;
    Thai_abundance_explained(sp)=sum(Thai_abundance_chia_full(:,sp))./5000;
    b_a_real=Thai_abundance_chia_full(:,sp)./5000;
    %sum(Thai_abundance_chia_full(:,sp));
    b_a=b_a_real;
    mac_arthur_loop;
    k_in=sum(sign(c2b),2);
    k_out=sum(sign(c2b_out),2);
    k_out_w=sum(c2b_out,2);
    Thai_k_in(:,sp)=k_in;
    Thai_k_out(:,sp)=k_out;
    Thai_k_out_w(:,sp)=k_out_w;
    Thai_metabolome_frac_consumable(sp)=sum((k_in>0).*Thai_metabolome_chia_full(:,sp))./sum(Thai_metabolome_chia_full(:,sp));
    Thai_metabolome_n_met_consumable(sp)=sum((k_in>0).*(Thai_metabolome_chia_full(:,sp)>0));
    Thai_metabolome_total_1(sp)=sum((k_in>0).*Thai_metabolome_chia_full(:,sp));
    Thai_metabolome_total_2(sp)=sum((k_in==0).*Thai_metabolome_chia_full(:,sp));
end;
%
ib=find(sum(Thai_abundance_chia_full_norm,2)); whos ib
im=find(sum(Thai_metabolome_chia_full,2)); whos im
CGobj = clustergram(Thai_metabolome_chia_full(im,:), 'Standardize','Row', 'RowPDist', 'correlation', 'ColumnPDist', 'correlation',  'RowLabels',names_all(im),'ColumnLabels',Thai_names_all');
CGobj = clustergram(Thai_abundance_chia_full_norm(ib,:), 'Standardize','Row', 'RowPDist', 'correlation', 'ColumnPDist', 'correlation',  'RowLabels',names_all(ib),'ColumnLabels',Thai_names_all');