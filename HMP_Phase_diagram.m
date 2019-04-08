%% This function is used to generate the phase diagram of f and uu_max

%% calculating the HMP data from Veronika's new mappings
% read Veronika's new mapping of microbes in 380 patients to Chia namespace
[a,b]=xlsread('../abundance_table_matched.xlsx');
abundance_table_full=a(:,3:end);
i1=find(a(:,1)); whos i1;
bugs_2_microbes_full=a(:,1);
for m=1:380; 
    a1=sparse(bugs_2_microbes_full(i1),ones(size(i1)), abundance_table_full(i1,m), 2244,1); 
    abundance_chia_full(1:2244,m)=a1; 
end;


%%
% assuming that uu_max is the same for all patients.

%f_list = [0.1,0.3,0.5,0.7,0.9];
f_list = [0.9];
uu_max_list = [4];
% f_list = [0.1];
% uu_max_list = [2];

%load('Thai.mat')
load('microbiome_PageRank.mat')

%f1 = waitbar(0,'1','Name','Approximating pi...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

%uu_max = 4;
diet_all = zeros(2244,length(f_list),length(uu_max_list),380);
pred_abun = zeros(2244,length(f_list),length(uu_max_list),380);
corrS_all = zeros(length(f_list),length(uu_max_list),380);
corrP_all = zeros(length(f_list),length(uu_max_list),380);
pvalS_all = zeros(length(f_list),length(uu_max_list),380);
pvalP_all = zeros(length(f_list),length(uu_max_list),380);
c_layer_not_used_all = zeros(2244,length(f_list),length(uu_max_list),380);
c_layer_all = zeros(2244,length(f_list),length(uu_max_list),380);
b_pred_all = zeros(2244,length(f_list),length(uu_max_list),380);

for kk = 1:length(f_list)
    for ll = 1:length(uu_max_list)
        
        
        f_byproduct = f_list(kk);
        uu_max = uu_max_list(ll);
        
        for pa=1:380
            
            %waitbar(((ll)+5*(kk-1))/25,f1);
            
            b_real = abundance_chia_full(:,pa) / sum(abundance_chia_full(:,pa));

            %f = 0.9;  % define "leakage" f
            %f_byproduct = 0.5;  % define "leakage" f
            % i_selfish calculates all microbe ids which
            % do not have any byproduct. You can set their f to 0.
            deg = reshape(sum(full(bp_out),1), [2244,1]);
            i_selfish = i_b_all(find(deg(i_b_all) == 0));
            f = f_byproduct .* ones(2244,1);
            f(i_selfish) = 0.0;
            %uu_max = 10;

            % A_in, A_out matrix [c2b_real, bp_out] and [mu_matrix, ctot_matrix] is calculated
            [c2b_real, bp_out] = Ain_out(b_real, i_all_filt, j_all_filt, v_all_filt);
            [mu_matrix, ctot_matrix, mu_layer] = mu(f, c2b_real, bp_out, uu_max);

            % i_fit_diet containins labels of diets we want to optimize over
            %i_met_patient = find(Thai_metabolome_chia_full_norm(:,pa));
            %i_fit_diet = [i_polysach_diet; i_sugar_diet; 2170];
            %i_fit_diet = aaa;
            %i_fit_diet = i_met_patient;
            %i_fit_diet = [i_polysach_diet; 2170];
            %i_fit_diet = [i_metabolome; 2170];
            i_list  = unique(j_all_filt(find((v_all_filt==7) .* (i_all_filt~=2170)))); 
            i_fit_diet = [i_list([3:(end-2),end]); i_sugar_diet([3,7])];

            % This following part is used to remove all nutrients in the diet
            % not used by the patient:
            x_full_init = zeros(2244,1);
            x_full_init(i_fit_diet) = 1;
            i_diet_used = find((sum(c2b_real,2) .* x_full_init));
            i_fit_diet = i_diet_used;


            lb = zeros(size(i_fit_diet,1),1)'; % lb is lower bounds of variables
            ub = ones(size(i_fit_diet,1),1)' * 100; % ub is upper bounds of variables

            x0 = zeros(length(i_fit_diet),1);  % x0 is initial conditions of variables
            x0(:) = 0.1;
            x0 = x0 /sum(x0);
            A = [];  
            b = [];  % A and b specifies the constrait domain by treating A.x < b
            Aeq = ones(size(i_fit_diet,1),1)';
            beq = [1];  % Aeq and beq specifies the constrait domain by treating Aeq.x = beq
            %nonlcon = [];  % 
            funct = @(x)pred_error_func(x,f, b_real, mu_matrix, c2b_real, bp_out,i_fit_diet);
            %options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','iter','TolFun',1e-4,'TolX',1e-4);
            options = optimoptions(@lsqnonlin,'Algorithm','Trust-region-reflective','Display','iter','TolFun',1e-4,'TolX',1e-4);
            x = lsqnonlin(funct,x0,lb, ub, options);

            x_full = zeros(2244,1);
            x_full(i_fit_diet) = x;
            %predicted_b_a_layer = zeros(2244,2244);
            c_layer = zeros(2244,uu_max);
            c_layer_not_used = zeros(2244,uu_max);
            c_layer_not_used_now = zeros(2244,1);
            for ii = 1:uu_max
                %predicted_b_a_layer = predicted_b_a_layer + mu_layer(:,:,ii);
                %pred_error_layer = 
                c_layer(:,ii) = mu_layer(:,:,ii) * x_full;
                i_layer_not_used = find(sign(c_layer(:,ii)) - sign(sum(c2b_real,2) .* c_layer(:,ii)));
                %disp(i_layer_not_used)
                c_layer_not_used(i_layer_not_used,ii) = c_layer(i_layer_not_used,ii); 
                c_layer_not_used_now = c_layer_not_used_now + c_layer_not_used(:,ii);
                c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used_now;
                %c_layer_not_used(i_layer_not_used,ii) = x_full(i_layer_not_used); 
                %c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used(:,ii);
            end
            
            %%{
            b_pred = mu_matrix * x_full;
            b_pred_all(:,kk,ll,pa) = b_pred;
            i_common = find(sign(b_real) .* sign(b_pred));
            if length(i_common) ~= 0
                [corrS,pvalS] = corr(b_real(i_common), b_pred(i_common),'type','Spearman');
                [corrP,pvalP] = corr(b_real(i_common), b_pred(i_common),'type','Pearson');
            end
            %%}
            
            c_layer_not_used_all(:,kk,ll,pa) = c_layer_not_used(:,uu_max);
            c_layer_all(:,kk,ll,pa) = c_layer(:,uu_max);

            diet_all(:,kk,ll,pa) = x_full;
            predicted_b_a = mu_matrix * x_full;
            pred_abun(:,kk,ll,pa) = predicted_b_a;
            
            %%{
            corrS_all(kk,ll,pa) = corrS;
            corrP_all(kk,ll,pa) = corrP;
            pvalS_all(kk,ll,pa) = pvalS;
            pvalP_all(kk,ll,pa) = pvalP; 
            %%}

        end  
    end
end

%%
pa = 27;
b_real = abundance_chia_full(:,pa) / sum(abundance_chia_full(:,pa));
b_pred = b_pred_all(:,1,1,pa);
figure;
plot(b_pred +1e-6, b_real+1e-6, 'ko')
hold on
plot([1e-6, 1],[1e-6, 1],'k-')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1e-6 1 1e-6 1])
xlabel('Predicted abundance from optimized diet','FontSize',15,'Fontweight','Bold')
ylabel('real abundance','FontSize',15,'Fontweight','Bold')

%% Plot predicted vs real abundance of a single patient
for kk=1:length(patient_bad_ids);
    pa=patient_bad_ids(kk);
    b_real = abundance_chia_full(:,pa) / sum(abundance_chia_full(:,pa));
    b_pred = b_pred_all(:,1,1,pa);
    i_common = find(sign(b_real) .* sign(b_pred));
    [corrP,pvalP] = corr(b_pred(i_common), b_real(i_common),'type','Pearson');

    figure;
    plot(b_pred(i_common), b_real(i_common), 'ko')
    hold on
    plot([1e-6, 1],[1e-6, 1],'k-')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    axis([1e-6 1 1e-6 1])
    xlabel('Predicted abundance from optimized diet','FontSize',15,'Fontweight','Bold')
    ylabel('real abundance','FontSize',15,'Fontweight','Bold')
end;
%% Find bacteria that are constantly overestimated or underestimated
b_pred_all_reshaped = reshape(b_pred_all, 2244, 380);
% names_all(find(sum(b_pred_all_reshaped,2)))
% names_all(find(sum(abundance_chia_full,2)))
bact_error = log10(b_pred_all_reshaped + 1e-6) - log10(abundance_chia_full ./ repmat(sum(abundance_chia_full,1), 2244, 1) + 1e-6);

sign_total = sum(sign(bact_error),2);
figure;
hist(sign_total(find(sign_total~=0)));

error_total = sum(bact_error,2) ./ sum(sign(abundance_chia_full),2);
figure;
hist(error_total(find(error_total~=0)),60);

disp('Sign underestimated:')
names_all(find(sign_total < -100))

disp('Underestimated by 2 order:')
names_all(find(error_total < -2))



%% Find the properties of underestimated bacteria which is prevalent and 2 order under 
% -> The story is that majority cannot eat anything. But most are exact matches 
% in Veronika's matching table

i_bad = find((sign_total < -100) .* (error_total < -2));
i1=find((v_all_filt==2)+(v_all_filt==5));
lambda_in_vec=ones(size(i1));
lambda_in=sparse(i_all_filt(i1),j_all_filt(i1),lambda_in_vec,2244,2244);

num_in_links = sum(lambda_in(:,i_bad),1);

figure;
scatter(error_total(i_bad) ,num_in_links, 'ko');

%% Find the properties of overestimated bacteria which is prevalent and 2 order under
i_bad = find((sign_total > 100) .* (error_total > 1.5));
i1=find((v_all_filt==2)+(v_all_filt==5));
lambda_in_vec=ones(size(i1));
lambda_in=sparse(i_all_filt(i1),j_all_filt(i1),lambda_in_vec,2244,2244);

num_in_links = sum(lambda_in(:,i_bad),1);

figure;
scatter(error_total(i_bad) ,num_in_links, 'ko');

for kk = 1:length(i_bad)
    disp([names_all(i_bad(kk)), 'eats metabolites as follows:'])
    disp(names_all(find(lambda_in(:,i_bad(kk)))))
end




%% Layer assignment
uu_max = 4;
met_layer_assign = zeros(2244,380);
bact_layer_assign = zeros(2244,380);
c_lbyl = zeros(2244,uu_max,380);
c_left = zeros(2244,uu_max,380);
b_lbyl = zeros(2244,uu_max,380);
for pa = 1:380
    disp(pa);
    b_real = abundance_chia_full(:,pa) / sum(abundance_chia_full(:,pa));
    f_byproduct = 0.9;
    f = f_byproduct .* ones(2244,1);
    f(i_selfish) = 0.0;
    f_mul = repmat(f, 1, 2244);
    [c2b_real, bp_out] = Ain_out(b_real, i_all_filt, j_all_filt, v_all_filt);
    [mu_matrix, ctot_matrix, mu_layer] = mu(f, c2b_real, bp_out, uu_max);
    
    %c_layer = c_layer_all(:,5,5,pa);
    %x_full = zeros(2244,1);
    x_full = diet_all(:,1,1,pa);
    %[a,b] = sort(x_full(i_fit_diet),'descend');  
    %x_full(i_fit_diet(b(5:end))) = 0;   % keep only top 4 metabolites
    %predicted_b_a_layer = zeros(2244,2244);
    c_layer = zeros(2244,uu_max);
    b_layer = zeros(2244,uu_max);
    c_layer_cum = zeros(2244,uu_max);
    b_layer_cum = zeros(2244,1);
    c_layer_not_used = zeros(2244,uu_max);
    c_layer_not_used_now = zeros(2244,1);
    %b_layer = zeros(2244,uu_max);
    
    b_layer_cum = mu_matrix * x_full;
    
    for ii = 1:uu_max
        %predicted_b_a_layer = predicted_b_a_layer + mu_layer(:,:,ii);
        %pred_error_layer = 
        c_layer(:,ii) = mu_layer(:,:,ii) * x_full;
        b_layer(:,ii) = (1-f_mul) .* c2b_real' * c_layer(:,ii);
        %b_layer(:,ii) = ;
        if 1==1
            c_layer_cum(:,ii) = c_layer(:,ii);
            i_layer_not_used = find(sign(c_layer(:,ii)) - sign(sum(c2b_real,2) .* c_layer(:,ii)));
            %disp(i_layer_not_used)
            c_layer_not_used(i_layer_not_used,ii) = c_layer(i_layer_not_used,ii); 
            c_layer_not_used_now = c_layer_not_used_now + c_layer_not_used(:,ii);
            c_layer_cum(:,ii) = c_layer_cum(:,ii) + c_layer_not_used_now;
        end
        %c_layer_not_used(i_layer_not_used,ii) = x_full(i_layer_not_used); 
        %c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used(:,ii);
    end
    
    %c_layer = reshape(c_layer, 2244, 5);
    i_nonzero_metabolites = find(sum(c_layer,2));
    %figure;
    %imagesc(c_layer(i_nonzero_metabolites,:))
    %imagesc(log10(c_layer(i_nonzero_metabolites,:)+1e-5));
    %colorbar()
    [max_con, max_pos] = max(c_layer(i_nonzero_metabolites,:)');  % peak
    %concentration assignment
    %[max_con, max_pos] = max(sign(c_layer(i_nonzero_metabolites,:))'); % breadth-first assignment
    % c_layer-weight assignment:
    if 1==0
        c_layer_total = sum(c_layer(i_nonzero_metabolites,:),2);
        numLayer = [1,2,3,4,5];
        numLayer = repmat(numLayer,length(i_nonzero_metabolites),1);
        max_pos = sum(c_layer(i_nonzero_metabolites,:) .* numLayer,2) ./ c_layer_total;
    end
    if 1==0
        for ii = 1:max(max_pos)
            i_location = find(max_pos==ii);
            [aa,bb] = sort(log10(c_layer(i_nonzero_metabolites(i_location),ii)),'descend');
            i_location = i_location(bb);
            disp(['layer',num2str(ii)]);
            disp(names_all(i_nonzero_metabolites(i_location)));
            for jj = 1:length(i_location)
                disp([names_all(i_nonzero_metabolites(i_location(jj))), log10(c_layer(i_nonzero_metabolites(i_location(jj)),ii)), sum(i_common == i_nonzero_metabolites(i_location(jj))) > 0]);
                met_layer_assign(i_nonzero_metabolites(i_location(jj)),pa) = ii;
            end
        end
    end
    met_layer_assign(i_nonzero_metabolites,pa) = max_pos';
    
    i_nonzero_bacts = find(sum(b_layer,2));
    [max_bact_con, max_bact_pos] = max(b_layer(i_nonzero_bacts,:)');
    bact_layer_assign(i_nonzero_bacts,pa) = max_bact_pos';
    
    c_lbyl(:,:,pa) = c_layer;
    b_lbyl(:,:,pa) = b_layer;
    c_left(:,:,pa) = c_layer_not_used;
end

%% Layer by layer activity analysis (metabolic and bacterial)
% metabolites produced at each layer
c1=reshape(sum(c_lbyl,1),4,380);
cdiets = reshape(diet_all(:,1,1,:),2244,380);
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
figure; bar(horzcat(c1pp,c2p,b1p), 'stacked'); legend('byproducts','leftover','bacteria');
figure; bar(horzcat(c1pp./sum(c1pp),c2p./sum(c2p),b1p./sum(b1p)), 'grouped'); legend('byproducts','leftover','bacteria');

%%
i_chinese = [1:157]';
i_american = [158:295]';
i_european = [296:380]';

%Chinese
c1p = sum(c1_norm(:,i_chinese),2)/length(i_chinese);
c2p=sum(c2_norm(:,i_chinese),2)/length(i_chinese);
b1p=sum(b1_norm(:,i_chinese),2)/length(i_chinese);
b1p_c_std=std(b1_norm(:,i_chinese),0,2)/sqrt(length(i_chinese));
%c2p_c_std=std(c2_norm(:,i_chinese),0,2)/sqrt(length(i_chinese));
b1p_c=vertcat([0],b1p);
c1pp=c1p-c2p;
c2p(4)=c2p(4)+c1pp(4);
c1pp(4)=0;
c2p_c = c2p;
%figure; bar(horzcat(c1pp,c2p,b1p_c), 'stacked'); legend('byproducts','leftover','bacteria'); title('Chinese');

%%American
c1p = sum(c1_norm(:,i_american),2)/length(i_american);
c2p=sum(c2_norm(:,i_american),2)/length(i_american);
b1p=sum(b1_norm(:,i_american),2)/length(i_american);
b1p_a_std=std(b1_norm(:,i_american),0,2)/sqrt(length(i_american));
%c2p_a_std=std(c2_norm(:,i_american),0,2)/sqrt(length(i_american));
b1p_a=vertcat([0],b1p);
c1pp=c1p-c2p;
c2p(4)=c2p(4)+c1pp(4);
c1pp(4)=0;
c2p_a = c2p;
%figure; bar(horzcat(c1pp,c2p,b1p_a), 'stacked'); legend('byproducts','leftover','bacteria'); title('American');

%%European
c1p = sum(c1_norm(:,i_european),2)/length(i_european);
c2p=sum(c2_norm(:,i_european),2)/length(i_european);
b1p_e_std=std(b1_norm(:,i_european),0,2)/sqrt(length(i_european));
%c2p_e_std=std(c2_norm(:,i_european),0,2)/sqrt(length(i_european));
b1p=sum(b1_norm(:,i_european),2)/length(i_european);
b1p_e=vertcat([0],b1p);
c1pp=c1p-c2p;
c2p(4)=c2p(4)+c1pp(4);
c1pp(4)=0;
c2p_e = c2p;
%figure; bar(horzcat(c1pp,c2p,b1p_e), 'stacked'); legend('byproducts','leftover','bacteria'); title('European');

%%
figure; %plot(b1p_a(2:end), 'ko-'); hold on;
errorbar([2 3 4], b1p_a(2:end), b1p_a_std, 'ko-'); hold on;
errorbar([2 3 4], b1p_c(2:end), b1p_c_std, 'd-'); 
errorbar([2 3 4], b1p_e(2:end), b1p_e_std, 's-');
legend('American', 'Chinese', 'European');
%plot(b1p_c(2:end), 'd-');
%plot(b1p_e(2:end), 's-');


figure; %plot(b1p_a(2:end), 'ko-'); hold on;
%errorbar([2 3 4], c2p_a(2:end), c2p_a_std, 'ko-'); hold on;
%errorbar([2 3 4], c2p_c(2:end), c2p_c_std, 'd-'); 
%errorbar([2 3 4], c2p_e(2:end), c2p_e_std, 's-');
plot([2 3 4], c2p_a(2:end), 'ko-'); hold on;
plot([2 3 4], c2p_c(2:end), 'd-'); 
plot([2 3 4], c2p_e(2:end), 's-');
legend('American', 'Chinese', 'European');

%%
figure; bar(sum(c2p_a(2:end)), sum(c2p_e(2:end)),sum(c2p_c(2:end))); legend('American', 'European','Chinese');
%% Print metabolome
i_layer = 3;
[a,b] = sort(mean(c_left_norm(:,i_layer,i_european),3),'descend');
names_all(b(1:20))
log10(a(1:20))

[a,b] = sort(mean(c_left_norm(:,i_layer,i_american),3),'descend');
names_all(b(1:20))
log10(a(1:20))


%% Calculating z-scores of the metabolomes
all_metabolome=reshape(c_layer_all, 2244,380);
mean_metabolome = mean(all_metabolome,2);
std_metabolome = std(all_metabolome,0, 2);
mean_matrix = repmat(mean_metabolome,1,380);
std_matrix = repmat(std_metabolome,1,380);

mean_metabolome_american = mean(all_metabolome(:,i_american),2);
mean_metabolome_european = mean(all_metabolome(:,i_european),2);
mean_metabolome_chinese = mean(all_metabolome(:,i_chinese),2);
% mean_matrix_american = repmat(mean_metabolome_american,1,380);
% mean_matrix_european = repmat(mean_metabolome_european,1,380);
% mean_matrix_chinese = repmat(mean_metabolome_chinese,1,380);

std_matrix_american = std(all_metabolome(:,i_american),0,2)./sqrt(length(i_american));
std_matrix_european = std(all_metabolome(:,i_european),0,2)./sqrt(length(i_european));
std_matrix_chinese = std(all_metabolome(:,i_chinese),0,2)./sqrt(length(i_chinese));

zscore_american_european = (mean_metabolome_american - mean_metabolome_european) ./ (sqrt(std_matrix_american.^2 + std_matrix_european.^2));
zscore_american_chinese = (mean_metabolome_american - mean_metabolome_chinese) ./ (sqrt(std_matrix_american.^2 + std_matrix_chinese.^2));
zscore_european_chinese = (mean_metabolome_european - mean_metabolome_chinese) ./ (sqrt(std_matrix_european.^2 + std_matrix_chinese.^2));

%%
zscore_chinese = (all_metabolome(:,i_chinese) - mean_matrix(:,i_chinese)) ./ (std_matrix(:,i_chinese));
zscore_american = (all_metabolome(:,i_american) - mean_matrix(:,i_american)) ./ (std_matrix(:,i_american));
zscore_european = (all_metabolome(:,i_european) - mean_matrix(:,i_european)) ./ (std_matrix(:,i_european));
%zscore_thai = (Thai_metabolome_chia_full - mean_thai)./(std_thai);
zscore_all= (all_metabolome - mean_matrix) ./ (std_matrix);

%%
[pca_coeff, pca_score, pca_latent, pca_tsquared, pca_explained] = pca(log10(all_metabolome + 1e-6)');
figure; plot( pca_score(i_chinese,1),  pca_score(i_chinese,2),'o');hold on;
plot( pca_score(i_american,1),  pca_score(i_american,2),'o');
plot( pca_score(i_european,1),  pca_score(i_european,2),'o');
legend('Chinese','American','European');

%% Find the reason behind two clusters of the PCA 
i_left = find(pca_score(:,1) < -5);
i_right = find(pca_score(:,1) >= -5);




%% Find the correlation between metabolites in metabolome with BMI of patients
i_nonzeros=find(sum(all_metabolome,2));
corr_all = zeros(length(i_nonzeros),1);
p_all = zeros(length(i_nonzeros),1);
i_known_bmi_patients=find(HMP_BMI);
for ii=1:length(i_nonzeros);
    [cor,p]=corr(log10(all_metabolome(i_nonzeros(ii),i_known_bmi_patients)+1e-6)',HMP_BMI(i_known_bmi_patients), 'type','Pearson');
    corr_all(ii) = cor;
    p_all(ii) = p;
    if p<0.05;
        disp([names_all(i_nonzeros(ii)),cor, p]);
    end;
end;

[a,b] = sort(p_all);
for ii = 1:10;
    disp([names_all(i_nonzeros(b(ii))),corr_all(b(ii)), p_all(b(ii))]);
end;

%% linear regression over metabolome
%{
ii = 1;
while p_all(b(ii)) < 0.05;
    if ii == 1
        X = log10(all_metabolome(i_nonzeros(b(ii)),i_known_bmi_patients)+1e-6)';
    else
        X = [X, log10(all_metabolome(i_nonzeros(b(ii)),i_known_bmi_patients)+1e-6)'];
    end
    ii = ii + 1;
end
%}

ii_list = [1,5];
for jj = 1:length(ii_list);
    ii = ii_list(jj);
    if ii == 1
        X = log10(all_metabolome(i_nonzeros(b(ii)),i_known_bmi_patients)+1e-6)';
    else
        X = [X, log10(all_metabolome(i_nonzeros(b(ii)),i_known_bmi_patients)+1e-6)'];
    end
    ii = ii + 1;
end

y = HMP_BMI(i_known_bmi_patients);

mdl = fitlm(X,y)

figure;
scatter(y,mdl.Fitted,'ko');

corr(y,mdl.Fitted, 'type','Pearson')

% aa = all_metabolome(i_nonzeros(b(ii_list)),i_known_bmi_patients);
% patients_new(find(sum(aa, 1)==0))
% names_all(i_nonzeros(b(ii_list)))


%% Find the correlation between bacterial metagemomics with BMI of patients
abundance_chia_full_norm = abundance_chia_full ./ repmat(sum(abundance_chia_full,1),2244,1);
i_nonzeros=find(sum(abundance_chia_full_norm,2));
corr_all = zeros(length(i_nonzeros),1);
p_all = zeros(length(i_nonzeros),1);
i_known_bmi_patients=find(HMP_BMI);
for ii=1:length(i_nonzeros);
    [cor,p]=corr(abundance_chia_full_norm(i_nonzeros(ii),i_known_bmi_patients)',HMP_BMI(i_known_bmi_patients), 'type','Spearman');
    corr_all(ii) = cor;
    p_all(ii) = p;
    if p<0.05;
        %disp([names_all(i_nonzeros(ii)),cor, p]);
    end;
end;

[a,b] = sort(p_all);
for ii = 1:10;
    disp([names_all(i_nonzeros(b(ii))),corr_all(b(ii)), p_all(b(ii))]);
end;

%% linear regression over bacterial metagenomics
ii = 1;
%{
while p_all(b(ii)) < 0.05;
    if ii == 1
        X = log10(abundance_chia_full_norm(i_nonzeros(b(ii)),i_known_bmi_patients)+1e-6)';
    else
        X = [X, log10(abundance_chia_full_norm(i_nonzeros(b(ii)),i_known_bmi_patients)+1e-6)'];
    end
    ii = ii + 1;
end
%}

for ii = 1:19
    if ii == 1
        X = log10(abundance_chia_full_norm(i_nonzeros(b(ii)),i_known_bmi_patients)+1e-6)';
    else
        X = [X, log10(abundance_chia_full_norm(i_nonzeros(b(ii)),i_known_bmi_patients)+1e-6)'];
    end
    ii = ii + 1;
end

y = HMP_BMI(i_known_bmi_patients);

mdl = fitlm(X,y)

figure;
scatter(y,mdl.Fitted,'ko');

corr(y,mdl.Fitted, 'type','Pearson')



%% Sort the p-value of correlation in increasing order and scatter plot the most 
% significant metabolite with patient's BMI.
i_nonzeros = find(sum(all_metabolome, 2));
[a, b] = sort(p_all);
names_all(i_nonzeros(b))
corr_all(b)
%%
figure;
semilogx(all_metabolome(i_nonzeros(b(5)),i_known_bmi_patients),HMP_BMI(i_known_bmi_patients),'ko')


%%
i_nonzeros = find(sum(all_metabolome,2));
%CGobj = clustergram(all_metabolome(i_nonzeros,:)', 'Standardize','Column', 'RowLabels',patients_new','ColumnLabels',names_all(i_nonzeros),'ColumnPDist','correlation', 'RowPDist','correlation');
CGobj = clustergram(log10(all_metabolome(i_nonzeros,:))', 'RowLabels',patients_new','ColumnLabels',names_all(i_nonzeros),'ColumnPDist','correlation', 'RowPDist','correlation');


%%
b_lbyl1 = b_lbyl(:,1:end-1,:);
c_lbyl1 = c_lbyl(:,2:end,:);
c_left1 = c_left(:,2:end,:);

c1=reshape(sum(c_lbyl,1),4,380);
%cdiets = reshape(diet_all(:,1,1,:),2244,380);

norm_factor = repmat(reshape(sum(cdiets,1), 1, 1, 380), [2244,3,1]);

b_norm = b_lbyl1./norm_factor;

c_left_norm = c_left1./norm_factor;
c_lbyl_norm = c_lbyl1./norm_factor;
c_left_norm(:,end,:) = c_lbyl_norm(:,end,:);

%% alpha and gamma diversity for microbes
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
figure; 
bar(2:4,[alpha_D2_mean,beta_D2,gamma_D'./4])
% [AX,H1,H2] =plotyy(2:4,[alpha_D2_mean,beta_D2], 2:4,gamma_D', 'bar', 'bar');
% set(H1,'FaceColor','r') % a
% set(H2,'FaceColor','b') % b
% hLgnd=legend([H1(1) H1(2) H2(1)],'alpha','beta','gamma');
%% alpha and gamma diversity for metabolites using leftovers
for pa=1:380; 
    b2p=c_left_norm(:,:,pa)./repmat(sum(c_left_norm(:,:,pa),1),2244,1,1); 
    alpha_D_met(1:3,pa)=1./sum(b2p.^2,1); 
    alpha_D_shannon_met(1:3,pa)=exp(-sum(b2p.*log(max(b2p,1e-10)),1)); 
end;
alpha_D2_met_mean = 1 ./ mean(1./alpha_D_met,2);

b_norm_total = sum(c_left_norm,3);
b2p = b_norm_total ./ repmat(sum(b_norm_total,1),2244,1);
gamma_D_met(1:3) = 1./sum(b2p.^2,1);
gamma_D_shannon_met(1:3)=exp(-sum(b2p.*log(max(b2p,1e-10)),1)); 
%
beta_D_met = gamma_D_met' ./ mean(alpha_D_met,2);
beta_D2_met = gamma_D_met' .* mean(1./alpha_D_met,2);

beta_D_shannon_met = gamma_D_shannon_met' ./ mean(alpha_D_shannon_met,2);
beta_D2_shannon_met = gamma_D_shannon_met' .* mean(1./alpha_D_shannon_met,2);

for pa=1:380; 
    b2p=sum(c_left_norm(:,:,pa),2) ./ sum(sum(c_left_norm(:,:,pa),2));
    alpha_D_1l_met(pa)=1./sum(b2p.^2,1); 
end;
%% alpha, beta and gamma diversity for diets
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

%%
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
%bar(horzcat(alpha_all, beta_all, gamma_all, beta_bact_all)')

%%
[pca_coeff, pca_score, pca_latent, pca_tsquared, pca_explained] = pca(Deff_p');
figure; plot( pca_score(:,1),  pca_score(:,2),'ko');

%%
pa = 24;
c1p = c1_norm(:,pa);
c2p = c2_norm(:,pa);
b1p =b1_norm(:,pa);
b1p=vertcat([0],b1p);
c1pp=c1p-c2p;
figure; bar(horzcat(c1pp,c2p,b1p), 'stacked'); legend('byproducts','leftover','bacteria');


%%
mean_met_layer_assign = sum(met_layer_assign,2) ./ max(sum(met_layer_assign>0,2),1);

i_nonzero = find(mean_met_layer_assign);
std_met_layer_assign = zeros(2244,1);

for kk = 1:length(i_nonzero)
    i_nonzero_layers = find(met_layer_assign(i_nonzero(kk),:));
    std_met_layer_assign(i_nonzero(kk)) = std(met_layer_assign(i_nonzero(kk),i_nonzero_layers));
end

%std_met_layer_assign = std(met_layer_assign,0,2);

i_nonzero = find(mean_met_layer_assign);

[x, y] = sort(mean_met_layer_assign(i_nonzero),'ascend');

for kk = 1:length(i_nonzero)
    disp([names_all(i_nonzero(y(kk))), mean_met_layer_assign(i_nonzero(y(kk))), std_met_layer_assign(i_nonzero(y(kk)))]);
    %disp([names_all(i_nonzero(y(kk))), mean_met_layer_assign(i_nonzero(y(kk)))]);
end

figure;
scatter(mean_met_layer_assign(i_nonzero(y)),std_met_layer_assign(i_nonzero(y)),'ko')


%%
pa = 24; %24;
metabolome = Thai_metabolome_chia_full_norm(:,pa);
b_real = Thai_abundance_chia_full_norm(:,pa);
f_byproduct = 0.9;
f = f_byproduct .* ones(2244,1);
f_mul = repmat(f, 1, 2244);
[c2b_real, bp_out] = Ain_out(b_real, i_all_filt, j_all_filt, v_all_filt);
[mu_matrix, ctot_matrix, mu_layer] = mu(f, c2b_real, bp_out, uu_max);

c_layer = zeros(2244,uu_max);
b_layer = zeros(2244,uu_max);
c_layer_cum = zeros(2244,uu_max);
c_layer_not_used = zeros(2244,uu_max);
c_layer_not_used_now = zeros(2244,1);

x_full = diet_all(:,1,1,pa);
for ii = 1:uu_max
    %predicted_b_a_layer = predicted_b_a_layer + mu_layer(:,:,ii);
    %pred_error_layer = 
    c_layer(:,ii) = mu_layer(:,:,ii) * x_full;
    b_layer(:,ii) = (1-f_mul) .* c2b_real' * c_layer(:,ii); 
    i_layer_not_used = find(sign(c_layer(:,ii)) - sign(sum(c2b_real,2) .* c_layer(:,ii)));
    disp(i_layer_not_used)
    c_layer_not_used(i_layer_not_used,ii) = c_layer(i_layer_not_used,ii); 
    c_layer_not_used_now = c_layer_not_used_now + c_layer_not_used(:,ii);
    c_layer_cum(:,ii) = c_layer(:,ii);
    c_layer_cum(:,ii) = c_layer_cum(:,ii) + c_layer_not_used_now;
    %c_layer_not_used(i_layer_not_used,ii) = x_full(i_layer_not_used); 
    %c_layer(:,ii) = c_layer(:,ii) + c_layer_not_used(:,ii);
end

ii=5;
i_common = find(sign(c_layer_cum(:,ii)) .* sign(metabolome));
a = met_layer_assign(i_common,pa);
for kk = 1:length(a)
    legend_names{kk} = num2str(a(kk));
end

b_layer_cum = mu_matrix * x_full;
i_bact_common = find(b_layer_cum);
figure;
%p1=loglog(c_layer_cum(i_common,ii) / sum(c_layer_cum(i_common,ii)),metabolome(i_common)/sum(metabolome(i_common)),'o')
scatter(b_layer_cum(i_bact_common), b_real(i_bact_common), 100, bact_layer_assign(i_bact_common,pa), 'filled');
%legend(names_all(i_common));
hold on;
p2=loglog([1e-5, 1],[1e-5, 1],'k-')
axis([1e-5, 1, 1e-5, 1])
set(gca,'XScale','log')
set(gca,'YScale','log')
colorbar()

figure;
%p1=loglog(c_layer_cum(i_common,ii) / sum(c_layer_cum(i_common,ii)),metabolome(i_common)/sum(metabolome(i_common)),'o')
scatter(c_layer_cum(i_common,ii) / sum(c_layer_cum(i_common,ii)),metabolome(i_common)/sum(metabolome(i_common)), 100, a(:), 'filled');
%legend(names_all(i_common));
hold on;
p2=loglog([1e-5, 1],[1e-5, 1],'k-')
axis([1e-5, 1, 1e-5, 1])
set(gca,'XScale','log')
set(gca,'YScale','log')
colorbar()

%%
pred_mean_metabolome = sum(c_layer_all(:,1,1,:),4) ./ max(sum(c_layer_all(:,1,1,:)>0,4),1);

Thai_mean_metabolome = sum(Thai_metabolome_chia_full_norm,2) ./ max(sum(Thai_metabolome_chia_full_norm>0,2),1);

i_common = find(sign(pred_mean_metabolome) .* sign(obs_mean_metabolome));

%{
ii=5;
i_common = find(sign(c_layer_cum(:,ii)) .* sign(metabolome));
a = met_layer_assign(i_common,pa);
for kk = 1:length(a)
    legend_names{kk} = num2str(a(kk));
end
%}


figure;
%p1=loglog(c_layer_cum(i_common,ii) / sum(c_layer_cum(i_common,ii)),metabolome(i_common)/sum(metabolome(i_common)),'o')
scatter(pred_mean_metabolome(i_common) / sum(pred_mean_metabolome(i_common)),obs_mean_metabolome(i_common)/sum(obs_mean_metabolome(i_common)), 100, mean_met_layer_assign(i_common), 'filled');
%legend(names_all(i_common));
hold on;
p2=loglog([1e-5, 1],[1e-5, 1],'k-')
axis([1e-5, 1, 1e-5, 1])
set(gca,'XScale','log')
set(gca,'YScale','log')
colorbar()



%legend([p1,p2],cellstr(legend_names),'line');


%% histogram
figure;
hist(mean_met_layer_assign(i_nonzero(y)))



%% Cluster the diet from the best f and best layer
diet_to_cluster = zeros(2244,380);
for pa=1:380
    diet_to_cluster(:,pa) = diet_all(:,1,1,pa);
end
%i_fit_diet = [i_polysach_diet; i_sugar_diet; 2170];
i_list  = unique(j_all_filt(find((v_all_filt==7) .* (i_all_filt~=2170)))); 
i_fit_diet = [i_list([3:(end-2),end]); i_sugar_diet([3,7])];
cd=diet_to_cluster(i_fit_diet,:);
%i_fit_diet_filt = i_fit_diet([1:2, 4:end]);
i_fit_diet_filt = i_fit_diet(find(std(cd_filt,0,2) ~= 0));
cd_filt = diet_to_cluster(i_fit_diet_filt,:);
%CGobj = clustergram(cd, 'Standardize','Column', 'RowLabels',names_all(i_fit_diet),'ColumnLabels',Thai_names_all');
CGobj = clustergram(cd_filt, 'Standardize','Row', 'RowLabels',names_all(i_fit_diet_filt),'ColumnLabels',patients_new','ColumnPDist','correlation', 'RowPDist','correlation');
%%
% mean k_in and k_out in layers
k_in=sum(sign(lambda_in), 1);
k_out=sum(sign(bp_out), 1);
sum(k_out>0)
max(k_out)
max(k_in)
mean(k_out(i_b_all))
mean(k_in(i_b_all))
b1=mean(b_norm, 3);
for m=1:3; b1(:,m)=b1(:,m)./sum(b1(:,m)); end;
for m=1:3; k_out_mean(m)=sum(k_out'.*b1(:,m)); k_in_mean(m)=sum(k_in'.*b1(:,m)); end;
k_out_mean
k_in_mean