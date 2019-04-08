%% calculating the HMP data from Veronika's new mappings
% read Veronika's new mapping of microbes in 380 patients to Chia namespace
clear all;

[a,b]=xlsread('./abundance_table_matched.xlsx');
abundance_table_full=a(:,3:end);
i1=find(a(:,1)); whos i1;
bugs_2_microbes_full=a(:,1);
for m=1:380; 
    a1=sparse(bugs_2_microbes_full(i1),ones(size(i1)), abundance_table_full(i1,m), 2244,1); 
    abundance_chia_full(1:2244,m)=a1; 
end;

%% Running the phase diagram
f_list = [0.9];   % a list of byproduct fraction f calculated
numLayer_max_list = [4]; % a list of number of layers N_l calculated

load('chia_network_new.mat')

intake_all = zeros(2244,length(f_list),length(numLayer_max_list),380);
b_pred_all = zeros(2244,length(f_list),length(numLayer_max_list),380);
corrS_all = zeros(length(f_list),length(numLayer_max_list),380);
corrP_all = zeros(length(f_list),length(numLayer_max_list),380);
pvalS_all = zeros(length(f_list),length(numLayer_max_list),380);
pvalP_all = zeros(length(f_list),length(numLayer_max_list),380);
c_layer_not_used_all = zeros(2244,length(f_list),length(numLayer_max_list),380);
c_layer_all = zeros(2244,length(f_list),length(numLayer_max_list),380);
c_layer_cum_all = zeros(2244,length(f_list),length(numLayer_max_list),380);
c_lbyl = zeros(2244,numLayer_max_list(1),380);
c_left = zeros(2244,numLayer_max_list(1),380);
b_lbyl = zeros(2244,numLayer_max_list(1),380);

for kk = 1:length(f_list)
    for ll = 1:length(numLayer_max_list)
        
        f_byproduct = f_list(kk);
        numLayer_max = numLayer_max_list(ll);
        
        for pa=1:380
            disp(pa)
            b_real = abundance_chia_full(:,pa) / sum(abundance_chia_full(:,pa));

            % i_selfish calculates all microbe ids which
            % do not produce any byproduct. You can set their f to 0.
            f = f_byproduct .* ones(2244,1);
            f(i_selfish) = 0.0;
            f_mul = repmat(f, 1, 2244);

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
            %options = optimoptions(@lsqnonlin,'Algorithm','Trust-region-reflective','Display','iter','TolFun',1e-4,'TolX',1e-4);
            options = optimoptions(@lsqnonlin,'Algorithm','Trust-region-reflective','TolFun',1e-4,'TolX',1e-4);
            x = lsqnonlin(funct,x0,lb, ub, options);

            % predicted metabolome vs real metabolome
            x_full = zeros(2244,1);
            x_full(i_intake) = x;
            b_layer = zeros(2244,numLayer_max);
            c_layer = zeros(2244,numLayer_max);
            c_layer_cum = zeros(2244,numLayer_max);
            c_layer_not_used = zeros(2244,numLayer_max);
            c_layer_not_used_now = zeros(2244,1);

            %b_layer_cum = m2b_total * x_full;

            for ii = 1:numLayer_max
                % m2m_layer transform the nutrient intake into nutrient at each layer
                c_layer(:,ii) = m2m_layer(:,:,ii) * x_full;
                c_layer_cum(:,ii) = c_layer(:,ii);
                b_layer(:,ii) = (1-f_mul) .* m2b' * c_layer(:,ii);
                i_layer_not_used = find(sign(c_layer(:,ii)) - sign(sum(m2b,2) .* c_layer(:,ii)));
                c_layer_not_used(i_layer_not_used,ii) = c_layer(i_layer_not_used,ii); 
                c_layer_not_used_now = c_layer_not_used_now + c_layer_not_used(:,ii);
                c_layer_cum(:,ii) = c_layer_cum(:,ii) + c_layer_not_used_now;
            end
    
            c_lbyl(:,:,pa) = c_layer;
            b_lbyl(:,:,pa) = b_layer;
            c_left(:,:,pa) = c_layer_not_used;


            b_pred = m2b_total * x_full;
            b_pred_all(:,kk,ll,pa) = b_pred;
            i_common = find(sign(b_real) .* sign(b_pred));
            if length(i_common) ~= 0
                [corrS,pvalS] = corr(b_real(i_common), b_pred(i_common),'type','Spearman');
                [corrP,pvalP] = corr(b_real(i_common), b_pred(i_common),'type','Pearson');
            end
            
            c_layer_not_used_all(:,kk,ll,pa) = c_layer_not_used(:,numLayer_max);
            c_layer_all(:,kk,ll,pa) = c_layer(:,numLayer_max);
            c_layer_cum_all(:,kk,ll,pa) = c_layer_cum(:,numLayer_max);

            intake_all(:,kk,ll,pa) = x_full;
            %predicted_b_a = m2b_total * x_full;
            %pred_abun(:,kk,ll,pa) = predicted_b_a;
            
            corrS_all(kk,ll,pa) = corrS;
            corrP_all(kk,ll,pa) = corrP;
            pvalS_all(kk,ll,pa) = pvalS;
            pvalP_all(kk,ll,pa) = pvalP; 

        end  
    end
end


%save('HMP_prediction.mat','*_all','c_lbyl','b_lbyl','c_left')