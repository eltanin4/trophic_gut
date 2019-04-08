clear all;
f_list = [0.1,0.3,0.5,0.7,0.9];   % a list of byproduct fraction f calculated
numLayer_max_list = [2,3,4,5,10]; % a list of number of layers N_l calculated

load('Thai.mat')
load('chia_network_new.mat')

intake_all = zeros(2244,length(f_list),length(numLayer_max_list),41);
pred_abun_all = zeros(2244,length(f_list),length(numLayer_max_list),41);
corrS_all = zeros(length(f_list),length(numLayer_max_list),41);
corrP_all = zeros(length(f_list),length(numLayer_max_list),41);
pvalS_all = zeros(length(f_list),length(numLayer_max_list),41);
pvalP_all = zeros(length(f_list),length(numLayer_max_list),41);
c_layer_not_used_all = zeros(2244,length(f_list),length(numLayer_max_list),41);
c_layer_all = zeros(2244,length(f_list),length(numLayer_max_list),41);

for kk = 1:length(f_list)
    for ll = 1:length(numLayer_max_list)
        
        f_byproduct = f_list(kk); % define "leakage fraction" f
        numLayer_max = numLayer_max_list(ll); % define maximum number of trophic layers
        
        for pa=1:41
            % The simulation is the same as the optimization over an individual 
            disp([kk,ll,pa])

            b_real = Thai_abundance_chia_full_norm(:,pa);  % load the real metagenomic 
                                                           % data of the child pa.

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

            metabolome = Thai_metabolome_chia_full_norm(:,pa);
            i_common = find(sign(c_layer(:,numLayer_max)) .* sign(metabolome));
            if length(i_common) ~= 0
                [corrS,pvalS] = corr(c_layer(i_common,numLayer_max), metabolome(i_common),'type','Spearman');
                [corrP,pvalP] = corr(c_layer(i_common,numLayer_max), metabolome(i_common),'type','Pearson');
            end
            
            c_layer_not_used_all(:,kk,ll,pa) = c_layer_not_used(:,numLayer_max);
            c_layer_all(:,kk,ll,pa) = c_layer(:,numLayer_max);

            intake_all(:,kk,ll,pa) = x_full;
            predicted_b_a = m2b_total * x_full;
            pred_abun_all(:,kk,ll,pa) = predicted_b_a;
            
            corrS_all(kk,ll,pa) = corrS;
            corrP_all(kk,ll,pa) = corrP;
            pvalS_all(kk,ll,pa) = pvalS;
            pvalP_all(kk,ll,pa) = pvalP; 
            
        end  
    end
end

save('Calibration_of_the_model.mat','*all','numLayer_max_list','f_list')