% This file "Data_processing.m" involves the processing of the data to the
% desired form wanted by the "main_code.m"


%% The number of patients in metabolome data is different 
% from the number of patients in abundance data, a common
% intersection of patients is needed to be found
[metabolome_Thai_num,metabolome_Thai_txt,metabolome_Thai_raw]  = xlsread('metabolome_matched_thai_modified_by_Tong.xlsx');
abundance_Thai  = readtable('abundance_matched_thai.txt');
% abundance_Thai.Properties.VariableNames(9:end)'
% metabolome_Thai_txt(1,6:end)'
%
length(abundance_Thai.Properties.VariableNames(9:end)')
length(metabolome_Thai_txt(1,6:end)')

for ii = 1:length(abundance_Thai.Properties.VariableNames(9:end)')
    if ii <= 42
        disp([abundance_Thai.Properties.VariableNames(8+ii) metabolome_Thai_txt(1,5+ii)])
    else
        disp([abundance_Thai.Properties.VariableNames(8+ii) ''])
    end
end

disp('Intersection as follows:')

valid_index1 = [1:12, 14:28, 30:39, 41, 43:45]; % It is manually found
valid_index2 = [1:15, 17:42]; % It is manually found

for ii = 1:length(valid_index1)
    disp([abundance_Thai.Properties.VariableNames(8+valid_index1(ii)) metabolome_Thai_txt(1,5+valid_index2(ii))])
end

for ii = 1:length(valid_index1)
    disp(strcmp(metabolome_Thai_txt(1,5+valid_index2(ii)),abundance_Thai.Properties.VariableNames(8+valid_index1(ii))))
end

Thai_patients_names_all = abundance_Thai.Properties.VariableNames(8+valid_index1);

%% Construct the data in the format we are familiar with
%abundance_Thai=table2array(abundance_Thai(:,8+valid_index1));
Thai_abundance_chia_full = zeros(2244,41);
for ii = 1:size(abundance_Thai,1)
    if abundance_Thai(ii,5).Chia_id ~= 0
        Thai_abundance_chia_full(abundance_Thai(ii,5).Chia_id,:) = Thai_abundance_chia_full(abundance_Thai(ii,5).Chia_id,:) + table2array(abundance_Thai(ii,8+valid_index1));
    end
end
Thai_abundance_chia_full_norm = zeros(2244,41);
for ii = 1:size(Thai_abundance_chia_full,2)
    Thai_abundance_chia_full_norm(:,ii) = Thai_abundance_chia_full(:,ii) / sum(Thai_abundance_chia_full(:,ii));
end



Thai_metabolome_chia_full = zeros(2244,41);
for ii = 1:size(metabolome_Thai_num,1)
    if metabolome_Thai_num(ii,1) ~= 0
        Thai_metabolome_chia_full(metabolome_Thai_num(ii,1),:) = Thai_metabolome_chia_full(metabolome_Thai_num(ii,1),:) + metabolome_Thai_num(ii,2+valid_index2);
    end
end
Thai_metabolome_chia_full_norm = zeros(2244,41);
for ii = 1:size(Thai_metabolome_chia_full,2)
    Thai_metabolome_chia_full_norm(:,ii) = Thai_metabolome_chia_full(:,ii) / sum(Thai_metabolome_chia_full(:,ii));
end
