all_fnames = SPM.xY.P

this_group = {};
this_indiv = {};
for i = 1:length(all_fnames)
this_fname = strsplit(all_fnames{i},'/');
this_group{end+1} = this_fname{end-3};
this_indiv{end+1} = this_fname{end-2};
end

try
all_real_diags = {};
for j = 1:length(Participant)
all_real_diags{end+1} = Participant{j}.diag;
end
all_real_groups = unique(all_real_diags,'stable');
end
all_groups = unique(this_group,'stable');

for j = 1:length(all_groups)
disp(['There are ' num2str(nnz(strcmp(this_group,all_groups{j}))/13) ' ' all_groups{j} ' in the design matrix'])
try
disp(['There are ' num2str(nnz(strcmp(all_real_diags,all_real_groups{j}))) ' ' all_real_groups{j} ' in reality'])
end
end

all_indivs = unique(this_indiv,'stable');
for j = 1:length(all_indivs)
disp(['There are ' num2str(nnz(strcmp(this_indiv,all_indivs{j}))) ' files for ' all_indivs{j} ' in the design matrix'])
end

%all_fnames{strcmp(this_indiv,all_indivs{j})}