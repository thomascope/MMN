%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4899 $)
%-----------------------------------------------------------------------
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.type = 'cfg_files';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.name = '.mat Filename';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.tag = 'matname';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.filter = 'mat';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.ufilter = '.*';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.dir = '';
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.num = [1 1];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.check = [];
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.help = {'The name of the .mat file to load.'};
matlabbatch{1}.menu_cfg{1}.menu_entry{1}.conf_files.def = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_const.type = 'cfg_const';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_const.name = 'All Variables';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_const.tag = 'allvars';
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_const.val = {true};
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_const.check = [];
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_const.help = {'Load all variables found in the .mat file.'};
matlabbatch{2}.menu_cfg{1}.menu_entry{1}.conf_const.def = [];
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.type = 'cfg_entry';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.name = 'Variable Name';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.tag = 'varname';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.strtype = 's';
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.extras = [];
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.num = [1 Inf];
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.check = @(job)cfg_load_vars('check','isvarname',job);
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.help = {'Enter the name of the variable to be loaded from the .mat file.'};
matlabbatch{3}.menu_cfg{1}.menu_entry{1}.conf_entry.def = [];
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.type = 'cfg_repeat';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.name = 'Specified Variables';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.tag = 'varnames';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1) = cfg_dep;
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tname = 'Values Item';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).tgt_spec = {};
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).sname = 'Entry: Variable Name (cfg_entry)';
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.values{1}(1).src_output = substruct('()',{1});
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.num = [1 Inf];
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.forcestruct = false;
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.check = [];
matlabbatch{4}.menu_cfg{1}.menu_struct{1}.conf_repeat.help = {'Enter a list of variable names to be loaded from the .mat file.'};
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.type = 'cfg_choice';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.name = 'Variables to load';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.tag = 'loadvars';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{1}(1) = cfg_dep;
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{1}(1).tname = 'Values Item';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{1}(1).tgt_spec = {};
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{1}(1).sname = 'Const: All Variables (cfg_const)';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{1}(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{1}(1).src_output = substruct('()',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{2}(1) = cfg_dep;
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{2}(1).tname = 'Values Item';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{2}(1).tgt_spec = {};
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{2}(1).sname = 'Repeat: Specified Variables (cfg_repeat)';
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{2}(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.values{2}(1).src_output = substruct('()',{1});
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.check = [];
matlabbatch{5}.menu_cfg{1}.menu_struct{1}.conf_choice.help = {'Choose whether all variables or a list of variables with specified names should be loaded.'};
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.type = 'cfg_exbranch';
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.name = 'Load Variables from .mat File';
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.tag = 'load_vars';
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1) = cfg_dep;
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tname = 'Val Item';
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).tgt_spec = {};
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).sname = 'Files: .mat Filename (cfg_files)';
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{1}(1).src_output = substruct('()',{1});
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1) = cfg_dep;
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tname = 'Val Item';
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).tgt_spec = {};
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).sname = 'Choice: Variables to load (cfg_choice)';
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.val{2}(1).src_output = substruct('()',{1});
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.prog = @(job)cfg_load_vars('run',job);
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.vout = @(job)cfg_load_vars('vout',job);
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.check = [];
matlabbatch{6}.menu_cfg{1}.menu_struct{1}.conf_exbranch.help = {'This function loads variables from a .mat file and passes them on as a dependency. It can load either all variables from a .mat file or a list of specified variables. In the first case, it will return a single struct variable. The variable names in the .mat file will become field names in this struct. If a list of variable names is given, each of them will be loaded into a separate variable.'};
