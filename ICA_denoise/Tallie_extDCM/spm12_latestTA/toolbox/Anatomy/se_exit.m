global defaults

fg = se_figure('getWin','Graphics');
delete(fg)
clear st MAP SPM xSPM CLUSTER displayType group defaults;


try; delete(st.FX); end
try
    spm_defaults = defaults.oldDefaults;
catch
    spm_defaults;
end
spm_figure('getWin','Graphics');