function DoDCMTA1_ecd2(n)


F = {'/imaging/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/forwardmodels/ftd/meg11_0179/fmraedfffmeg11_0179.mat'};

% open pool:
% P = parpool(cbupool(length(F))); % P = matlabpool(cbupool(length(F))); %P = parpool(2);
% assignin('base','P',P)

for KK = n %1:21
    %%
    for k = 1%:length(F)
        %try
            dodcmTA1_ecd(F{k},KK)
        %catch err
        %end
    end

end

% close pool:
%parpool close force CBU_Cluster % matlabpool close force CBU_Cluster
%delete(P)


end

