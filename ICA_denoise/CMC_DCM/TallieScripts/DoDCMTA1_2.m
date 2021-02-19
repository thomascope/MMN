function DoDCMTA1(n)


F = FilesDCMTA;
F = F{1};

% open pool:
P = parpool(cbupool(length(F))); % P = matlabpool(cbupool(length(F))); %P = parpool(2);
assignin('base','P',P)

for KK = n %1:21
    %%
    parfor k = 1:length(F)
        try
            dodcmTA1(F{k},KK)
        catch err
        end
    end

end

% close pool:
%parpool close force CBU_Cluster % matlabpool close force CBU_Cluster
delete(P)


end

