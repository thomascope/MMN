all_MMN_0_400 = all_MMN(:,:,26:126);

Ns = size(squeeze(all_MMN_0_400(i,group_inds==grp,:)),2);
X0 = spm_dctmtx(Ns,1);
R = speye(Ns) - X0*X0';
%R_Tuk = R*sparse(diag((tukeywin(Ns,.3)+.025)))*R;
R_Han = R*sparse(diag(hanning(Ns)))*R;
R_mean = R*sparse(diag(Ns))*R;

all_demeaned_MMN = [];
all_Tukey_MMN = [];
all_Hanning_MMN = [];
for i = 1:size(all_MMN_0_400,1)
for j = 1:size(all_MMN_0_400,2)
all_demeaned_MMN(i,j,:) = R_mean*squeeze(all_MMN_0_400(i,j,:));
%all_Tukey_MMN(i,j,:) = R_Tuk*squeeze(all_MMN_0_400(i,j,:));
all_Hanning_MMN(i,j,:) = R_Han*squeeze(all_MMN_0_400(i,j,:));
end
end

figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(all_demeaned_MMN(i,group_inds==grp,:)),0.2,cmap(grp,:),D{1}.time(26:126),1,1);
end
title('Demeaned','FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time(26:126),zeros(1,length(D{1}.time(26:126))),'k--','LineWidth',2)
ylabel('Mismatch Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
legend(linehandle,groups)
end
end

figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(all_Hanning_MMN(i,group_inds==grp,:)),0.2,cmap(grp,:),D{1}.time(26:126),1,1);
end
title('Hanning','FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time(26:126),zeros(1,length(D{1}.time(26:126))),'k--','LineWidth',2)
ylabel('Mismatch Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
legend(linehandle,groups)
end
end

% figure;
% for i = 1:8
% subplot(2,4,i)
% hold on
% for grp = 1:length(groups)
% linehandle(grp) = stdshade_TEC(squeeze(all_Tukey_MMN(i,group_inds==grp,:)),0.2,cmap(grp,:),D{1}.time(26:126),1,1);
% end
% title('Tukey','FontSize',34)
% xlabel('Time (s)')
% xlim([-0.1 0.500])
% plot(D{1}.time(26:126),zeros(1,length(D{1}.time(26:126))),'k--','LineWidth',2)
% ylabel('Mismatch Response (AU)')
% set(gca,'FontWeight','bold')
% set(gca,'LineWidth',2)
% if i==4
% legend(linehandle,groups)
% end
% end

figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(all_MMN_0_400(i,group_inds==grp,:)),0.2,cmap(grp,:),D{1}.time(26:126),1,1);
end
title('Raw','FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time(26:126),zeros(1,length(D{1}.time(26:126))),'k--','LineWidth',2)
ylabel('Mismatch Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
    legend(linehandle,groups)
end
end