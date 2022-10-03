function visualize(tbl,target_site,membrane_width,npc_centroids,npc_site)

npc_list = unique(tbl{:,'npc_index'});

col = lines(max(npc_list));

figure('position',[0,0,600,300]);
clf;

subplot(1,2,1);
hold on;
axis equal;

if any(tbl{:,'npc_index'} == 0)
    plot(tbl{tbl{:,'npc_index'} == 0, 'x [nm]'}, ...
         tbl{tbl{:,'npc_index'} == 0, 'y [nm]'}, ...
         'x', 'color', ones(1,3)*0.7);
    npc_list_start = 2;
else
    npc_list_start = 1;
end

for i_npc = npc_list(npc_list_start:end)'
    plot(tbl{tbl{:,'npc_index'} == i_npc, 'x [nm]'}, ...
         tbl{tbl{:,'npc_index'} == i_npc, 'y [nm]'}, ...
         'x', 'color', col(i_npc,:));
end

xlim([0,membrane_width]);
ylim([0,membrane_width]);
xlabel('x (nm)');
ylabel('y (nm)');

box on;

if nargin > 3
    plot(squeeze(npc_centroids(1,1,:)),squeeze(npc_centroids(2,1,:)),'ko');
end

if nargin > 4
    plot(reshape(npc_site(1,:,:),1,[]),reshape(npc_site(2,:,:),1,[]),'kx');
end

subplot(1,2,2);
if exist('ash2','file') == 2
    img = ash2(tbl{:,'x [nm]'},tbl{:,'y [nm]'},10,[0,membrane_width],[0,membrane_width])';
else
    img = histcounts2(tbl{:,'x [nm]'},tbl{:,'y [nm]'},[0,membrane_width],[0,membrane_width],'binwidth',10)';
end

imshow(img,[]);
set(gca,'YDir','normal');

sgtitle(upper(target_site));

formatfig(gcf,'presentation');

end