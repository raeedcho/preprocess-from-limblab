function h = arrayWaveLayout(cds)

% Loop along units and plot mean of all waveforms in the appropriate location
unit_colors = {'k','b','r','g','c','m','y'};

arrayList = unique({cds.units.array});
h = cell(size(arrayList));

for array_idx = 1:numel(arrayList)
    h{array_idx} = figure;
    subplot(10,10,1);

    idx = find(strcmpi({cds.units.array},arrayList{array_idx}));

    for i = idx
        if cds.units(i).ID > 0 && cds.units(i).ID ~= 255
            row = cds.units(i).rowNum;
            col = cds.units(i).colNum;

            subplot(10,10,10*(row-1)+col); hold all;
            plot(mean(cds.units(i).spikes.wave,1),'-','LineWidth',1,'Color',unit_colors{cds.units(i).ID});
        end
    end
    for i = 1:10
        for j = 1:10
            subplot(10,10,10*(i-1)+j);
            set(gca,'Box','off','TickDir','out','XTick',[],'YTick',[]);
            ax = gca; 
            ax.XAxis.Visible = 'off';
            ax.YAxis.Visible = 'off';
            axis('square');
        end
    end
end

disp('Done.');