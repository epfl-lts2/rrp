function [ ] = set_baw_color(  )
%SET_BAW_COLOR Set all color to black and white
    
% for imagesc
colormap(flipud(gray))
% for plot
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',1.5)
    % add other symbol if more than one plot
    for ii=1:length(hline)
        switch mod(ii,4)
            case 1
                set(hline(ii),'LineStyle','-') 
                set(hline(ii),'color','k');
            case 2
                set(hline(ii),'LineStyle',':') 
                set(hline(ii),'color','k');
            case 3
                set(hline(ii),'LineStyle','-') 
                set(hline(ii),'color',[0.5 0.5 0.5]);

            case 0
                set(hline(ii),'LineStyle','-') 
                set(hline(ii),'color',[0.75 0.75 0.75]);
        end
    end

end

