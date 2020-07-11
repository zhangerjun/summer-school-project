function []=plot_overlayed_images(image1,image2,mask,options)

if isfield(options,'outline')
    if options.outline
        outline_mask=options.outline_mask;
    end
else
    options.outline=0;
end

if isfield(options,'rotate')
    angle=options.rotate;
    image1=imrotate(image1,angle);
    image2=imrotate(image2,angle);
    mask=imrotate(mask,angle);
    outline_mask=imrotate(outline_mask,angle);
end

if isfield(options,'fliplr')
    if options.fliplr
        image1=fliplr(image1);
        image2=fliplr(image2);
        mask=fliplr(mask);
        outline_mask=fliplr(outline_mask);
    end
end


if isfield(options,'subplot')
    if options.subplot
        subplot_tight(options.subplot_size(1),...
            options.subplot_size(2),...
            options.subplot_index,...
            options.subplot_margins);
    else
        figure;
    end
else
    figure;
    options.subplot=0;
end

% Create two axes
if options.subplot
    ax1=gca;
else
    ax1 = axes;
end

%imagesc(ax1,imadjust(image1,[0 1],[0 1]));
imagesc(ax1,image1);
if isfield(options,'caxis1')
    if ~isempty(options.caxis1)        
        caxis(options.caxis1)
    end
end

if options.subplot
    %create second axes with same position as first
    ax2 = axes;
    ax2.Position = ax1.Position;
else
    ax2 = axes;
end



h=imagesc(ax2,image2);
set(h,'AlphaData',mask);
set(h,'AlphaDataMapping','scaled');

if options.outline  
    hold on;
    visboundaries(outline_mask,'color','red','linewidth',1)
end

%imagesc(ax2,dwi_structure.(scanid{i}).mask.img(:,:,10))

% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%remove labels from bottom axis
ax1.XTickLabels=[];
ax1.YTickLabels=[];
ax1.XTick = [];
ax1.YTick = [];


% Give each one its own colormap
colormap(ax1,'gray');
if isfield(options,'colormap2')
    if ~isempty(options.colormap2)
        colormap(ax2,options.colormap2);
    end
else
    colormap(ax2,'jet');
end

if isfield(options,'colorbar') && options.colorbar
    %add colorbar for FA
    %c=colorbar(ax2,'Position',[.88 .11 .0675 .815]);        
    %ylabel(c,options.colorbar_label)
    
    c=colorbar(ax2,'location','north','color',[1-eps 1 1]); 
    %shrink_colorbar(c)
    %if options.colorbar_label_on
    %    ylabel(c,options.colorbar_label,'color',[1-eps 1 1])
    %end
    if isfield(options,'colorbar_tick_on')
        if options.colorbar_tick_on
            set(c,'Ticks',[options.cmin options.cmax])
            set(c,'YTickLabel',options.colorbar_tick)
        end
    end
end

if isfield(options,'cmin') && isfield(options,'cmax')   
    caxis([options.cmin options.cmax])  
end
    
if isfield(options,'title') && options.title
   title(ax2,options.title_label,'FontSize',options.FontSize)
   axis off;box off;
end


if isfield(options,'x_limits') && isfield(options,'y_limits')
    if ~isempty(options.x_limits) && ~isempty(options.y_limits)
        ax1.XLim=options.x_limits;
        ax1.YLim=options.y_limits;
        ax2.XLim=options.x_limits;
        ax2.YLim=options.y_limits;
    end
end





end





   