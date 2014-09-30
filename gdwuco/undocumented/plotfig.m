function [ ] = plotfig( save_name,param )
%PLOTFIG Ploting figures with optimal size for paper.
%   Usage: plotfig(save_name);
%          plotfig(save_name,param);
%   
%   Input parameters:
%       save_name: name to save the figure
%       param   : optional parameters
%
%   *param* a Matlab structure containing the following fields:
%
%   * *param.pathfigure* : path to the folder to save the figures
%   * *param.legendlocation* : location of the figure (default 'Best');
%   * *param.position* : position and size of the figure 
%     (default [100 100 600 400])
%   * *param.labelsize* : Size of the label (default 12)
%   * *param.titlesize* : Size of the title (default 14)
%   * *param.titleweight*: Weight of the title (default 'normal')
%   * *param.save*: Save the figure (default 1)


% Nathanael Perraudin
% 19 June 2013

if nargin<2
    param=struct;
end

% Optional parameters
if ~isfield(param, 'pathfigure'), param.pathfigure = 'figures/'; end
if ~isfield(param, 'position'), param.position = [100 100 600 400]; end
if ~isfield(param, 'labelsize'), param.labelsize = 14; end
if ~isfield(param, 'titlesize'), param.titlesize = 16; end
if ~isfield(param, 'titleweight'), param.titlesize = 'normal'; end
if ~isfield(param, 'save'), param.save =1; end
if ~isfield(param, 'baw'), param.baw =0; end

if param.baw
    set_baw_color();
end



% set the axes
try  %#ok<TRYNC>
    set(gcf, 'Position', param.position);
    set(gcf,'PaperPositionMode','auto');
end

% set the title
try %#ok<TRYNC>
    set(gca,'FontSize',param.labelsize)
end
try %#ok<TRYNC>
    h=get(gca,'Title');
    set(h,'FontSize',param.titlesize);
    set(h,'FontWeight',param.titleweight);
end
try %#ok<TRYNC>
    h=get(gca,'xlabel');
    set(h,'FontSize',param.labelsize);
end
try %#ok<TRYNC>
    h=get(gca,'ylabel');
    set(h,'FontSize',param.labelsize);
end
drawnow;

% save the results


if param.save
    if ~isdir(param.pathfigure)
           mkdir(param.pathfigure);
    end
    filename=strcat(param.pathfigure,save_name);
    print('-dpng','-zbuffer','-r300',[filename,'.png']);
    %hgsave([filename,'.fig']);
end
end

