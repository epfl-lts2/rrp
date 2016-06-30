function ai_plot_graph(G,transitions,hole_interval,param)
%AI_PLOT_GRAPH Plot the transitions graph
%   Usage: ai_plot_graph(G, transition, sthole, finhole, fs);
%          ai_plot_graph(G, transition, sthole, finhole, fs, param);
%   
%   Input parameters:
%       G           : graph of transitions
%       transitions : transitions to be highlighted
%       hole_interval : start and end position of the hole (in seconds)
%       param       : optional structure of parameters;
%   
%   This function plots the graph of transitions and highlights the to
%   selected transitions with arrows.
%

% Authors: Nathanael Perraudin, Nicki Hollighaus
% Date   : June 2016

if nargin<4
    param = struct;
end
time = G.time;
fs = G.fs;
if isfield(param,'subgraph') && isfield(param.subgraph,'time_in')
    sthole = param.subgraph.time_in; % start of the hole
    finhole = param.subgraph.time_fin; % final of the hole    
else
    sthole = hole_interval(1)*fs*G.ratio; % start of the hole
    finhole = hole_interval(2)*fs*G.ratio; % final of the hole
end
sig = ones(G.N,1);
sig(time>sthole & time < finhole) = 0;
if isfield(param,'subgraph')
    sig(time<sthole & time > (sthole - param.subgraph.time_range * fs * G.ratio)) = 2;
    sig(time>finhole & time < (finhole + param.subgraph.time_range * fs * G.ratio)) = 2;
end


gsp_plot_signal(G,sig);

hold on
%     plot([G.coords(jumps(1,1),1);G.coords(jumps(1,2),1)'],...
%         [G.coords(jumps(1,1),2);G.coords(jumps(1,2),2)'],...
%         G.plotting.edge_style, 'LineWidth',2,...
%         'Color','k');
% 
%     plot([G.coords(jumps(2,1),1);G.coords(jumps(2,2),1)'],...
%         [G.coords(jumps(2,1),2);G.coords(jumps(2,2),2)'],...
%         G.plotting.edge_style, 'LineWidth',2,...
%         'Color','k');
quiver(G.coords(transitions(1,1),1),G.coords(transitions(1,1),2),...
    G.coords(transitions(1,2),1) - G.coords(transitions(1,1),1),...
    G.coords(transitions(1,2),2) - G.coords(transitions(1,1),2),...
    G.plotting.edge_style, 'LineWidth',3,...
    'Color','k','Autoscale','off','MaxHeadSize',0.5);
quiver(G.coords(transitions(2,1),1),G.coords(transitions(2,1),2),...
    G.coords(transitions(2,2),1) - G.coords(transitions(2,1),1),...
    G.coords(transitions(2,2),2) - G.coords(transitions(2,1),2),...
    G.plotting.edge_style, 'LineWidth',3,...
    'Color','k','Autoscale','off','MaxHeadSize',0.5);
colorbar off

end