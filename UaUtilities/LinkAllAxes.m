

function LinkAllAxes(FigHandle)


%%
%
% An attempt to ensure all axes of a given figure are linked and in the same position.
%
% Reason: I find that despite having already linked the axes, the still do not always align after resizing the figure. Not
% sure what I am doing wrong, if anything, but then re-linking using this utility seems to do the job.
%
%%

if nargin==0 || isempty(FigHandle)
    FigHandle=gcf;
end

axHandles = findall(FigHandle, 'Type', 'axes') ;

linkaxes(axHandles) ; % this should do it for all the axes, by default this works for xy 

axHandles(1).Visible="on"; % it is better to have this on for all axes, in that way one can see if the axes are really overlapping and identical.
                           % afterwards this can be turned off

for I=2:numel(axHandles)
    axHandles(I).Position = axHandles(1).Position;
    axHandles(I).Visible="on";
end

%%