
function Axis00(ax)

%
% Sets the x and y axis at (0,0)
%
% 


if nargin==0
    ax=gca;
end

ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

end

