function SetLabels(CtrlVar,x,y,z)

% SetLabels("km","km","m")

xlabel("x ("+x+")")

if nargin>3
    ylabel("y ("+y+")")
end

if nargin==4
    zlabel("z ("+z+")")
end

end