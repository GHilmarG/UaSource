function ColorIndex=Variable2ColorIndex(Variable)

% simple linear mapping from values of a variable to a linear index
% having same number of elements as the variable
%
% Useful for color mapping.
%



if (max(Variable)-min(Variable))<eps
    
    ColorIndex=Variable*0+1;
    
else
    
    % This is, according to the matlab help pages, exactly how matlab maps between values and colors in the colormap by default.
    
    ColorIndex = fix((Variable-min(Variable))/(max(Variable)-min(Variable))*numel(Variable))+1;
    ColorIndex(ColorIndex<1) = 1;
    ColorIndex(ColorIndex>numel(Variable)) = numel(Variable);
    
end


end