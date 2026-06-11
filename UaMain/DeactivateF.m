




function Fnew=DeactivateF(CtrlVar,MUA,Fold,k)

%%
%
%  Creates a new F from a previous one using the mapping
%
%   Fnew=Fold(k)
%
% This can be used, for example, after element deactivation to create a corresponding "deactivated" F field.
%
%
% SEE ALSO: DeactitivateMUAelements.m
%%


Fnew=UaFields ;

Fields=fieldnames(Fold) ;

%% loop over all fields
for I=1:numel(Fields)

    if isempty(Fold.(Fields{I}))

        Fnew.(Fields{I})=[];

    elseif isnumeric(Fold.(Fields{I}))



        if isscalar(Fold.(Fields{I}))

            Fnew.(Fields{I})=Fold.(Fields{I}) ;

        else

            Fnew.(Fields{I})=Fold.(Fields{I})(k) ;

        end

    end

end

% What to do about GF?  Best to let the user re-calculate GF as this might include some adjustments to s and b. But I can
% cover a very typical case here which is to remap just the GF.node variable.

if ~isempty(Fold.GF.node)
    Fnew.GF.node=Fold.GF.node(k);
end




end



