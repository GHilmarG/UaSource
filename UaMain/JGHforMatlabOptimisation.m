function varargout=JGHforMatlabOptimisation(func,p)
%%
% Wrapper around the cost function for the MATLAB optimisation toolbox (fmincon, fminunc). (23 November 2026)
%
% fmincon (interior-point) may evaluate the cost function at trial points with extreme parameter values (for example when
% the bounds CtrlVar.Cmin, Cmax are very wide), at which the forward problem cannot be solved (eg the Cholesky factorisation in
% the velocity solve fails). An error then stops the inversion. Instead, J=Inf is returned, and fmincon rejects the trial point
% and tries a shorter step (the interior-point algorithm can recover from Inf or NaN values of the objective).
%
%   varargout=JGHforMatlabOptimisation(func,p)
%
try
    [varargout{1:nargout}]=func(p);
catch ME
    fprintf("JGHforMatlabOptimisation: the cost function could not be evaluated at a trial point (%s: %s); J=Inf returned.\n",ME.identifier,ME.message)
    varargout{1}=Inf;
    if nargout>1
        varargout{2}=nan(size(p));
    end
    for k=3:nargout
        varargout{k}=[];
    end
end
end
