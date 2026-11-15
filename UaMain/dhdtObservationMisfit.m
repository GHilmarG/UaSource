function [J,dJdu,dJdv,dJdh,Mats]=dhdtObservationMisfit(CtrlVar,MUA,F,Meas,Request)

%%
%
%   [J,dJdu,dJdv,dJdh,Mats]=dhdtObservationMisfit(CtrlVar,MUA,F,Meas,Request)
%
% dh/dt misfit term for the "-points-" and "-blocks-" formulations (CtrlVar.Inverse.DataMisfit.dhdt.Method), in which the
% modelled dh/dt is brought to the data:
%
% $$ J = \frac{1}{2} \mathbf{r}^T \mathbf{W} \mathbf{r} , \qquad \mathbf{r} = \mathbf{A} \dot{\mathbf{h}} - \mathbf{d} , \qquad \mathbf{W}=\mathrm{diag}(1/\sigma_k^2) $$
%
% where hdot is the modelled dh/dt at a set of points (dhdtAtPoints.m), and
%
%   "-points-" : the points are the measurement locations, and A is the identity (no averaging)
%   "-blocks-" : the points are the integration points, and A (Meas.dhdtA) averages over the blocks (footprints) of the measurements
%
% Errors are treated as uncorrelated, and there is no normalisation with the area of the domain. With g = A' W r, the gradients are
%
%   dJ/du = Gu' g ,  dJ/dv = Gv' g ,  dJ/dh = Gh' g
%
% Request (optional) defines which matrices are returned in Mats:
%
%   ""           : none (default)
%   "Jacobians"  : Mats.AGu, Mats.AGv, Mats.AGh (=A*Gu etc), and Mats.w (=1./sigma.^2), for the Gauss-Newton Hessian terms
%                  (AG)' W (AG) in Jqq.m, JhdotBB.m and JhdotBq.m
%   "Second"     : as "Jacobians", and in addition the residual terms Mats.Rhu, Mats.Rhv (see dhdtAtPoints.m) for JhdotBq.m
%
% see also: dhdtAtPoints.m, EvaluateJhdotAndDerivatives.m, CreateDhdtPointMeasurements.m, CreateDhdtBlockMeasurements.m
%
%%

narginchk(4,5)

if nargin<5 || isempty(Request)
    Request="";
end

Method=dhdtMisfitMethod(CtrlVar);

switch Method
    case "-points-"
        P=Meas.dhdtP;
        A=[];
    case "-blocks-"
        P=ShapeFunctionsAtIntegrationPoints(CtrlVar,MUA);
        A=Meas.dhdtA;
    otherwise
        error('dhdtObservationMisfit:Method','This function is not used for CtrlVar.Inverse.DataMisfit.dhdt.Method=%s',Method)
end

d=Meas.dhdtData(:);
w=1./Meas.dhdtErr(:).^2;

[hdot,Gu,Gv,Gh]=dhdtAtPoints(MUA,F,P);

if isempty(A)
    r=hdot-d;
    g=w.*r;
else
    r=A*hdot-d;
    g=A'*(w.*r);
end

J=full(r'*(w.*r))/2;

dJdu=full(Gu'*g);
dJdv=full(Gv'*g);
dJdh=full(Gh'*g);

Mats=[];

if Request=="Jacobians" || Request=="Second"
    if isempty(A)
        Mats.AGu=Gu ; Mats.AGv=Gv ; Mats.AGh=Gh ;
    else
        Mats.AGu=A*Gu ; Mats.AGv=A*Gv ; Mats.AGh=A*Gh ;
    end
    Mats.w=w;
end

if Request=="Second"
    [~,~,~,~,Mats.Rhu,Mats.Rhv]=dhdtAtPoints(MUA,F,P,g);
end

end
