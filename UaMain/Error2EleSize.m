function h=Error2EleSize(CtrlVar,e,e0,hMin,hMax,p)

%%
%   Assuming I want same errors everywhere, i.e. equidistribution of errors, and
%   assuming the errors vary with ele size h as
%
%  Errors = e h^p
%
%  where e is some indirect estimate of errors
%
%  then
%
%   (e+e0) (h-hMin)^p= K
%
%    h=  (K/(e+e0))^(1/p)  + hMin
%
%    hMax= (K/e0))^(1/p)) + hMin
%
%    K=e0 (hMax-hMin)^p
%
%    ->  h =  (e0/(e+e0))^(1/p) (hMax-hMin)  + hMin
%
%
%  For e=e0
%
% h =  (1/2)^(1/p) (hMax-hMin)  + hMin
%
% for p=1 : h(e0) =  (hMax+hMin)/2
% for p=2 : h(e0) =  1/4 (hMax-hMin)  + hMin
%
%
%%

if isempty(hMin)
    hMin=CtrlVar.MeshSizeMin;
end


if isempty(hMax)
    hMax=CtrlVar.MeshSizeMax;
end


if isempty(p)
    p=1;
end

e=e(:) ;

h=hMin+(e0./(e+e0)).^(1/p).*(hMax-hMin);
h=h(:);


end