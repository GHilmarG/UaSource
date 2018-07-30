function [UserVar,s,b,S,B,alpha]=GetGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

% 
% Wrapper around DefineGeometry.m
%
% Calls DefineGeometry and does some error checks on the returned values.
%
% Only passes on s, b, B and S (as returned by the call to DefineGeomety) if these variables 
% are contained in the string FieldsToBeDefined.
%
% Otherwise returns empty fields. 
% 
%


error(' No longer to be used')

nOut=nargout;
if nOut~=6
    error('Ua:GetGeometry','Need 6 output arguments')
end


s=[] ; b=[] ; B=[] ; S=[] ; 

if nargin<5 || isempty(FieldsToBeDefined)
    FieldsToBeDefined='sbSB';
end





[UserVar,sTemp,bTemp,STemp,BTemp,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined);

% if FieldsToBeDefined==""
%     return
% end


% some error checks
errorStruct.identifier = 'GetGeometry:NaNinInput';

if contains(FieldsToBeDefined,'s')
    if any(isnan(sTemp))
        errorStruct.message = 'nan in s';
        error(errorStruct)
    end
    s=sTemp;
end

if contains(FieldsToBeDefined,'b')
    if any(isnan(bTemp))
        errorStruct.message = 'nan in b';
        error(errorStruct)
    end
    b=bTemp;
end

if contains(FieldsToBeDefined,'S')
    if any(isnan(STemp))
        errorStruct.message = 'nan in S';
        error(errorStruct)
    end
    S=STemp;
end

if contains(FieldsToBeDefined,'B')
    if any(isnan(BTemp))
        errorStruct.message = 'nan in B';
        error(errorStruct)
    end
    B=BTemp;
end







end