function [UserVar,F]=GetGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined)

%
% Wrapper around DefineGeometry.m
%
% Calls DefineGeometry and does some error checks on the returned values.
%
% Only updates on s, b, B and S fields if these variables
% are contained in the string FieldsToBeDefined.
%
% Otherwise returns input fields.
%
%


nOut=nargout;
if nOut~=2
    error('Ua:GetGeometry','Need 2 output arguments')
end


if nargin<5 || isempty(FieldsToBeDefined)
    FieldsToBeDefined="-s-b-S-B-rho-rhow-g-";
end

% Note: GF can not be calculated without knowing both the geometrical variables (s,b,S,B) returned by
% DefineGeometry.m and the densities (rho, rhow) returned by DefineDensities.m

if ~isempty(F.GF)
    
    if ~isequal(numel(F.GF.node),MUA.Nnodes)
        error('InternalUAerror: numel(F.GF.node) ~= MUA.Nnodes')
    end
end


F.x=MUA.coordinates(:,1) ;  F.y=MUA.coordinates(:,2) ; 


if isfile("DefineGeometryAndDensities.m")
    
    [UserVar,sTemp,bTemp,STemp,BTemp,rhoTemp,rhowTemp,gTemp]=DefineGeometryAndDensities(UserVar,CtrlVar,MUA,F,FieldsToBeDefined);
    
    if contains(FieldsToBeDefined,'rho')
        [F.rho,F.rhow,F.g]=TestDensityInputValues(CtrlVar,MUA,rhoTemp,rhowTemp,gTemp);
        F.rho=double(F.rho) ; F.rhow=double(F.rhow) ; F.g=double(F.g) ; 
    end
    
    F.alpha=CtrlVar.alpha ;
    
else
    
    nArgs=nargin('DefineGeometry');
    
    switch nArgs
        case 5
            [UserVar,sTemp,bTemp,STemp,BTemp,F.alpha]=DefineGeometry(UserVar,CtrlVar,MUA,CtrlVar.time,FieldsToBeDefined);
        case 6
            [UserVar,sTemp,bTemp,STemp,BTemp,F.alpha]=DefineGeometry(UserVar,CtrlVar,MUA,CtrlVar.time,FieldsToBeDefined,F);
        otherwise
            error('DefineGeometry must have either 5 or 6 inputs arguments.')
    end
    
    
    [UserVar,F]=GetDensities(UserVar,CtrlVar,MUA,F);
    
end

% some error checks
errorStruct.identifier = 'GetGeometry:NaNinInput';

if contains(FieldsToBeDefined,'s')
    
    if any(isnan(sTemp))
        errorStruct.message = 's returned by DefineGeometry  contains nan.';
        error(errorStruct)
    end
    
    if ~isfinite(sTemp)
        errorStruct.message = 's returned by DefineGeometry  not a finite number.';
        error(errorStruct)
    end
    
    if isempty(sTemp)
        errorStruct.message = 's returned by DefineGeometry is empty.';
        error(errorStruct)
    end
    
    F.s=double(sTemp);
    if numel(F.s)==1 ; F.s=F.s+zeros(MUA.Nnodes,1) ; end
end

if contains(FieldsToBeDefined,'b')
    if any(isnan(bTemp))
        errorStruct.message = 'nan in b';
        error(errorStruct)
    end
    
    if ~isfinite(bTemp)
        errorStruct.message = 'b returned by DefineGeometry  not a finite number.';
        error(errorStruct)
    end
    
    if isempty(bTemp)
        errorStruct.message = 'b returned by DefineGeometry is empty.';
        error(errorStruct)
    end
    
    
    F.b=double(bTemp);
    if numel(F.b)==1 ; F.b=F.b+zeros(MUA.Nnodes,1) ; end
end

if contains(FieldsToBeDefined,'S')
    if any(isnan(STemp))
        errorStruct.message = 'nan in S';
        error(errorStruct)
    end
    
    
    if ~isfinite(STemp)
        errorStruct.message = 'S returned by DefineGeometry  not a finite number.';
        error(errorStruct)
    end
    
    if isempty(STemp)
        errorStruct.message = 'S returned by DefineGeometry is empty.';
        error(errorStruct)
    end
    
    
    F.S=double(STemp);
    if numel(F.S)==1 ; F.S=F.S+zeros(MUA.Nnodes,1) ; end
end

if contains(FieldsToBeDefined,'B')
    
    if any(isnan(BTemp))
        errorStruct.message = 'nan in B';
        error(errorStruct)
    end
    
    
    if ~isfinite(BTemp)
        errorStruct.message = 'B returned by DefineGeometry  not a finite number.';
        error(errorStruct)
    end
    
    if isempty(BTemp)
        errorStruct.message = 'B returned by DefineGeometry is empty.';
        error(errorStruct)
    end
    
    F.B=double(BTemp);
    if numel(F.B)==1 ; F.B=F.B+zeros(MUA.Nnodes,1) ; end
end






if contains(FieldsToBeDefined,'s')|| contains(FieldsToBeDefined,'b')
    F.h=F.s-F.b;
end


% Generally the min ice-thickness constraint should be enforced using the active set
% method. However, when the user defines s and b it seems reasonable to expect that
% initial user-defined ice distribution to be consistent with user-defined min ice
% thickness.
%

if contains(FieldsToBeDefined,'s')|| contains(FieldsToBeDefined,'b')
    CtrlVar.ResetThicknessToMinThickness=1;
end



switch CtrlVar.Calculate.Geometry
    
    case "bs-FROM-hBS"
        
        [F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);
        
    case "bh-FROM-sBS"
        
        
        [F.b,F.h,F.GF]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,F.B,F.S,F.rho,F.rhow) ;
        
    otherwise
        
        error('which case')
        
end



end