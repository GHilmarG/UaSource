function [UserVar,F]=GetSeaIceParameters(UserVar,CtrlVar,MUA,F)

narginchk(4,4)
nargoutchk(2,2)

if CtrlVar.IncludeMelangeModelPhysics

    InputFile="DefineSeaIceParameters.m";
    TestIfInputFileInWorkingDirectory(InputFile) ;
    NargIn=nargin(InputFile);
    NargOut=nargout(InputFile);

    if NargIn ~= 4 || NargOut ~= 9

        fprintf(" The function DefineSeaIceParameters must be on the form: \n ")
        fprintf("[UserVar,uo,vo,Co,mo,ua,va,Ca,ma]=DefineSeaIceParameters(UserVar,CtrlVar,MUA,F) \n ")
        error("Ua:GetSeaIceParameters","DefineSeaIceParameters.m must have 4 input and 9 output arguments")

    end



    [UserVar,F.uo,F.vo,F.Co,F.mo,F.ua,F.va,F.Ca,F.ma]=DefineSeaIceParameters(UserVar,CtrlVar,MUA,F);



    %% Test values

    if isempty(F.uo)
        error('Ua:GetSeaIceParameters','Ocean velocity component uo is empty.')
    end
    
    
    if isempty(F.vo)
        error('Ua:GetSeaIceParameters','Ocean velocity component vo is empty.')
    end
    
    
    if all(F.mo==0)
        
        error('Ua:GetSeaIceParameters','Ocean/ice stress exponent (mo) is zero, but must have a non-zero value. ')
        
    end
    
    
    if all(F.ma==0)
        
        error('Ua:GetSeaIceParameters','atmo/ice stress exponent (ma) is zero, but must have a non-zero value. ')
        
    end
    
    
    if isempty(F.Co)
        error('Ua:GetSeaIceParameters','Ocean/ice slipperiness coefficient Co is empty.')
    end
    
    
    if isempty(F.mo)
        error('Ua:GetSeaIceParameters','Ocean/ice drag stress-exponent mo is empty.')
    end
    
    
    
    if isempty(F.Ca)
        error('Ua:GetSeaIceParameters','Atmo/ice slipperiness coefficient Ca is empty.')
    end
    
    
    
    if isempty(F.ma)
        error('Ua:GetSeaIceParameters','Atmo/ice drag stress-exponent ma is empty.')
    end
    
    CtrlVar.SlidingLaw="Weertman" ;  
    
    % The "sliding law" used here is only applied over the floating section, and provides a link between the basal drag and ocean
    % velocities over the floating parts of the domain. This drag is calculated using a Weertman type law. Hence, here the test
    % needs to be done for a Weertman sliding law. The sliding law used for the grounded parts can, and in general will, be
    % different. The sliding law for the grounded section is selected, as always, by defining the value of the variable
    % CtrlVar.SlidingLaw appropriately in DefineInitialInputs.m
                                     

    [F.Co,F.mo]=TestSlipperinessInputValues(CtrlVar,MUA,F.Co,F.mo);
    [F.Ca,F.ma]=TestSlipperinessInputValues(CtrlVar,MUA,F.Ca,F.ma);
    
    
    if isscalar(F.uo)
        F.uo=F.uo+zeros(MUA.Nnodes,1);
    end
    
    if isscalar(F.vo)
        F.vo=F.vo+zeros(MUA.Nnodes,1);
    end
    
    if isscalar(F.ua)
        F.ua=F.ua+zeros(MUA.Nnodes,1);
    end
    
    if isscalar(F.va)
        F.va=F.va+zeros(MUA.Nnodes,1);
    end
    
    
    if all(F.Co == 0)
        
        error('GetSeaIceParameters:CoIszero','Co can not be zero')
        
    end
    
    if all(F.Ca== 0)
        
        error('GetSeaIceParameters:CaIsZero','Ca can not be zero')
        
    end
    
    
    
end



end