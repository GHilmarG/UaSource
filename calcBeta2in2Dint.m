function [beta2int,Dbeta2Duuint,Dbeta2Dvvint,Dbeta2Duvint] = calcBeta2in2Dint(uint,vint,Cint,mint,Heint,CtrlVar,uoint,voint)
    
    % calculates beta^2 and D beta^2/Du at integration points
    % must be called within an integration point loop over all elements
    % Sets all values to zero where floating using input values of Heint 
    %

    %Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin; % for higher order elements it is possible for Cint to be less than zero  even for C strictly positiv
    %beta2int=Cint.^(-1/m).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1/m-1) ;
    
    if nargin<7
        U=uint; 
        V=vint;
    else
        U=uint-uoint;
        V=vint-voint;
    end
    
    % To do: 
    %
    % Allow for different slipperiness parameters for 
    % beta2= Heint 
    
    
    beta2int=(Cint+CtrlVar.Czero).^(-1./mint).*(sqrt(U.*U+V.*V+CtrlVar.SpeedZero^2)).^(1./mint-1) ;
    
    beta2int=Heint.*beta2int;
    
    
    
    if nargout>1
        
        % The directional derivative is
        % D beta^2(u,v)[Delta u, Delta v]= (1/m-1) C^(-1/m) (u^2+v^2)^((1-3m)/2m)  (u \Delta u + v \Delta v)
        %
        
        %Dbeta2int=(1/m-1).*Cint.^(-1/m).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1-3*m)/(2*m));
        Dbeta2int=(1./mint-1).*(Cint+CtrlVar.Czero).^(-1./mint).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1-3*mint)./(2*mint));
        
        
        Dbeta2int=Dbeta2int.*Heint;
        
        Dbeta2Duuint=Dbeta2int.*U.*U;
        Dbeta2Dvvint=Dbeta2int.*V.*V;
        Dbeta2Duvint=Dbeta2int.*U.*V;
        
        if any(isnan(Dbeta2int)) ; save TestSave  ;  error(' NaN in Dbeta2int ' ) ; end
        
    end
    
    if any(isnan(beta2int)) ; save TestSave  ;  error(' NaN in beta2int ' ) ; end
    
end