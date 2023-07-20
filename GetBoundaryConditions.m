
function [UserVar,BCs]=GetBoundaryConditions(UserVar,CtrlVar,MUA,BCs,F)


narginchk(5,5)
nargoutchk(2,2)

% Only allow the use of the (not so new anymore) DefineBoundaryConditions.m file

InputFile="DefineBoundaryConditions.m" ; TestIfInputFileInWorkingDirectory(InputFile) ;
InputFileNargIn=nargin(InputFile);
InputFileNargOut=nargout('DefineBoundaryConditions');

if CtrlVar.InfoLevel>=10

    fprintf(' Using DefineBoundaryConditions.m to define boundary conditions \n')
    
end



switch InputFileNargOut
    
    case 1
        
        if InputFileNargIn==5
            
            BCs=DefineBoundaryConditions(UserVar,CtrlVar,MUA,F,BCs);

        else

            BCs=DefineBoundaryConditions(UserVar,CtrlVar,MUA,BCs,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.GF);
                
        end
        
    case 2
        
        if InputFileNargIn==5
        
            [UserVar,BCs]=DefineBoundaryConditions(UserVar,CtrlVar,MUA,F,BCs);
        
        else

            [UserVar,BCs]=DefineBoundaryConditions(UserVar,CtrlVar,MUA,BCs,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.GF);
        
        end
        
end


switch lower(CtrlVar.FlowApproximation)
    
    case 'sstream'
        
        % sstream only uses the ud boundary conditions, so set ub to empty.
        
        BCs.udFixedNode=[];
        BCs.udFixedValue=[];
        BCs.vdFixedNode=[];
        BCs.vdFixedValue=[];
        BCs.udTiedNodeA=[];
        BCs.udTiedNodeB=[];
        BCs.vdTiedNodeA=[];
        BCs.vdTiedNodeB=[];
        BCs.udvdFixedNormalNode=[];
        BCs.udvdFixedNormalValue=[];
        
    case 'ssheet'
        
        % ssheet can now have both ud and ub boundary conditions, so do not set to empty
        
        %         BCs.ubFixedNode=[];
        %         BCs.ubFixedValue=[];
        %         BCs.vbFixedNode=[];
        %         BCs.vbFixedValue=[];
        %         BCs.ubTiedNodeA=[];
        %         BCs.ubTiedNodeB=[];
        %         BCs.vbTiedNodeA=[];
        %         BCs.vbTiedNodeB=[];
        %         BCs.ubvbFixedNormalNode=[];
        %         BCs.ubvbFixedNormalValue=[];
        
end


% do some basic checks, will not pick everything, but should find the most typical mistakes.

if islogical(BCs.hFixedNode)
   error('GetBoundaryConditions:h','BCs.hFixedNode is a logical variable, but must be numeric. \n')
end

if islogical(BCs.ubFixedNode)
   error('GetBoundaryConditions:ub','BCs.ubFixedNode is a logical variable, but must be numeric. \n')
end

if islogical(BCs.vbFixedNode)
   error('GetBoundaryConditions:vb','BCs.vbFixedNode is a logical variable, but must be numeric. \n')
end

if islogical(BCs.udFixedNode)
   error('GetBoundaryConditions:ub','BCs.udFixedNode is a logical variable, but must be numeric. \n')
end

if islogical(BCs.vdFixedNode)
   error('GetBoundaryConditions:vb','BCs.vdFixedNode is a logical variable, but must be numeric. \n')
end



if numel(BCs.ubFixedNode) ~=  numel(BCs.ubFixedValue)
    error('GetBoundaryConditions:ub','Number of fixed ub nodes not equal to number of ub fixed values! \n')
end


if numel(BCs.vbFixedNode) ~=  numel(BCs.vbFixedValue)
    error('GetBoundaryConditions:vb','Number of fixed vb nodes not equal to number of vb fixed values! \n')
end

if numel(BCs.hFixedValue) ~=  numel(BCs.hFixedNode)
    error('GetBoundaryConditions:h','Number of fixed h nodes not equal to number of h fixed values! \n')
end

if numel(BCs.LSFFixedNode) ~=  numel(BCs.LSFFixedValue)
    error('GetBoundaryConditions:LSF','Number of fixed LSF nodes not equal to number of LSF fixed values! \n')
end

%% Check for duplicate nodes in boundary conditions and do some basic corrections (although is this not really the job of the modeller...?)
%  Also, I am not checking for all possible such cases and combinations.

if numel(BCs.ubTiedNodeA) ~= numel(BCs.ubTiedNodeB)
    
    fprintf('# ub A ties =%i \n ',  numel(BCs.ubTiedNodeA))
    fprintf('# ub B ties =%i \n ',  numel(BCs.ubTiedNodeB))
    
    error("Ua:BoundaryConditions","Number of A and B ub ties not the same")
    
end

if numel(BCs.vbTiedNodeA) ~= numel(BCs.vbTiedNodeB)
    
    fprintf('# vb A ties =%i \n ',  numel(BCs.vbTiedNodeA))
    fprintf('# vb B ties =%i \n ',  numel(BCs.vbTiedNodeB))
    
    error("Ua:BoundaryConditions","Number of A and B vb ties not the same")
    
end

%% This is a questionable approach.
%  Even if a node is twice in list A, it can be unique if it is tied
%  to two different freedome of degrees
% if ~isempty(BCs.ubTiedNodeA)
%     % eliminate dublications in nodal ties
%     % this needs to be done for both A and B ties, and then
%     % dublicates in A deleted from B, and duplicates in B deleted from A.
%     [BCs.ubTiedNodeA,ia]=unique(BCs.ubTiedNodeA);
%     BCs.ubTiedNodeB=BCs.ubTiedNodeB(ia) ;
%
%     [BCs.ubTiedNodeB,ia]=unique(BCs.ubTiedNodeB);
%     BCs.ubTiedNodeA=BCs.ubTiedNodeA(ia) ;
%
% end
%
% if ~isempty(BCs.vbTiedNodeA)
%     % eliminate dublications in nodal ties
%     [BCs.vbTiedNodeA,ia]=unique(BCs.vbTiedNodeA);
%     BCs.vbTiedNodeB=BCs.vbTiedNodeB(ia) ;
%
%     [BCs.vbTiedNodeB,ia]=unique(BCs.vbTiedNodeB);
%     BCs.vbTiedNodeA=BCs.vbTiedNodeA(ia) ;
%
% end
%%

if ~isempty(BCs.ubFixedNode)
    [BCs.ubFixedNode,ia]=unique(BCs.ubFixedNode);
    BCs.ubFixedValue=BCs.ubFixedValue(ia) ;
end

if ~isempty(BCs.vbFixedNode)
    [BCs.vbFixedNode,ia]=unique(BCs.vbFixedNode);
    BCs.vbFixedValue=BCs.vbFixedValue(ia) ;
end

if ~isempty(BCs.hFixedNode)
    [BCs.hFixedNode,ia]=unique(BCs.hFixedNode);
    BCs.hFixedValue=BCs.hFixedValue(ia) ;
end

if ~isempty(BCs.vbFixedNode) && ~isempty(BCs.vbTiedNodeA)
    % If a degree of freedom is both 'fixed' and 'tied', only 'fix' it and get rid of the 'tie'
    [BCs.vbTiedNodeA,ia]=setdiff(BCs.vbTiedNodeA,BCs.vbFixedNode) ;
    BCs.vbTiedNodeB=BCs.vbTiedNodeB(ia);
    
end

if ~isempty(BCs.ubFixedNode) && ~isempty(BCs.ubTiedNodeA)
    % If a degree of freedom is both 'fixed' and 'tied', only 'fix' it and get rid of the 'tie'
    [BCs.ubTiedNodeA,ia]=setdiff(BCs.ubTiedNodeA,BCs.ubFixedNode) ;
    BCs.ubTiedNodeB=BCs.ubTiedNodeB(ia);
    
end
%%

if CtrlVar.doplots && CtrlVar.PlotBCs
    fig=FindOrCreateFigure("Boundary Conditions");
    clf(fig)
    hold off
    PlotBoundaryConditions(CtrlVar,MUA,BCs,'k');
end




end





