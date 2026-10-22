function  [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

%  Fully-Implicit solution for uvh (SSTREAM) and h (SSHEET)

narginchk(8,8)
nargoutchk(4,5)






switch lower(CtrlVar.FlowApproximation)


    case "sstream"

        switch CtrlVar.uvhSolver


            case "uvh_GaussNewton_fsolve"
                % This is a Newton-Gauss approach, and was found to be much slower than the "Root-Finding Newton"
                % approach. This is positioned INSIDE the active-set loop. A valuable experience was gained, inasmuch it was
                % established that the Newton approach is far better in terms of performance.
                [UserVar,RunInfo,F1,l1,BCs1]=uvh_GaussNewton_fsolve(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

            case "Root-Finding Newton"

                [UserVar,RunInfo,F1,l1,BCs1]=SSTREAM_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

            otherwise
                
                error("uvh2D:UnknownSolver","Unknown CtrlVar.uvhSolver=%s",CtrlVar.uvhSolver)
        end

        %  [UserVar,RunInfo,F1,l1,BCs1]=SSTREAM_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

        F1.ud=zeros(MUA.Nnodes,1)  ; F1.vd=zeros(MUA.Nnodes,1);

    case 'ssheet'

        [UserVar,RunInfo,F1,l1,BCs1]=SSHEET_TransientImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

    case "uvhprescribed"


        [UserVar,RunInfo,F1,l1,BCs1]=uvhPrescibed(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

    otherwise

        error('what case')
        
end


end


