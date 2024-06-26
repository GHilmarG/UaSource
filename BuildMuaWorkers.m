




function MUAworkers=BuildMuaWorkers(CtrlVar,MUA,MUAworkers)


if ~CtrlVar.Parallel.BuildWorkers  % Somtimes building workers is suppressed.
                                   % For example when adapting the mesh repeatedly, there is no need to build workers for each
                                   % mesh. Generally, only build workers ahead of a uv of uvh solve

    
    return

end



if isempty(CtrlVar.Parallel.uvhAssembly.spmd.nWorkers)
    poolobj = gcp ;
    CtrlVar.Parallel.uvhAssembly.spmd.nWorkers=poolobj.NumWorkers;
end

nW=CtrlVar.Parallel.uvhAssembly.spmd.nWorkers;

%%





%% create element lists for each partition
% Use round to make sure that there are exactly nW partitions
% and ensure that all elements are included.

Partition=cell(nW,1);
N=round(MUA.Nele/nW) ; i1=1 ; i2=N;
for iWorker=1:(nW-1)
    Partition{iWorker}=i1:i2 ;
    i1=i2+1 ;
    i2=i2+N ;
end

i2=MUA.Nele;
Partition{nW}=i1:i2 ;

% outside of spmd  M is  composite
% inside of spmd M is struct on each worker

MUA.dM=[] ;

% Have I already build the workers for this mesh?
% I test this by seeing if the previous partition is the same 

% Not clear to me how to check if composite is in correct state.
% The only option seems to be just to try to access it and see if it produces an error.

BuildWorkers= false ; 


if isempty(MUAworkers) || numel(MUAworkers)==0  % if all empty, rebuild

    BuildWorkers= true;

else

    try
        isCompoositeInCorrectState=MUAworkers{1}.Nnodes==MUA.Nnodes ;
    catch
        isCompoositeInCorrectState=false;  % if not in correct state, rebuild

    end

    if ~isCompoositeInCorrectState
        BuildWorkers=true;
    else  % Compoosite is tecknically in correct state, but may not contain the right elements

        % check if all partitions already contain the correct elements, if so then there is no need to build the workers anew
        spmd (0,nW)

            T=all(MUAworkers.Partition==Partition{spmdIndex}) ;
            Tand=spmdReduce(@and,T,1) ;

        end

       
        if ~Tand{1}
             % Does not contain the correct elements, so rebuild
            BuildWorkers=true;
        end
    end

    % Turns out the spmd version above is much faster
    % if isCompoositeInCorrectState
    %
    %     areWorkersAlreadyBuild=true;
    %     for iWorker=1:nW
    %
    %         % % do a simple test first
    %         % if numel(MUAworkers{iWorker}.Partition)~=numel(Partition{iWorker})
    %         %     areWorkersAlreadyBuild=false ;
    %         %     break
    %         % end
    %
    %         % now check if exacly the same elements are part of each partition
    %         if ~all(MUAworkers{iWorker}.Partition==Partition{iWorker})
    %             % if sum(MUAworkers{iWorker}.Partition)~=sum(Partition{iWorker})
    %             areWorkersAlreadyBuild=false ;  % if new partition is not equal to previous one, rebuild
    %
    %             break
    %         end
    %
    %     end
    %
    %     if ~areWorkersAlreadyBuild
    %         BuildWorkers=true;
    %     end
    %
    % end

end






if BuildWorkers
    MUA.workers=[];   % Make sure to get rid of this composite, within the MUA structure ahead of the spmd call
    % First of all this field is no longer needed as it has be judged to be in an incorrect state and needs to be recreated.
    % Second, keeping it means that this composite, ie MUA.workers, is used within the body of the spmd statement and this is not
    % allowed in matlab and will create a warning.
    spmd (0,nW)

        % Build M directly on the workers to avoid communication

        MUAworkers.nod=MUA.nod;
        MUAworkers.nip=MUA.nip;

        MUAworkers.Nnodes=MUA.Nnodes;
        MUAworkers.points=MUA.points;
        MUAworkers.weights=MUA.weights;

        MUAworkers.coordinates=MUA.coordinates;
        %
        MUAworkers.Partition=Partition{spmdIndex} ;
        MUAworkers.connectivity=MUA.connectivity(Partition{spmdIndex},:);
        MUAworkers.Nele=numel(Partition{spmdIndex});
        MUAworkers.Deriv=MUA.Deriv(Partition{spmdIndex},:,:,:);
        MUAworkers.DetJ=MUA.DetJ(Partition{spmdIndex},:);
        MUAworkers.EleAreas=MUA.EleAreas(Partition{spmdIndex}) ;


    end
end



end