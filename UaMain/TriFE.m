function tri=TriFE(connectivity)

%% Creates 3 node triangle connectivity from FE connectivity.
%
% tri=TriFE(connectivity)
%
% On return tri is a 3-node connectivity matrix
% that uses all nodes in connectivity and is based on identical nodal numbering sceme
% as connectivity.
%
% If on input connectivity is that of 3-node elements, same connectivity is returned.
% If on input connectivity is that of 6 or 10 node elements, each element is split into
% 4 or 9 elements, respectively.
%


[Nele,nod]=size(connectivity) ;


switch nod
    
    case 3
        
        tri=connectivity;
        
    case 6
        
        tri=zeros(4*Nele,3);
        tri(1:Nele,:)=connectivity(:,[1 2 6]);
        tri(Nele+1:2*Nele,:)=connectivity(:,[2 3 4]);
        tri(2*Nele+1:3*Nele,:)=connectivity(:,[4 5 6]);
        tri(3*Nele+1:4*Nele,:)=connectivity(:,[2 4 6]);
        
    case 10
        
        tri=zeros(9*Nele,3);
        tri(1:Nele,:)=connectivity(:,[1 2 9]);
        tri(Nele+1:2*Nele,:)  =connectivity(:,[2 10 9]);
        tri(2*Nele+1:3*Nele,:)=connectivity(:,[2 3 10]);
        tri(3*Nele+1:4*Nele,:)=connectivity(:,[3 5 10]);
        tri(4*Nele+1:5*Nele,:)=connectivity(:,[3 4 5]);
        tri(5*Nele+1:6*Nele,:)=connectivity(:,[5 6 10]);
        tri(6*Nele+1:7*Nele,:)=connectivity(:,[6 7 8]);
        tri(7*Nele+1:8*Nele,:)=connectivity(:,[8 10 6]);
        tri(8*Nele+1:9*Nele,:)=connectivity(:,[8 9 10]);
        
    otherwise
        error('Ua:TriFE','what case?')
end


end