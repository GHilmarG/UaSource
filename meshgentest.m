
function [coo,connectivity,TR]=meshgentest
node=[-1 -1 ; 1 -1 ; 1 1 ; -1 1];

hdata.fun = @hfun1;

options.output=false;
[coo,connectivity] = mesh2d(node,[],hdata,options);


% create triangulation object
TR = TriRep(connectivity, coo(:,1),coo(:,2));


triplot(TR) ; hold on
% to find coordinates of TR number 3  use;
tri=TR.X(TR(3,:),:);
scatter(tri(:,1),tri(:,2),'MarkerFaceColor','y');

tri=TR.X(TR(1,:),:);
scatter(tri(:,1),tri(:,2),'MarkerFaceColor','c');

tri=TR.X(TR(2,:),:);
scatter(tri(:,1),tri(:,2),'MarkerFaceColor','g');

FB=freeBoundary(TR) ; 

x=TR.X(:,1) ; y=TR.X(:,2);
plot(x(FB),y(FB),'r','LineWidth',3)

% find nodes along x=1
ind=(x(FB(:,1))  == 1 & x(FB(:,2))==1);
NeumannEdge=[FB(ind,1) FB(ind,2)];

plot(x(NeumannEdge),y(NeumannEdge),'g','LineWidth',3)
NeumannEle=edgeAttachments(TR,NeumannEdge)

for II=1:length(TR.Triangulation)
   
    centre=mean(TR.X(TR(II,:),:));
    
    text(centre(1),centre(2),num2str(II))
end
    
    
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = hfun1(x,y)

% User defined size function for square

%h = 0.01 + 0.1*sqrt( (x-0.25).^2+(y-0.75).^2 );
h = 0.5 +abs(x)*0;

end      % hfun1()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

