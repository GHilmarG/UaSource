

clear all ; close all
ndim=1; nip=6; nod=4;


[points,weights]=sampleEdge('line',nip,ndim);
fun=zeros(nip,nod) ; der=zeros(nip,nod) ;

for Iint=1:nip
	fun(Iint,:)=shape_funEdge(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
	der(Iint,:)=shape_derEdge(Iint,ndim,nod,points);
end


figure ; hold on
for I=1:nod
	hold on
	plot(points,fun(:,I),'b')
	plot(points,der(:,I),'r')
end
hold off

[points,weights]=sample('line',nip,ndim);
fun=zeros(nip,nod) ; der=zeros(nip,nod) ;

for Iint=1:nip
	fun(Iint,:)=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
	der(Iint,:)=shape_der(Iint,ndim,nod,points);
end


figure ; hold on
for I=1:nod
	hold on
	plot(points,fun(:,I),'b')
	plot(points,der(:,I),'r')
end

