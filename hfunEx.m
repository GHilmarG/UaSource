
function h = hfunEx(x,y,MeshSize,MeshSizeMin,MeshSizeMax)
	

	
	
	if numel(MeshSize)== 1 ; h0=MeshSize ; end
	
	if numel(MeshSize) >1
		
		
		F = TriScatteredInterp(MeshSize(:,1),MeshSize(:,2),MeshSize(:,3),'natural');
		%h=h0/2+abs(x)*0.;
		h = F(x,y);
		h(h<MeshSizeMin)=MeshSizeMin ; h(h>MeshSizeMax)=MeshSizeMax;
		
	else
		h =  h0+abs(x)*0.;
	end
	
end

