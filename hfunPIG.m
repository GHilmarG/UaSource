

function h = hfunPIG(x,y,hfunOpt)
	
	
	
	if numel(hfunOpt) >1
		
		F = TriScatteredInterp(hfunOpt(:,1),hfunOpt(:,2),hfunOpt(:,3),'natural');
		%h=h0/2+abs(x)*0.;
		h = h0*F(x,y);
		h(h<MeshSizeMin)=MeshSizeMin ; h(h>MeshSizeMax)=MeshSizeMax;
		
	elseif numel(hfunOpt) ==1
	    h =  hfunOpt(1)+abs(x)*0.;
	end
		
	
end

