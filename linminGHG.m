function [funmin,xmin]=linminGHG(func,dfunc,xtriplet,ftriplet,tol)
    
    xtriplet=xtriplet(:); ftriplet=ftriplet(:);  x=xtriplet(2) ; f=ftriplet(2);
    
    disp(' linmin ghg ')
    

    
    [xtriplet,ix] = sort(xtriplet) ;
    ftriplet=ftriplet(ix) ;

    a=xtriplet(1); b=xtriplet(2); c=xtriplet(3);
    fa=ftriplet(1); fb=ftriplet(2); fc=ftriplet(3);
    [ xnew ] = parabolamin(a,b,c,fa,fb,fc);
    
    % dx=dfunc(x);
    
    fnew=func(xnew);
   
  
      
  %  A=[A ; 0 1 0 ]; b=[b ; dx ];
 
  
        
    
    disp([' xtriplet : ',num2str(a),', ',num2str(b),' , ',num2str(c)])
    disp([' ftriplet : ',num2str(fa),', ',num2str(fb),', ',num2str(fc)])
    disp(['     xnew : ',num2str(xnew),', f new : ',num2str(fnew,10)])
  %  disp(['     dxnew : ',num2str(dxnew),', dx: ',num2str(dx)])
    
  
    
    
    
    [xttt,ix] = sort([xtriplet ; xnew]);
    fttt=[ftriplet ; fnew ];
    
    a=xttt(1); b=xttt(2); c=xttt(3);
    fa=fttt(ix(1)); fb=fttt(ix(2)); fc=fttt(ix(3));
   
    [ xnew ] = parabolamin(a,b,c,fa,fb,fc);
    fnew=func(xnew);
        
    disp([' xtriplet : ',num2str(a),', ',num2str(b),' , ',num2str(c)])
    disp([' ftriplet : ',num2str(fa),', ',num2str(fb),', ',num2str(fc)])
    disp(['     xnew : ',num2str(xnew),', f new : ',num2str(fnew,10)])
   % disp(['     dxnew : ',num2str(dxnew),', dx: ',num2str(dx)])
    
    funmin=fnew ; xmin=xnew ;
    
    
    
    
    
    return
    
    
    
    
    
    
    
    
    
    
end