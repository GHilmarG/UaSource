function [xGL,yGL,nxGL,nyGL]=FindNiceLookingGLforPlottingPurposes(Method,DTxy,GF,CtrlVar,xBoundary,yBoundary)
    
    % finds a `nice looking' GL good for plotting purposes
    
    
    switch Method
        
        case 'Contouring'
            
            % maps GF.node onto a regular grid, finds the contours of GF on that gird
            % assumes the user is only interested in one GL and returns the longest contour line if many are found
            
            
            [GLx,GLy,GLxUpper,GLyUpper,GLxLower,GLyLower] = FindGL(DTxy,GF.node,CtrlVar,xBoundary,yBoundary);
                
            glx=GLx ; gly=GLy;
            % are the multiple lines? if so then take the longest one only
            if any(isnan(glx))
                temp=sort(find(isnan(glx)|isnan(gly)));
                temp=[0;temp;numel(glx)+1];
                [temp2,I]=max(diff(temp));
                
                n1=temp(I)+1 ;
                n2=temp(I+1)-1 ;
                
                xGL=glx(n1:n2); yGL=gly(n1:n2);
            else
                xGL=glx; yGL=gly;
            end
            
                        
            
        case 'ElementWise'

            
            CtrlVar.GLthreshold=0.1; % it's fine the have GLthreshold=0.5 in SSS2d, but for plotting purposes
            % it is better to define the grounding line close to where it is fully floating
            % i.e. with a value close to 0
            %%
            GF = GL2d(B,S,h,rhow,rho,connectivity,CtrlVar);
            [GLgeo,GLinfo,xGL,yGL]=GLgeometry(connectivity,coordinates,GF,CtrlVar);
            % pick out individual grounding lines (this must be improved by lining up the GL edges
            % assuming that if there are more than one GL they are split up by a NaN in (xGL,yGL)
            %%
            if any(isnan(xGL))
                ind=find(isnan(xGL));
                xGL=xGL(1:ind(1)-1); yGL=yGL(1:ind(1)-1);
            end
            
            [xGL,yGL] = Arrange2dPos(xGL,yGL);
            
        otherwise
            error('what case')
    end
    
    [xGL,yGL,nxGL,nyGL]=SplineLine(xGL,yGL,CtrlVar);  % create GL curve
    
end
