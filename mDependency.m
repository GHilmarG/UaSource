
cd('C:\cygwin64\home\Hilmar\ghg\Ua\Source\mFiles')
[fList,pList] = matlab.codetools.requiredFilesAndProducts('Ua.m');

%%
for I=1:numel(fList)
    fprintf('%s \n',fList{I})
end


%%
%for I=1:numel(fList)
%    delete(fList{I})  % end of the world as we know it
%end