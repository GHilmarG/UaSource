

function RenameFileRecursively(OldName,NewName)
    
    %%
    %
    %  RenameFileRecursivly(OldName,NewName)
    %
    % Change files names in current folder and all subfolders.
    %
    % Example:
    % To change all filenames from UaOutputs.m to DefineOutput.m in current folder and all
    % subfolders: 
    %
    %   RenameFileRecursively("UaOutputs.m","DefineOutputs.m")
    %
    %   RenameFileRecursively("Ua2D_InitialUserInput.m","DefineInitialInputs.m")
    %
    % Note: You might also want to change the function name within the m-files.
    % On unix you could do something like:
    %
    %   find . -name "DefineOutputs.m" -exec sed -i 's/UaOutputs/DefineOutputs/g' {} +
    %
    %
    %  find . -name "DefineInitialInputs.m" -exec sed -i 's/Ua2d_InitialUserInput/DefineInitialInputs/g' {} +
    %
    %
    % 
    %%

    Curdir=pwd;
    Files = dir("**/"+OldName);
    for I=1:numel(Files)
        cd(Files(I).folder)
        fprintf("Renaming %s to %s in Folder %s\n"',OldName,NewName,pwd) ;
        movefile(OldName,NewName)
    end
    cd(Curdir)
end