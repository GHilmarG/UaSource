function TestIfInputFileInWorkingDirectory(InputFile)



if ~exist(fullfile(cd,InputFile),'file')
    
    fprintf("-------------  The input-file %s not found in the working directory (%s).\n",InputFile,pwd)
    fprintf("-------------- This input-file is required.\n")
    error("Ua2D:InputFileNotFound","Required input file not found in working directory.")
    
end


end