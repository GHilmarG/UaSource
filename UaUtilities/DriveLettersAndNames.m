function [driveletters,drivenames]=DriveLettersAndNames



r=java.io.File.listRoots;
for n=1:length(r)
     driveletters{n} = char(r(n).toString);
     drivenames{n}=DriveLetter2DriveName(driveletters{n});
end

end