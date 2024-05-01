function OutputString=RemoveSomeUnwantedCharactersFromString(InputString)

OutputString=InputString;

OutputString=replace(OutputString,"---","-");
OutputString=replace(OutputString,"--","-");
OutputString=replace(OutputString,"\-","\");
OutputString=replace(OutputString,"/-","/");
 





end