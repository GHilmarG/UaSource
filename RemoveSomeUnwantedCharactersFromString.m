function OutputString=RemoveSomeUnwantedCharactersFromString(InputString)

OutputString=InputString;

OutputString=replace(OutputString,"---","-");
OutputString=replace(OutputString,"--","-");
OutputString=replace(OutputString,"\-","\");
OutputString=replace(OutputString,"/-","/");
OutputString=replace(OutputString,"2.5","2k5");







end