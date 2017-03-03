function MUA=LoadAntarcticMUAexample

UaHomeDirectory=getenv('UaHomeDirectory');
cd(UaHomeDirectory)

cd UaUtilities\


load('MUA_Antarctica.mat','MUA')


end