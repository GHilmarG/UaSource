load Ex1aGLpos GLpos GLtime GLstd


figure 
plot(GLtime,GLpos/1000,'-o')
hold on
plot(GLtime,(GLpos-GLstd)/1000,'r')
plot(GLtime,(GLpos+GLstd)/1000,'r')


xlabel('t (a)')
ylabel('GL position (km)')



