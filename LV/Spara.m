function spara

global x konc cf cl lf ll rf rl dt npx f0 t0 tmax vtot dmol ru sc inito2 gb diff Model;

[newmatfile, newpath] = uiputfile('*.mat', 'Spara'); 

save(strcat(newpath, newmatfile),'x','konc','cf','cl','lf','ll','rf','rl','dt','npx','f0','t0','tmax','vtot','dmol','ru','sc','inito2','gb','diff','Model'); 
 