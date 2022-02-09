function ladda

global x konc cf cl lf ll rf rl dt npx f0 t0 tmax vtot dmol ru sc inito2 gb diff Model;

[oldmatfile, oldpath] = uigetfile('*.mat', 'Ladda'); 

load(strcat(oldpath, oldmatfile));   