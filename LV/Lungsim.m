
clear all
global dt npx f0 t0 tmax vtot tid xmesh konc dmol ru sc inito2 gb diff Model aa x ngs a cf cl lf ll rf rl;
 
dt = 0.01;
npx = 400;
f0 = 0.2;
t0 = 5;
tmax = 5;
vtot = 500;
dmol = 0.225;
ru = 0.153;
sc = ru/dmol;
inito2 = 0;
gb = 9;
diff = 3;
Model = 1;
cf = 1;
cl = 1;
lf = 1;
ll = 1;
rf = 1;
rl = 1;
double(dt);
double(npx);
double(t0);
double(tmax);
double(vtot);
double(dmol);
double(ru);
double(sc);
double(inito2);

Granssnitt 