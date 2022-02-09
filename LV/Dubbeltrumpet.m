%--------------------------------------------------------------------------------------
% MORPHOMETRICAL ROUTINE
% Routine to generate lung morphometricaldata as input for 
% the convection diffuion equation (According to Weibel).
%--------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------
% GL  = Length of a branch generation (cm)
% GA  = Total crossection per branch-generation (cm^2)
% GV  = Total volume per branch-generation (cm^3)
% SAC = Volume of a SAC of a generation (cm^3)
% ALV = Volume of the alveoli of a generation (cm^3)
% TV  = Total volume of tube system 
%
% Vector index i represents (i-1) generation
%--------------------------------------------------------------------------------------
global npxtot npx xmesh a aa cf cl lf ll rf rl gb;

GL=zeros(24,1);
GA=zeros(24,1);
SAC=zeros(24,1);
ALV=zeros(24,1);

GL(1)=12;
GL(2)=4.76;
GL(3)=1.9;
GL(4)=0.65;
GL(5)=1.085;
GL(6)=0.915;
GL(7)=0.769;
GL(8)=0.65;
GL(9)=0.547;
GL(10)=0.462;
GL(11)=0.393;
GL(12)=0.333;
GL(13)=0.282;
GL(14)=0.231;
GL(15)=0.197;
GL(16)=0.171;
GL(17)=0.141;
GL(18)=0.121;
GL(19)=0.100;
GL(20)=0.085;
GL(21)=0.071;
GL(22)=0.060;
GL(23)=0.050;
GL(24)=0.043;

GA(1)=2.54;
GA(2)=2.34;
GA(3)=2.16;
GA(4)=1.44;
GA(5)=1.86;
GA(6)=2.25;
GA(7)=2.88;
GA(8)=3.88;
GA(9)=5.08;
GA(10)=6.97;
GA(11)=9.9;
GA(12)=14;
GA(13)=21.2;
GA(14)=31.6;
GA(15)=51.5;
GA(16)=82;
GA(17)=135;
GA(18)=219;
GA(19)=376;
GA(20)=664;
GA(21)=1218;
GA(22)=2225;
GA(23)=4045;
GA(24)=8091;

%--------------------------------------------------------------------------------------
% SAC=0 in the Weibel model
%--------------------------------------------------------------------------------------

ALV(1)=0;
ALV(2)=0;
ALV(3)=0;
ALV(4)=0;
ALV(5)=0;
ALV(6)=0;
ALV(7)=0;
ALV(8)=0;
ALV(9)=0;
ALV(10)=0;
ALV(11)=0;
ALV(12)=0;
ALV(13)=0;
ALV(14)=0;
ALV(15)=0;
ALV(16)=0;
ALV(17)=0;
ALV(18)=3.9;
ALV(19)=13.1;
ALV(20)=39.4;
ALV(21)=138;
ALV(22)=273;
ALV(23)=551;
ALV(24)=938;

GV = GL(1:24).*GA(1:24);
TV = sum(GV);

double(GL);
double(GA);
double(SAC);
double(ALV);
double(GV);
double(TV);


dx = zeros(24,1);
xmesh = zeros(4500,1);
a = zeros(4500,1);
aa = zeros(4500,1);
ra = zeros(4500,1);
fac = zeros(4500,1);

double(dx);
double(xmesh);
double(a);
double(aa);
double(ra);
double(fac);

factot = sum(ALV);
cc = TV/npx;
ngs = ceil(GV(1:24)/cc);
dx = GL(1:24)./ngs(1:24);
npx = sum(ngs);

double(npx);
double(factot);
double(cc);
double(ngs);
double(dx);

%----------------------------------------------------------------------------------
% gb  = The generation were the branch point is.
%       The branch point is in the beginning of that generation
% ngc = Number of generations before the branch point
% ngl = Number of generations in the left branch
% ngr = Number of generations in the right branch
% GLC = Length of a branch generation before the branch point (cm)
% GLL = Length of a branch generation in the left branch (cm)
% GLR = Length of a branch generation in the right branch (cm)
% GAC = Crossection in the tube before the branch per generation (cm^2)
% GAL = Crossection in the left branch per generation (cm^2)
% GAR = Crossection in the right branch per generation (cm^2)
%----------------------------------------------------------------------------------
%gb = 9;
ngc = gb;
ngl = 24-gb+1;
ngr = 24-gb+1;		
GLC(1:ngc) = GL(1:ngc);
GLL(1:ngl) = GL(gb:24);
GLR(1:ngr) = GL(gb:24);	
GAC(1:ngc) = GA(1:ngc);
GAL(1:ngl) = GA(gb:24)/2;
GAR(1:ngr) = GA(gb:24)/2;	
GVC(1:ngc) = GLC(1:ngc).*GAC(1:ngc);
GVL(1:ngl) = GLL(1:ngl).*GAL(1:ngl);
GVR(1:ngr) = GLR(1:ngr).*GAR(1:ngr);
SACC(1:ngc) = zeros(ngc,1);
SACL(1:ngl) = zeros(ngl,1);
SACR(1:ngr) = zeros(ngr,1);
ALVC(1:ngc) = ALV(1:ngc);
ALVL(1:ngl) = ALV(gb:24)/2;
ALVR(1:ngr) = ALV(gb:24)/2;  


%----------------------------------------------------------------------------------
% NPX   = Number of desired points in x direction (=axial tube)
% 
% CC    = Unit volume based on total volume of the tubes from
%         the generation considered
%
% DX    = X-mesh length step at a certain generation
%
% A     = Summed cross sectional area at a node point based
%         on tubes
%
% AA    = As A, but sacs and alveolar volume divide by tube
%	       length to achieve their contribution at the cross 
%         sectional area
%
% RA    = Radius per tube with alveoli at a certain generation
%
% FAC   = Factor, which represents the alveolar contribution
%         at a node point to the total alveolar volume 
%         (FACTOT), and is needed to generate the appropriate
%         sink strength in the convection diffusion equation
% 
% XMESH = The length to node 0
%
% NGS	  = Number of steps for generation 0 to 23   
%--------------------------------------------------------------------------------------
dxc = zeros(ngc,1); 
dxl = zeros(ngl,1); 
dxr = zeros(ngr,1); 
xmeshc = zeros(4500,1);
xmeshl = zeros(4500,1);
xmeshr = zeros(4500,1);
ac = zeros(4500,1);
aac = zeros(4500,1);
rac = zeros(4500,1);
facc = zeros(4500,1);
al = zeros(4500,1);
aal = zeros(4500,1);
ral = zeros(4500,1);
facl = zeros(4500,1);
ar = zeros(4500,1);
aar = zeros(4500,1);
rar = zeros(4500,1);
facr = zeros(4500,1);


double(dxc);
double(dxl);
double(dxr);
double(xmeshc);
double(xmeshl);
double(xmeshr);
double(ac);
double(al);
double(ar);
double(aac);
double(aal);
double(aar);
double(rac);
double(ral);
double(rar);
double(facc);
double(facl);
double(facr);

%npx = 400;
ngsc(1:gb) = ngs(1:gb);
ngsl(1:ngl) = ngs(gb:24);
ngsr(1:ngr) = ngs(gb:24);

ngsc(gb) = ceil(ngs(gb)/2);
ngsl(1) = floor(ngs(gb)/2);
ngsr(1) = floor(ngs(gb)/2);

if ngsc(gb) < 15
   ngsc(gb) = 15;
end

if ngsl(1) < 15
   ngsl(1) = 15;
end

if ngsr(1) < 15
   ngsr(1) = 15;
end

ngsl(2:ngl) = ceil(ngsl(2:ngl)/2);   
ngsr(2:ngr) = ceil(ngsr(2:ngr)/2); 

dxc=GLC(1:ngc)./ngsc(1:ngc);
dxl=GLL(1:ngl)./ngsl(1:ngl);
dxr=GLR(1:ngr)./ngsr(1:ngr);
npxc=sum(ngsc);
npxl=sum(ngsl); 
npxr=sum(ngsr);
npxtot=npxc+npxl+npxr;

double(npxc);
double(npxl);
double(npxr);
double(factot);
double(cc);
double(ngsc);
double(ngsl);
double(ngsr);
double(dxc);
double(dxl);
double(dxr);

jj = 0;

for i=1:ngc
   for j=1:ngsc(i)   
      if j+jj>1,
         xmeshc(jj+j)=xmeshc(jj+j-1)+dxc(i);
      else
         xmeshc(jj+j)=dxc(i);
      end   
      ac(jj+j)=GAC(i);  
      aac(jj+j)=GAC(i)+(ALVC(i)+SACC(i))/GLC(i);
      rac(jj+j)=sqrt((aac(jj+j)/(pi*2^(i-1))));   
      facc(jj+j)=ALVC(i)/(GVC(i)+SACC(i)+ALVC(i));
   end
   jj=jj+ngsc(i);
end

xmeshc=xmeshc(1:npxc);  % Get the right index 
ac=ac(1:npxc);
aac=aac(1:npxc);
rac=rac(1:npxc);
facc=facc(1:npxc);


jj = 0;

for i=1:ngl
   for j=1:ngsl(i)   
      if j+jj>1,
         xmeshl(jj+j)=xmeshl(jj+j-1)+dxl(i);
      else
         xmeshl(jj+j)=xmeshc(npxc)+dxl(i);
      end   
      al(jj+j)=GAL(i);  
      aal(jj+j)=GAL(i)+(ALVL(i)+SACL(i))/GLL(i);
      ral(jj+j)=sqrt((aal(jj+j)/(pi*2^(i-1))));   
      facl(jj+j)=ALVL(i)/(GVL(i)+SACL(i)+ALVL(i));
   end
   jj=jj+ngsl(i);
end

xmeshl=xmeshl(1:npxl);  % Get the right index 
al=al(1:npxl);
aal=aal(1:npxl);
ral=ral(1:npxl);
facl=facl(1:npxl);


jj = 0;

for i=1:ngr
   for j=1:ngsr(i)   
      if j+jj>1,
         xmeshr(jj+j)=xmeshr(jj+j-1)+dxr(i);
      else
         xmeshr(jj+j)=xmeshc(npxc)+dxr(i);
      end   
      ar(jj+j)=GAR(i);  
      aar(jj+j)=GAR(i)+(ALVR(i)+SACR(i))/GLR(i);
      rar(jj+j)=sqrt((aar(jj+j)/(pi*2^(i-1))));   
      facr(jj+j)=ALVR(i)/(GVR(i)+SACR(i)+ALVR(i));
   end
   jj=jj+ngsr(i);
end

xmeshr=xmeshr(1:npxr);  % Get the right index 
ar=ar(1:npxr);
aar=aar(1:npxr);
rar=rar(1:npxr);
facr=facr(1:npxr);

xmesh=[xmeshc; xmeshl; xmeshr]; 
a=[ac; al; ar];
aa=[aac; aal; aar];
ra=[rac; ral; rar];
fac=[facc; facl; facr];
cf=1;					
cl=npxc;
lf=npxc+1;				% Data for the first node in the left branch is in element fl 
ll=npxc+npxl;			% Data for the last node in the left branch is in element ll
rf=npxc+npxl+1;
rl=npxc+npxl+npxr;
vl=sum(GLL.*GAL); 	% The total volume of the left branch
vr=sum(GLR.*GAR);		% The total volume of the right branch
vb=vl+vr;				% The total volume of the two branches


Dubbelmain;          % Start the mainprogram

