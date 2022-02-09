%--------------------------------------------------------------------------------------
% MORPHOMETRICAL ROUTINE
% Routine to generate lung morphometricaldata as input for 
% the convection diffuion equation (According to Weibel).
%--------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------
% GL  = Length of a branch generation (cm)
% GA  = Total crosssection per branch-generation (cm^2)
% GV  = Total volume per branch-generation (cm^3)
% SAC = Volume of a SAC of a generation (cm^3)
% ALV = Volume of the alveoli of a generation (cm^3)
% TV  = Total volume of tube system 
%
% Vector index i represents (i-1) generation
%--------------------------------------------------------------------------------------
global npx xmesh a aa ngs;

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

%npx = 400; 
factot=sum(ALV);  
cc=TV/npx;
ngs=ceil(GV(1:24)/cc); 
dx=GL(1:24)./ngs(1:24);
npx=sum(ngs);  % New npx after adjustment

double(npx);
double(factot);
double(cc);
double(ngs);
double(dx);

jj = 0;

for i=1:24
   for j=1:ngs(i)   
      if j+jj>1,
         xmesh(jj+j)=xmesh(jj+j-1)+dx(i);
      else
         xmesh(jj+j)=dx(i);
      end   
      a(jj+j)=GA(i);  
      aa(jj+j)=GA(i)+(ALV(i)+SAC(i))/GL(i);
      ra(jj+j)=sqrt((aa(jj+j)/(pi*2^(i-1))));   
      fac(jj+j)=ALV(i)/(GV(i)+SAC(i)+ALV(i));
   end
   jj=jj+ngs(i);
end

xmesh=xmesh(1:npx);  % Get the right index 
a=a(1:npx);
aa=aa(1:npx);
ra=ra(1:npx);
fac=fac(1:npx);
Enkelmain;    % Start the mainprogram

