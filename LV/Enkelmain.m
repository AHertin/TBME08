%--------------------------------------------------------------------------------------

% Program that compute oxygen koncentration in the lungs.

% The geometry is building on Weibels scaled model.

%-------------------------------------------------------------------------------------- 



%--------------------------------------------------------------------------------------

% MAIN PROGRAM

% The Program solves the one dimensional convection diffusion 

% equation in a channel with variable crosssectional area by 

% means of an implicit, forward in time, upwind, central or 

% downwind discretisation, unequal distance step and a 

% sink/source term and/or dc/dx

%--------------------------------------------------------------------------------------



%--------------------------------------------------------------------------------------

% DT    = Time step (s)
% DMOL  = Molecular diffusion coefficeint (cm^2/s)
% D(K)  = Augmented diffusion coefficient at node k
% TMAX  = Time till which will be calculated (s) 
% T0    = Length of breathing cycle (s)
% TIME  = Current time (s)
% NPT   = Total number of time step
% RU    = Kinematic viscosity (cm^2/s)
% VTOT  = Tidal volume (ml)
% NPT   = Number of time steps 
% QS    = Q strength (cm^3/s)
% DCDX  = Strength (cm^3/s)
% VTOT  = Tidal volume (ml)
% DEFFF = 1 / 0 Variable or konstant diffusion koefficient
% UOLD  = Velocity
% UNEW  = Velocity
% COLD  = Koncentration
% CNEW  = Koncentration
% REY   = Reynolds number

%--------------------------------------------------------------------------------------



%--------------------------------------------------------------------------------------

% Parameters, constants and variables for the main program

%--------------------------------------------------------------------------------------

global inito2; 

inito2=inito2/100;

uold = zeros(npx,1);  

unew = zeros(npx,1);
cold = ones(npx,1);
cold = inito2*cold;

cnew = zeros(npx,1);
d = zeros(npx,1);
rey = zeros(npx,1);
alpha = zeros(npx,1);


double(uold);
double(unew);
double(cold);
double(cnew);
double(d);
double(rey);
double(alpha);

global dt npx t0 tmax vtot tid xmesh konc dmol ru sc;

%dt=0.01;            
%dmol=0.225;
%dmol=0;    							

ru=0.153;

%sc=ru/dmol;		
%sc=0;    							
%t0=5;               
%tmax=5;             
omega=2*pi/t0;
%vtot=500;           

npt=ceil(tmax/dt);
defff=0;								         
%qs=4;
qs=0;
dcdx=0;             
d(1:npx)=dmol;
alpha(1:npx)=ra(1:npx).*(sqrt(omega/ru));  


double(omega);
double(npt);
double(defff);
double(qs);


%--------------------------------------------------------------------------------------

% Set the initial condition

%--------------------------------------------------------------------------------------

reymax=0;  
time=0;

%--------------------------------------------------------------------------------------

% Start of the time loop

%--------------------------------------------------------------------------------------


%--------------------------------------------------------------------------------------

% DDA	    = If k is the current node dda(k)=a(k+1)-a(k-1)
% DXP     = If k is the current node dxp(k)=xmesh(k+1)-xmesh(k)
% DXM     = If k is the current node dxm(k)=xmesh(k)-xmesh(k-1)
% UPW     = Velocity
% DIAG1   = Diagonal 1 in the  M matrix
% DIAG2   = Diagonal 2 in the  M matrix
% DIAG3   = Diagonal 3 in the  M matrix
% UPWNEG  = Vectorindex for upw(k)<0     
% UPWZERO = Vecorindex for upw(k)=0
% UPWPOS  = Vecorindex for upw(k)>0
% RS      = Right hand side for the system of equations
% M       = The three diagonal matrix 

%--------------------------------------------------------------------------------------

dda = zeros(npx-1,1);
dxp = zeros(npx-1,1);
dxm = zeros(npx-1,1);
upw = zeros(npx-1,1);
diag1 = zeros(npx,1);
diag2 = zeros(npx,1);
diag3 = zeros(npx,1);
rs = zeros(npx,1);


double(dda);
double(dxp);
double(dxm);
double(upw);
double(diag1);
double(diag2);
double(diag3);
double(rs);


for i=2:npx-1                    
  dda(i) = a(i+1) - a(i-1);
  dxp(i) = xmesh(i+1) - xmesh(i); 
  dxm(i) = xmesh(i) - xmesh(i-1);
end   

dda(1) = a(2) - a(1);
dxp(1) = xmesh(2) - xmesh(1); 
dxm(1) = xmesh(1);


% konc anv f?r 3D-plott

konc=zeros(npx,npt/10); 
tid=(1:npt)';


for i=1:npt
   time=i*dt;
   
%--------------------------------------------------------------------------------------

% Depending on the sign of the velocity an uwind, downwind
% or central difference scheme is chosen

%--------------------------------------------------------------------------------------

  upw(1:npx-1)=uold(1:npx-1)-d(1:npx-1).*dda(1:npx-1)./(aa(1:npx-1).*(dxp(1:npx-1)+...
     dxm(1:npx-1)));

   uneg=find(upw(1:npx-1)<0);   

   upos=find(upw(1:npx-1)>=0); 
   %uzero=find(upw(1:npx-1)==0);

   %--------------------------------------------------------------------------------------

% U < 0

%--------------------------------------------------------------------------------------

   diag1(uneg) = -d(uneg)*2*dt.*a(uneg)./((dxp(uneg)+dxm(uneg)).* ...
      dxm(uneg).*aa(uneg));
   

   diag2(uneg) = 1 + 2*d(uneg)*dt.*a(uneg).*(1./dxp(uneg)+1./dxm(uneg))./ ...
      (aa(uneg).*(dxp(uneg)+dxm(uneg))) - uold(uneg)*dt./dxp(uneg) + ...
      dda(uneg).*d(uneg)*dt./(aa(uneg).*(dxp(uneg)+dxm(uneg)).*dxp(uneg));

   
   diag3(uneg) = uold(uneg)*dt./dxp(uneg) - 2*d(uneg)*dt.*a(uneg)./ ...
      ((dxp(uneg)+dxm(uneg)).*dxp(uneg).*aa(uneg)) - ... 
      dda(uneg).*d(uneg)*dt./(aa(uneg).*(dxp(uneg)+dxm(uneg)).*dxp(uneg));
   
    
%--------------------------------------------------------------------------------------

% U > or = 0

%--------------------------------------------------------------------------------------

   diag1(upos) = -uold(upos)*dt./dxm(upos) - 2*d(upos)*dt.*a(upos)./ ...
      ((dxp(upos)+dxm(upos)).*dxm(upos).*aa(upos)) + ...
      dda(upos).*d(upos)*dt./(aa(upos).*(dxp(upos)+dxm(upos)).*dxm(upos));

   diag2(upos) = 1 + 2*d(upos)*dt.*a(upos).*(1./dxp(upos)+1./dxm(upos))./ ...
      ((dxp(upos)+dxm(upos)).*aa(upos)) + uold(upos)*dt./dxm(upos) - ...
      dda(upos).*d(upos)*dt./(aa(upos).*(dxp(upos)+dxm(upos)).*dxm(upos));


   diag3(upos) = -2*d(upos)*dt.*a(upos)./((dxp(upos)+dxm(upos)).* ...
      dxp(upos).*aa(upos)); 

%--------------------------------------------------------------------------------------
% U = 0
%--------------------------------------------------------------------------------------
%   diag1(uzero)=-uold(uzero)*dt./(dxp(uzero)+dxm(uzero))- ...
%      2*d(uzero)*dt.*a(uzero)./((dxp(uzero)+dxm(uzero)).* ...
%      dxm(uzero).*aa(uzero))+dda(uzero).*d(uzero)*dt./ ...
%      (aa(uzero).*(dxp(uzero)+dxm(uzero)).^2);
%   
%   diag2(uzero)=1+2*dt*d(uzero).*a(uzero).*(1./dxp(uzero)+1./ ...
%      dxm(uzero))./((dxp(uzero)+dxm(uzero)).*aa(uzero));
%   
%   diag3(uzero)=uold(uzero)*dt./(dxp(uzero)+dxm(uzero))- ...
%      2*dt*d(uzero).*a(uzero)./((dxp(uzero)+dxm(uzero)).* ...
%      dxp(uzero).*aa(uzero))-dda(uzero).*d(uzero)*dt./ ...
%      (aa(uzero).*(dxp(uzero)+dxm(uzero)).^2);
   

%--------------------------------------------------------------------------------------

% The boundary points at the alveoli (dxp = npx)

%--------------------------------------------------------------------------------------

   diag1(npx) = -2*d(npx)*dt*a(npx)/(aa(npx)*(xmesh(npx)-xmesh(npx-1))^2);

   diag2(npx) = 1 + 2*d(npx)*dt*a(npx)/(aa(npx)*(xmesh(npx)-xmesh(npx-1))^2);

   diag3(npx) = uold(npx)*dt/(2*(xmesh(npx)-xmesh(npx-1))) - ...
      d(npx)*dt*a(npx)/(aa(npx)*(xmesh(npx)-xmesh(npx-1))^2);
     

%--------------------------------------------------------------------------------------

% For the right hand side there can be chosen out of two possibility's

% or a combination of those two)

% 1) dc/dx = 0 at the alveoli and a distributed sink term
% 2) dc/dx = A and no source/sink term

%--------------------------------------------------------------------------------------

   rs(1:npx-1)=cold(1:npx-1)-qs*dt*fac(1:npx-1)/factot;

%--------------------------------------------------------------------------------------

% dc/dx = F

% rs(npx)=c(npx)-2*dx*F*diag3(npx)

% F=-J/A*D = -4/30028*0.225 = 0.000592  For Weibel's scaled model

%--------------------------------------------------------------------------------------

   if dcdx>0.000001 
      rs(npx)=cold(npx)+diag3(npx)*2*(xmesh(npx)-xmesh(npx-1))*dcdx/(30028*d(npx));
   else
      rs(npx)=cold(npx)-qs*dt*fac(npx)/factot;
   end


   rs(1)=cold(1)-diag1(1)*0.2;	
   
%--------------------------------------------------------------------------------------

% Solve the system of equations

%--------------------------------------------------------------------------------------

    cnew = Tridiag(diag1, diag2, diag3, rs, npx);

%--------------------------------------------------------------------------------------

% Calculate the new velocity in the tube system
% Tidal volume (= 500 cm^3) = u*aa(1)=100*pi*sin(time/time2)
% Overwrite the old solution with the new one to advance to the new

% time step

%--------------------------------------------------------------------------------------

   unew(1)=(vtot*pi/t0)*sin(2*pi*time/t0)/aa(1);
   
   for j=2:npx

      unew(j)=unew(j-1)*aa(j-1)/aa(j);

   end   


   rey(2:npx)=abs(2*unew(2:npx).*ra(2:npx)/ru);
   
   double(rey);
   double(reymax);
   
   
   if (reymax<=max(rey))
      reymax=max(rey);

   end

   

%--------------------------------------------------------------------------------------

% Alpha greater then five with different in- and expiration fomulas

% for deff at turbulent flow

%--------------------------------------------------------------------------------------
%   for j=2:npx

%      if ((alpha(j)>5) & (rey(j)>(400*alpha(j))) & ((unew(j)>=0)))
%         d(j)=dmol*(1+0.37*rey(j)*sc);
%      end

%      

%      if ((alpha(j)>5) & (rey(j)>(400*alpha(j))) & ((unew(j)<0)))
%         d(j)=dmol*(1+1.08*rey(j)*sc);
%      end

%   
%      if ((alpha(j)>5) & (rey(j)>=6) & (rey(j)<(400*alpha(j))))
%         d(j)=dmol*(1+(0.666*((rey(j)*sc)^2.2))/((alpha(j)*sc^0.5)^2.5));
%      end

%      
%      if ((alpha(j)>5) & (rey(j)<6))
%         d(j)=dmol*(1+0.21*((rey(j)*sc)^2));
%      end

%      

%--------------------------------------------------------------------------------------

% 1 < alpha < 5 with different in- and expiration formulas for deff

% at turbulent flow

%--------------------------------------------------------------------------------------

%  	if ((alpha(j)>1) & (alpha(j)<=5) & (rey(j)>2000) & (unew(j)>=0))

%     	 d(j)=dmol*(1+0.37*rey(j)*sc);

%     end

%   

%    	if ((alpha(j)>1) & (alpha(j)<=5) & (rey(j)>2000) & (unew(j)<0))

%      	d(j)=dmol*(1+1.08*rey(j)*sc);

%   	end

%     

%   	if ((alpha(j)>1) & (alpha(j)<=5) & (rey(j)>=6) & (rey(j)<=2000))

%         d(j)=dmol*(1+(0.36*((rey(j)*sc)^1.1))/((alpha(j)*sc^0.5)^0.43));

%   	end

%   

%   	if ((alpha(j)>1) & (alpha(j)<=5) & (rey(j)<6))

%      	d(j)=dmol*(1+0.21*((rey(j)*sc)^2));

%   	end

%   

%   	if (alpha(j)<=1)

%      	d(j)=dmol;                         

%      	%d(j)=dmol*(1+0.21*(rey(j)*sc)^2);

%     end

%   end



   d(1)=d(2);

   rey(1)=rey(2);

   

   if defff==0

      d(1:npx)=dmol; 

   else

      d(1:npx)=d(1:npx)*defff;   

   end

   

   uold(1:npx)=unew(1:npx);

   cold(1:npx)=cnew(1:npx);



%--------------------------------------------------------------------------------------

% Save the koncentration in a matrix 

%--------------------------------------------------------------------------------------
   if mod(time,0.1)==0

      konc(:,i/10)=cnew(1:npx);      
   end
      
end


konc=100*konc;

global x;

xmax = max(xmesh);
xmax = xmax * ones(npx,1);
x = xmax - xmesh;

Enkelplott;


