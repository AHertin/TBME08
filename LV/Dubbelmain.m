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
% The first npxc elements are for the tube before the branch point the next npxl elements
% are for the left branch and the npxr last elements are for the right branch
%--------------------------------------------------------------------------------------
global inito2 uold npxtot dxp dxm M d uold unew;

inito2=inito2/100;

uold = zeros(npxtot,1);  
unew = zeros(npxtot,1);
cold = ones(npxtot,1);
cold = inito2*cold;
cnew = zeros(npxtot,1);
d = zeros(npxtot,1);
rey = zeros(npxtot,1);
alpha = zeros(npxtot,1);

double(uold);
double(unew);
double(cold);
double(cnew);
double(d);
double(rey);
double(alpha);

global dt npx t0 tmax vtot tid xmesh konc dmol ru sc M upw dda diff;

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
d(1:npxtot)=dmol; 
alpha(1:npxtot)=ra(1:npxtot).*(sqrt(omega/ru));  

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
dda = zeros(npxtot-1,1);
dxp = zeros(npxtot-1,1);
dxm = zeros(npxtot-1,1);
upw = zeros(npxtot-1,1);
diag1 = zeros(npxtot,1);
diag2 = zeros(npxtot,1);
diag3 = zeros(npxtot,1);

cl_clm1 = 0;	
cl_cl = 0;
cl_lf = 0;
cl_rf = 0;

lf_cl = 0;
lf_lf = 0;
lf_lfp1 = 0;

rf_cl = 0;
rf_rf = 0;
rf_rfp1 = 0;

ll_llm1 = 0;
ll_ll = 0;
ll_llp1 = 0;


rs = zeros(npxtot,1);

double(dda);
double(dxp);
double(dxm);
double(upw);
double(diag1);
double(diag2);
double(diag3);

double(cl_clm1);	
double(cl_cl);
double(cl_lf);
double(cl_rf);

double(lf_cl);
double(lf_lf);
double(lf_lfp1);

double(rf_cl);
double(rf_rf);
double(rf_rfp1);

double(ll_llm1);
double(ll_ll);
double(ll_llp1);

double(rs);


for i=2:npxtot-1                    
  dda(i) = a(i+1) - a(i-1);
  dxp(i) = xmesh(i+1) - xmesh(i); 
  dxm(i) = xmesh(i) - xmesh(i-1);
end   

dda(cl) = (a(lf)+a(rf)) - a(cl-1);					
dda(lf) = a(lf+1) - a(cl)*(a(lf)/(a(lf)+a(rf)));	
dda(rf) = a(rf+1) - a(cl)*(a(rf)/(a(lf)+a(rf)));	
dda(ll) = a(ll) - a(ll-1);  
dda(1) = a(2) - a(1);
dxp(cl) = (xmesh(lf)+xmesh(rf))/2 - xmesh(cl);
dxp(ll) = xmesh(ll) - xmesh(ll-1);	
dxp(1) = xmesh(2) - xmesh(1);
dxm(lf) = xmesh(lf) - xmesh(cl);
dxm(rf) = xmesh(rf) - xmesh(cl);
dxm(1) = xmesh(1);



% konc anv för 3D-plott
konc=zeros(npxtot,npt/10);
tid=(1:npt)';

for i=1:npt
   time=i*dt;
   
%--------------------------------------------------------------------------------------
% Depending on the sign of the velocity an uwind, downwind
% or central difference scheme is chosen
%--------------------------------------------------------------------------------------
  
  upw(1:npxtot-1)=uold(1:npxtot-1)-d(1:npxtot-1).*dda(1:npxtot-1)./(aa(1:npxtot-1).* ...
     (dxp(1:npxtot-1)+dxm(1:npxtot-1)));
      
   uneg=find(upw(1:npxtot-1)<0);   	
   upos=find(upw(1:npxtot-1)>=0); 
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
   

%----------------------------------------------------------------------------------------
% The continuity conditions for the nodes at both sides of 
% the branch point 
%----------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------
% Equations for row cl in the matrix, the first equation if uppwind scheme else the second one
% Ex: cl_clm1 = row cl, column cl-1
%----------------------------------------------------------------------------------------
if upw(cl-1) >= 0
   cl_clm1 = -2*d(cl)*a(cl)*dt/(aa(cl)*(dxp(cl)+dxm(cl))*dxm(cl)) + ...
      d(cl)*dda(cl)*dt/(aa(cl)*(dxp(cl)+dxm(cl))*dxm(cl)) - uold(cl)*dt/dxm(cl);
else
   cl_clm1 = -2*d(cl)*a(cl)*dt/(aa(cl)*(dxp(cl)+dxm(cl))*dxm(cl)); 
end


if upw(cl) >= 0 
   cl_cl = 1 + 2*d(cl)*a(cl)*dt*(1/dxp(cl)+1/dxm(cl))/(aa(cl)*(dxp(cl)+dxm(cl))) - ...
      d(cl)*dda(cl)*dt/(aa(cl)*(dxp(cl)+dxm(cl))*dxm(cl)) + uold(cl)*dt/dxm(cl);
else 
   cl_cl = 1 + 2*d(cl)*a(cl)*dt*(1/dxp(cl)+1/dxm(cl))/(aa(cl)*(dxp(cl)+dxm(cl))) + ...
      d(cl)*dda(cl)*dt/(aa(cl)*(dxp(cl)+dxm(cl))*dxp(cl)) - uold(cl)*dt/dxp(cl);
end


if upw(lf) >= 0
   cl_lf = -2*d(cl)*a(cl)*a(lf)*dxm(lf)*dt/(aa(cl)*(dxp(cl)+dxm(cl))*dxp(cl)* ...
      (a(lf)*dxm(lf)+a(rf)*dxm(rf)));
else 
   cl_lf = -2*d(cl)*a(cl)*a(lf)*dxm(lf)*dt/(aa(cl)*(dxp(cl)+dxm(cl))*dxp(cl)* ...
      (a(lf)*dxm(lf)+a(rf)*dxm(rf))) - d(cl)*dda(cl)*a(lf)*dxm(cl)*dt/ ...
      (aa(cl)*(dxp(cl)+dxm(cl))*dxp(cl)*(a(lf)*dxm(lf)+a(rf)*dxm(rf))) + ...
   uold(cl)*a(lf)*dxm(lf)*dt/(dxp(cl)*(a(lf)*dxm(lf)+a(rf)*dxm(rf)));
end


if upw(rf) >=0
   cl_rf = -2*d(cl)*a(cl)*a(rf)*dxm(rf)*dt/(aa(cl)*(dxp(cl)+dxm(cl))*dxp(cl)* ...
      (a(lf)*dxm(lf)+a(rf)*dxm(rf)));
else
   cl_rf = -2*d(cl)*a(cl)*a(rf)*dxm(rf)*dt/(aa(cl)*(dxp(cl)+ ...
      dxm(cl))*dxp(cl)*(a(lf)*dxm(lf)+a(rf)*dxm(rf))) - ...
      d(cl)*dda(cl)*a(rf)*dxm(rf)*dt/(aa(cl)*(dxp(cl)+ ...
      dxm(cl))*dxp(cl)*(a(lf)*dxm(lf)+a(rf)*dxm(rf))) + ...
   uold(cl)*a(rf)*dxm(rf)*dt/(dxp(cl)*(a(lf)*dxm(lf)+a(rf)*dxm(rf)));
end

%----------------------------------------------------------------------------------------
% Equations for row lf in the matrix, the first equation if uppwind scheme else the second one
%----------------------------------------------------------------------------------------
if upw(cl) >= 0
   lf_cl = -2*d(lf)*a(lf)*dt/(aa(lf)*(dxp(lf)+dxm(lf))*dxm(lf)) + ...
      d(lf)*dda(lf)*dt/(aa(lf)*(dxp(lf)+dxm(lf))*dxm(lf)) - uold(lf)*dt/dxm(lf);
else
   lf_cl = -2*d(lf)*a(lf)*dt/(aa(lf)*(dxp(lf)+dxm(lf))*dxm(lf));
end


if upw(lf) >= 0
   lf_lf = 1 + 2*d(lf)*a(lf)*dt*(1/dxp(lf)+1/dxm(lf))/(aa(lf)*(dxp(lf)+dxm(lf))) - ...
      d(lf)*dda(lf)*dt/(aa(lf)*(dxp(lf)+dxm(lf))*dxm(lf)) + uold(lf)*dt/dxm(lf);
else
   lf_lf = 1 + 2*d(lf)*a(lf)*dt*(1/dxp(lf)+1/dxm(lf))/(aa(lf)*(dxp(lf)+dxm(lf))) + ...
      d(lf)*dda(lf)*dt/(aa(lf)*(dxp(lf)+dxm(lf))*dxp(lf)) - uold(lf)*dt/dxp(lf);
end


if upw(lf+1) >= 0
   lf_lfp1 = -2*d(lf)*a(lf)*dt/(aa(lf)*(dxp(lf)+dxm(lf))*dxp(lf));
else
   lf_lfp1 = -2*d(lf)*a(lf)*dt/(aa(lf)*(dxp(lf)+dxm(lf))*dxp(lf)) - ...
      d(lf)*dda(lf)*dt/(aa(lf)*(dxp(lf)+dxm(lf))*dxp(lf)) + ...
      uold(lf)*dt/dxp(lf);
end

%----------------------------------------------------------------------------------------
% Equations for row rf in the matrix, the first equation if uppwind scheme else the second one
%----------------------------------------------------------------------------------------
if upw(cl) >= 0
   rf_cl = -2*d(rf)*a(rf)*dt/(aa(rf)*(dxp(rf)+dxm(rf))*dxm(rf)) + ...
      d(rf)*dda(rf)*dt/(aa(rf)*(dxp(rf)+dxm(rf))*dxm(rf)) - uold(rf)*dt/dxm(rf);
else
   rf_cl = -2*d(rf)*a(rf)*dt/(aa(rf)*(dxp(rf)+dxm(rf))*dxm(rf));
end


if upw(rf) >= 0
   rf_rf = 1 + 2*d(rf)*a(rf)*dt*(1/dxp(rf)+1/dxm(rf))/(aa(rf)*(dxp(rf)+dxm(rf))) - ...
      d(rf)*dda(rf)*dt/(aa(rf)*(dxp(rf)+dxm(rf))*dxm(rf)) + ...
      uold(rf)*dt/dxm(rf);
else
   rf_rf = 1 + 2*d(rf)*a(rf)*dt*(1/dxp(rf)+1/dxm(rf))/(aa(rf)*(dxp(rf)+dxm(rf))) + ...
      d(rf)*dda(rf)*dt/(aa(rf)*(dxp(rf)+dxm(rf))*dxp(rf)) - uold(rf)*dt/dxp(rf);
end


if upw(rf+1) >= 0
   rf_rfp1 = -2*d(rf)*a(rf)*dt/(aa(rf)*(dxp(rf)+dxm(rf))*dxp(rf));
else
   rf_rfp1 = -2*d(rf)*a(rf)*dt/(aa(rf)*(dxp(rf)+dxm(rf))*dxp(rf)) - ...
      d(rf)*dda(rf)*dt/(aa(rf)*(dxp(rf)+dxm(rf))*dxp(rf)) + ...
      uold(rf)*dt/dxp(rf);
end


%--------------------------------------------------------------------------------------
% The boundary points at the alveoli in the right branch 
%--------------------------------------------------------------------------------------
   diag1(npxtot) = -2*d(npxtot)*dt*a(npxtot)/(aa(npxtot)*(xmesh(npxtot)-xmesh(npxtot-1))^2);
   diag2(npxtot) = 1 + 2*d(npxtot)*dt*a(npxtot)/(aa(npxtot)*(xmesh(npxtot)-xmesh(npxtot-1))^2);
   diag3(npxtot) = uold(npxtot)*dt/(2*(xmesh(npxtot)-xmesh(npxtot-1))) - ...
      d(npxtot)*dt*a(npxtot)/(aa(npxtot)*(xmesh(npxtot)-xmesh(npxtot-1))^2);
   
%--------------------------------------------------------------------------------------
% The boundary points at the alveoli in the left branch, row ll in the matrix
%--------------------------------------------------------------------------------------
   ll_llm1 = -2*d(ll)*dt*a(ll)/(aa(ll)*(xmesh(ll)-xmesh(ll-1))^2);
   ll_ll = 1 + 2*d(ll)*dt*a(ll)/(aa(ll)*(xmesh(ll)-xmesh(ll-1))^2);
   ll_llp1 = 0;

%--------------------------------------------------------------------------------------
% For the right hand side there can be chosen out of two possibility's
% or a combination of those two)
% 1) dc/dx = 0 at the alveoli and a distributed sink term
% 2) dc/dx = A and no source/sink term
%--------------------------------------------------------------------------------------
   rs(1:npxtot-1)=cold(1:npxtot-1)-qs*dt*fac(1:npxtot-1)/factot;
   
%--------------------------------------------------------------------------------------
% dc/dx = F
% rs(npx)=c(npx)-2*dx*F*diag3(npx)
% F=-J/A*D = -4/30028*0.225 = 0.000592  For Weibel's scaled model
%--------------------------------------------------------------------------------------
   if dcdx>0.000001 
      rs(npxtot)=cold(npxtot)+diag3(npxtot)*2*(xmesh(npxtot)-xmesh(npxtot-1))*dcdx/(30028*d(npxtot));
      rs(ll)=cold(ll)+diag3(ll)*2*(xmesh(ll)-xmesh(ll-1))*dcdx/(30028*d(ll));
   else
      rs(npxtot)=cold(npxtot)-qs*dt*fac(npxtot)/factot;
      rs(ll)=cold(ll)-qs*dt*fac(ll)/factot;
   end
   
   rs(1)=cold(1)-diag1(1)*0.2;

%--------------------------------------------------------------------------------------
% Solve the system of equations
%--------------------------------------------------------------------------------------
	M = zeros(npxtot,npxtot);
   double(M);
   
	M = diag(diag1(2:npxtot),-1) + diag(diag2(1:npxtot)) + diag(diag3(1:npxtot-1),1);
 
%--------------------------------------------------------------------------------------
% Put the values from the continuity conditions for the m-node in matrix M
%--------------------------------------------------------------------------------------
   M(cl,cl-1) = cl_clm1;
   M(cl,cl) = cl_cl;
   M(cl,lf) = cl_lf;
   M(cl,rf) = cl_rf;
   
%--------------------------------------------------------------------------------------
% Put the values from the continuity conditions for the l-node in matrix M
%--------------------------------------------------------------------------------------
   M(lf,cl) = lf_cl;
   M(lf,lf) = lf_lf;
   M(lf,lf+1) = lf_lfp1;
   
%--------------------------------------------------------------------------------------
% Put the values from the continuity conditions for the r-node in matrix M
%--------------------------------------------------------------------------------------
   M(rf,cl) = rf_cl;
   M(rf,ll) = 0;
   M(rf,rf) = rf_rf;
   M(rf,rf+1) = rf_rfp1;
   
%--------------------------------------------------------------------------------------
% Put the values from the boundery conditions for the ll-node in matrix M
%--------------------------------------------------------------------------------------
   M(ll,ll-1) = ll_llm1;
   M(ll,ll) = ll_ll;
   M(ll,ll+1) = ll_llp1;
   
   cnew = M\rs;
   %cnew = inv(M)*rs;
   
%--------------------------------------------------------------------------------------
% Calculate the new velocity in the tube system
% Tidal volume (= 500 cm^3) = u*aa(1)=100*pi*sin(time/time2)
% Overwrite the old solution with the new one to advance to the new
% time step
%--------------------------------------------------------------------------------------
   unew(1)=(vtot*pi/t0)*sin(2*pi*time/t0)/aa(1);
   
   for j=2:cl
      unew(j)=unew(j-1)*aa(j-1)/aa(j);
   end   
      
   unew(lf)=(diff/(diff+1))*unew(cl)*aa(cl)/aa(lf);	
   
   for j=lf+1:ll
      unew(j)=unew(j-1)*aa(j-1)/aa(j);
   end   
   
   unew(rf)=(1/(diff+1))*unew(cl)*aa(cl)/aa(rf);     			
   
   for j=rf+1:rl
      unew(j)=unew(j-1)*aa(j-1)/aa(j);
   end   
   
   
   rey(2:npxtot)=abs(2*unew(2:npxtot).*ra(2:npxtot)/ru);
   
   
   double(rey);
   double(reymax);
   
   if (reymax<=max(rey))
      reymax=max(rey);
   end
   
%--------------------------------------------------------------------------------------
% Alpha greater then five with different in- and expiration fomulas
% for deff at turbulent flow
%--------------------------------------------------------------------------------------
%   for j=2:npxtot
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
%  		if ((alpha(j)>1) & (alpha(j)<=5) & (rey(j)>2000) & (unew(j)>=0))
%     	 	d(j)=dmol*(1+0.37*rey(j)*sc);
%     	end
%  
%    	if ((alpha(j)>1) & (alpha(j)<=5) & (rey(j)>2000) & (unew(j)<0))
%      	d(j)=dmol*(1+1.08*rey(j)*sc);
%   	end
% 
%  		if ((alpha(j)>1) & (alpha(j)<=5) & (rey(j)>=6) & (rey(j)<=2000))
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
%      end
%   end
   
   d(1)=d(2);
   rey(1)=rey(2);
   
   if defff==0
      d(1:npxtot)=dmol;
   else
      d(1:npxtot)=d(1:npxtot)*defff;   
   end
   
   uold(1:npxtot)=unew(1:npxtot);
   cold(1:npxtot)=cnew(1:npxtot);

%--------------------------------------------------------------------------------------
% Save the koncentration in a matrix
%--------------------------------------------------------------------------------------
   if mod(time,0.1)==0
      konc(:,i/10)=cnew(1:npxtot);       
   end
   
      
end

konc=100*konc;

global x;

xmax = max(xmesh);             
xmax = xmax * ones(npxtot,1);
x = xmax - xmesh;


Dubbelplott;







