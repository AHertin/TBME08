function CbEnkelplott(action)

global x konc t1 t2 t3 t4 t5 t6 t7 t8 x1 x2 x3 x4 proc1 proc2 proc3 proc4 tmax;


switch(action)
   
case 'tid1'							
   Str = get(gcbo,'String');				
   t1temp = ConvStr(Str);
   
   if(t1temp >= 0 & t1temp <= tmax)
      t1 = t1temp;
   else
   	t1Ptr = findobj(gcbf,'Tag','1Txt');
      set(t1Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(t1,4));	
   set(gcbo,'UserData',t1);
   
  	Str = get(gcbo,'String');				
     t1 = ConvStr(Str);
     
     
case 'tid2'							
   Str = get(gcbo,'String');				
   t2temp = ConvStr(Str);
   
   if(t2temp >= 0 & t2temp <= tmax)
      t2 = t2temp;
   else
   	t2Ptr = findobj(gcbf,'Tag','2Txt');
      set(t2Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(t2,4));	
   set(gcbo,'UserData',t2);
   
  	Str = get(gcbo,'String');				
   t2 = ConvStr(Str);
     
     
case 'tid3'							
   Str = get(gcbo,'String');				
   t3temp = ConvStr(Str);
   
   if(t3temp >= 0 & t3temp <= tmax)
      t3 = t3temp;
   else
   	t3Ptr = findobj(gcbf,'Tag','3Txt');
      set(t3Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(t3,4));	
   set(gcbo,'UserData',t3);
   
  	Str = get(gcbo,'String');				
	t3 = ConvStr(Str);
     

case 'tid4'							
   Str = get(gcbo,'String');				
   t4temp = ConvStr(Str);
   
   if(t4temp >= 0 & t4temp <= tmax)
      t4 = t4temp;
   else
   	t4Ptr = findobj(gcbf,'Tag','4Txt');
      set(t4Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(t4,4));	
   set(gcbo,'UserData',t4);
   
  	Str = get(gcbo,'String');				
	t4 = ConvStr(Str);
      
   
case 'tid5'							
   Str = get(gcbo,'String');				
   t5temp = ConvStr(Str);
   
   if(t5temp >= 0 & t5temp <= tmax)
      t5 = t5temp;
   else
   	t5Ptr = findobj(gcbf,'Tag','5Txt');
      set(t5Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(t5,4));	
   set(gcbo,'UserData',t5);
   
  	Str = get(gcbo,'String');				
	t5 = ConvStr(Str);
   
   
case 'tid6'							
   Str = get(gcbo,'String');				
   t6temp = ConvStr(Str);
   
   if(t6temp >= 0 & t6temp <= tmax)
      t6 = t6temp;
   else
   	t6Ptr = findobj(gcbf,'Tag','6Txt');
      set(t6Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(t6,4));	
   set(gcbo,'UserData',t6);
   
  	Str = get(gcbo,'String');				
	t6 = ConvStr(Str);
   
   
case 'xmin2d'							
   Str = get(gcbo,'String');				
   x1temp = ConvStr(Str);
   
   if(x1temp >= 0 & x1temp <= x2)
      x1 = x1temp;
   else
   	x1Ptr = findobj(gcbf,'Tag','X1Txt');
      set(x1Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(x1,4));	
   set(gcbo,'UserData',x1);
   
  	Str = get(gcbo,'String');				
	x1 = ConvStr(Str);

               
case 'xmax2d'							
   Str = get(gcbo,'String');				
   x2temp = ConvStr(Str);
   
   if(x2temp >= 0 & x2temp >= x1)
      x2 = x2temp;
   else
   	x2Ptr = findobj(gcbf,'Tag','X2Txt');
      set(x2Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(x2,4));	
   set(gcbo,'UserData',x2);
   
  	Str = get(gcbo,'String');				
	x2 = ConvStr(Str);

     
case 'procmin2d'
   Str = get(gcbo,'String');				
   proc1temp = ConvStr(Str);
  
   if(proc1temp <= proc2)
      proc1 = proc1temp;
   else
      proc1Ptr = findobj(gcbf,'Tag','Proc1Txt');
      set(proc1Ptr,'String','UserData');
   end   
      
   set(gcbo,'String',num2str(proc1,4));	
   set(gcbo,'UserData',proc1);
      
   Str = get(gcbo,'String');				
   proc1 = ConvStr(Str);
   
   
case 'procmax2d'
   Str = get(gcbo,'String');				
   proc2temp = ConvStr(Str);
  
   if(proc2temp >= proc1)
      proc2 = proc2temp;
   else
      proc2Ptr = findobj(gcbf,'Tag','Proc2Txt');
      set(proc2Ptr,'String','UserData');
   end   
      
   set(gcbo,'String',num2str(proc2,4));	
   set(gcbo,'UserData',proc2);
      
   Str = get(gcbo,'String');				
   proc2 = ConvStr(Str);
   
   
case 'visa1'
   figure
   plot(x(:),konc(:,10*t1),'r-',x(:),konc(:,10*t2),'b-',x(:),konc(:,10*t3),'g-',x(:),konc(:,10*t4),'y-',x(:),konc(:,10*t5),'m-',x(:),konc(:,10*t6),'k-');
   axis([x1 x2 proc1 proc2]);
   xlabel('Distance from lung periphery (cm)');
   ylabel('Gas conc (%)');
   legend
   
   
case 'tid13d'							
   Str = get(gcbo,'String');				
   t7temp = ConvStr(Str);
   
   if(t7temp >= 0 & t7temp <= tmax)
      t7 = t7temp;
   else
   	t7Ptr = findobj(gcbf,'Tag','7Txt');
      set(t7Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(t7,4));	
   set(gcbo,'UserData',t7);
   
  	Str = get(gcbo,'String');				
   t7 = ConvStr(Str);
   
   
 case 'tid23d'							
   Str = get(gcbo,'String');				
   t8temp = ConvStr(Str);
   
   if(t8temp >= 0 & t8temp <= tmax)
      t8 = t8temp;
   else
   	t8Ptr = findobj(gcbf,'Tag','8Txt');
      set(t8Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(t8,4));	
   set(gcbo,'UserData',t8);
   
  	Str = get(gcbo,'String');				
   t8 = ConvStr(Str);
  
  
case 'xmin3d'							
   Str = get(gcbo,'String');				
   x3temp = ConvStr(Str);
   
   if(x3temp >= 0 & x3temp <= x4)
      x3 = x3temp;
   else
   	x3Ptr = findobj(gcbf,'Tag','X3Txt');
      set(x3Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(x3,4));	
   set(gcbo,'UserData',x3);
   
  	Str = get(gcbo,'String');				
	x3 = ConvStr(Str);
  
  
case 'xmax3d'							
   Str = get(gcbo,'String');				
   x4temp = ConvStr(Str);
   
   if(x4temp >= 0 & x4temp >= x3)
      x4 = x4temp;
   else
   	x4Ptr = findobj(gcbf,'Tag','X4Txt');
      set(x4Ptr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(x4,4));	
   set(gcbo,'UserData',x4);
   
  	Str = get(gcbo,'String');				
	x4 = ConvStr(Str);
     

case 'visa2'
   index1 = find(x <= x4);
   index2 = find(x <= x3);
   figure
   mesh(t7+0.1:0.1:t8,x(index1(1):index2(1)),konc(index1(1):index2(1),10*t7+1:10*t8))
   xlabel('Time (s)');
   ylabel('Distance from lung periphery (cm)');
   zlabel('Gas conc (%)');
   
   
end   
   


%-------------------------------------------------------------------
% Underfunction wich converts string to a number.
%-------------------------------------------------------------------
function OutVal = ConvStr(InStr)

temp = str2num(InStr);              	% Converts string to number.

if (isempty(temp))							% String not a number,
   OutVal = get(gcbo,'UserData');			% use previous value.
else
   OutVal = temp;
end
%-------------------------------------------------------------------
% End of underfunction
%-------------------------------------------------------------------

