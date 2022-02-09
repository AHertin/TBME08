function CbDubbelplott(action)

global x konc cf cl lf ll rf rl t1 t2 x1 x2 proc1 proc2 tmax ;


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
   figure;
   plot(x(1:cl),konc(1:cl,10*t1),'r-',x(1:cl),konc(1:cl,10*t2),'b-',x(lf:ll),konc(lf:ll,10*t1),'r',x(lf:ll),konc(lf:ll,10*t2),'b',[x(cl) x(lf)],[konc(cl,10*t1) konc(lf,10*t1)],'r',[x(cl) x(lf)],[konc(cl,10*t2) konc(lf,10*t2)],'b');
   hold on
   plot(x(1:cl),konc(1:cl,10*t1),'r-',x(1:cl),konc(1:cl,10*t2),'b-',x(rf:rl),konc(rf:rl,10*t1),'r:',x(rf:rl),konc(rf:rl,10*t2),'b:',[x(cl) x(rf)],[konc(cl,10*t1) konc(rf,10*t1)],'r:',[x(cl) x(rf)],[konc(cl,10*t2) konc(rf,10*t2)],'b:');
   hold off
   axis([x1 x2 proc1 proc2]);
   xlabel('Distance from lung periphery (cm)');
   ylabel('Oxygen conc. (%)');
   legend('1','2');
   
   
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

