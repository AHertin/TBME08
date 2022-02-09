function CbLungsim(action)


global dt npx f0 t0 tmax vtot dmol ru sc inito2 gb diff Model;


switch(action)
      
case 'singeltrumpet'
   set(gcbo,'Value',1);
   set(gcbo,'UserData',1);
   
   DubbelPtr = findobj(gcbf,'Tag','DubbelRadio');
   set(DubbelPtr,'Value',0);
   
   ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
   set(ForgreningPtr,'String','');
   
   DifferensPtr = findobj(gcbf,'Tag','Differens');
   set(DifferensPtr,'String','');

   Model = 1;
   
case 'dubbeltrumpet'
   set(gcbo,'Value',1);
   set(gcbo,'UserData',1);
   
   SingelPtr = findobj(gcbf,'Tag','SingelRadio');
   set(SingelPtr,'Value',0);
   
   ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
   set(ForgreningPtr,'String',num2str(8,4));
   gb = 9;
   
   DifferensPtr = findobj(gcbf,'Tag','Differens');
   set(DifferensPtr,'String',num2str(3,4));
   
   Model = 2;
		
		
case 'andningsfrekvens'						
   Str = get(gcbo,'String');				
   f0temp = ConvStr(Str);
   
   if(f0temp > 0)
      f0 = f0temp;
      t0 = 60/(60*f0);
   else
      AndningsFrekvensPtr = findobj(gcbf,'Tag','AndningsFrekvens');
      set(AndningsFrekvensPtr,'String','UserData');
   end
   
   	set(gcbo,'String',num2str(f0,4));	
   	set(gcbo,'UserData',f0);
      
      Str = get(gcbo,'String');				
	   f0 = ConvStr(Str);
      
case 'simuleringstid'                    
   Str = get(gcbo,'String');				
   tmaxtemp = ConvStr(Str);						
   
   if(tmaxtemp >= 0.6)
   	tmax = tmaxtemp;
   else
      SimuleringsTidPtr = findobj(gcbf,'Tag','SimuleringsTid');
   	set(SimuleringsTidPtr,'String','UserData');
	end	
   
   set(gcbo,'String',num2str(tmax,4));	
   set(gcbo,'UserData',tmax);				
  
  	Str = get(gcbo,'String');				
	tmax = ConvStr(Str);

case 'tidalvolym'							
   Str = get(gcbo,'String');				
   vtottemp = ConvStr(Str);
   
   if(vtottemp >= 0)
      vtot = vtottemp;
   else
   	TidalVolymPtr = findobj(gcbf,'Tag','TidalVolym');
      set(TidalVolymPtr,'String','UserData');
   end
   
   set(gcbo,'String',num2str(vtot,4));	
   set(gcbo,'UserData',vtot);
   
  	Str = get(gcbo,'String');				
	vtot = ConvStr(Str);

case 'initialo2'
   Str = get(gcbo,'String');
   inito2temp = ConvStr(Str);
   
   if(inito2temp >= 0)
      inito2 = inito2temp;
   else
      InitialPtr = findobj(gcbf,'Tag','InitialO2');
   	set(InitialPtr,'String','UserData');
	end
    
   set(gcbo,'String',num2str(inito2,4));
   set(gcbo,'UserData',inito2);
   
  	Str = get(gcbo,'String');				
	inito2 = ConvStr(Str);
   
case 'diffusionskoefficient'
   Str = get(gcbo,'String');
   dmoltemp = ConvStr(Str);
   
   if(dmoltemp >= 0)
      dmol = dmoltemp;
   else
      DiffusionPtr = findobj(gcbf,'Tag','DiffusionsKoefficient');
   	set(DiffusionPtr,'String','UserData');
   end   
      
   set(gcbo,'String',num2str(dmol,4));
   set(gcbo,'UserData',dmol);
   
   Str = get(gcbo,'String');				
	dmol = ConvStr(Str);

case 'forgrening'
   DubbelPtr = findobj(gcbf,'Tag','DubbelRadio');
   DubbelVal = get(DubbelPtr,'Value');
   
   if(DubbelVal == 1)
   	Str = get(gcbo,'String');
   	gbtemp = ConvStr(Str)+1;
   
 		if (mod(gbtemp,1)==0 & (gbtemp <=23) & (gbtemp >=0))
      	gb = gbtemp;
   	else
      	ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
   		set(ForgreningPtr,'String','UserData');
   	end
   
   	set(gcbo,'String',num2str(gb-1,4));
   	set(gcbo,'UserData',gb-1);
   
   	Str = get(gcbo,'String');
   	gb = ConvStr(Str)+1;
   
	else
      ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
   	set(ForgreningPtr,'String','');
	end
      
case 'differens'  
   DubbelPtr = findobj(gcbf,'Tag','DubbelRadio');
   DubbelVal = get(DubbelPtr,'Value');
   
   if(DubbelVal == 1)
		Str = get(gcbo,'String');
  		diff = ConvStr(Str);
   	set(gcbo,'String',num2str(diff,4));
      set(gcbo,'UserData',diff);
   else
     	DifferensPtr = findobj(gcbf,'Tag','Differens');
   	set(DifferensPtr,'String',''); 
   end
   
case 'aterstall'
   SingelPtr = findobj(gcbf,'Tag','SingelRadio');
   set(SingelPtr,'Value',1);
   Model = 1;
   
   DubbelPtr = findobj(gcbf,'Tag','DubbelRadio');
   set(DubbelPtr,'Value',0);
   
   AndningsFrekvensPtr = findobj(gcbf,'Tag','AndningsFrekvens');
   set(AndningsFrekvensPtr,'String',num2str(0.2,4));
   f0 = 0.2;
   t0 = 60/(60*f0);
   
   SimuleringsTidPtr = findobj(gcbf,'Tag','SimuleringsTid');
   set(SimuleringsTidPtr,'String',num2str(5,4));
   tmax = 5;
   
   TidalVolymPtr = findobj(gcbf,'Tag','TidalVolym');
   set(TidalVolymPtr,'String',num2str(500,4));
   vtot = 500; 
   
   InitialPtr = findobj(gcbf,'Tag','InitialO2');
   set(InitialPtr,'String',num2str(0,4));
   inito2 = 0;

   DiffusionPtr = findobj(gcbf,'Tag','DiffusionsKoefficient');
   set(DiffusionPtr,'String',num2str(0.225,4));
   dmol = 0.225;
   
   ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
   set(ForgreningPtr,'String','');
   gb = 9;
   
   DifferensPtr = findobj(gcbf,'Tag','Differens');
   set(DifferensPtr,'String','');
   diff = 3;
   
%    X4TxtPtr = findobj(gcbf,'Tag','X4Txt');
%    set(X4TxtPtr,'String','');
%    diff = 1.5;
    
%    edit30Ptr = findobj(gcbf,'Tag','edit30');
%    set(edit30Ptr,'String','1');
%    set(edit30Ptr,'UserData','1');
%    t1 = 1;
   
case 'starta'
   if Model == 1
      npx = 400;
      dt = 0.01;
      Enkeltrumpet; 
   end
   
   if Model == 2
      npx = 400; 
      dt = 0.01;
      Dubbeltrumpet;
   end
   
   
   case 'save'
      Spara;
      
      
   case 'load'
      Ladda;
      
      if Model == 1      
         SingelPtr = findobj(gcbf,'Tag','SingelRadio');
         set(SingelPtr,'Value',1);
      
         DubbelPtr = findobj(gcbf,'Tag','DubbelRadio');
         set(DubbelPtr,'Value',0);
      else
         SingelPtr = findobj(gcbf,'Tag','SingelRadio');
         set(SingelPtr,'Value',0);
      
         DubbelPtr = findobj(gcbf,'Tag','DubbelRadio');
         set(DubbelPtr,'Value',1);
      end
   
      AndningsFrekvensPtr = findobj(gcbf,'Tag','AndningsFrekvens');
      set(AndningsFrekvensPtr,'String',num2str(f0,4));
      t0 = 60/(60*f0);
   
      SimuleringsTidPtr = findobj(gcbf,'Tag','SimuleringsTid');
      set(SimuleringsTidPtr,'String',num2str(tmax,4));
   
      TidalVolymPtr = findobj(gcbf,'Tag','TidalVolym');
      set(TidalVolymPtr,'String',num2str(vtot,4));
   
      InitialPtr = findobj(gcbf,'Tag','InitialO2');
      set(InitialPtr,'String',num2str(inito2,4));

      DiffusionPtr = findobj(gcbf,'Tag','DiffusionsKoefficient');
      set(DiffusionPtr,'String',num2str(dmol,4));
   
      if Model == 1
         ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
         set(ForgreningPtr,'String','');
   
         DifferensPtr = findobj(gcbf,'Tag','Differens');
         set(DifferensPtr,'String','');
      else 
         ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
         set(ForgreningPtr,'String',num2str(gb-1,4));
   
         DifferensPtr = findobj(gcbf,'Tag','Differens');
         set(DifferensPtr,'String',num2str(diff,4));
      end

        
   case 'visa'
      if Model == 1      
         SingelPtr = findobj(gcbf,'Tag','SingelRadio');
         set(SingelPtr,'Value',1);
      
         DubbelPtr = findobj(gcbf,'Tag','DubbelRadio');
         set(DubbelPtr,'Value',0);
      else
         SingelPtr = findobj(gcbf,'Tag','SingelRadio');
         set(SingelPtr,'Value',0);
      
         DubbelPtr = findobj(gcbf,'Tag','DubbelRadio');
         set(DubbelPtr,'Value',1);
      end
   
      AndningsFrekvensPtr = findobj(gcbf,'Tag','AndningsFrekvens');
      set(AndningsFrekvensPtr,'String',num2str(f0,4));
      t0 = 60/(60*f0);
   
      SimuleringsTidPtr = findobj(gcbf,'Tag','SimuleringsTid');
      set(SimuleringsTidPtr,'String',num2str(tmax,4));
   
      TidalVolymPtr = findobj(gcbf,'Tag','TidalVolym');
      set(TidalVolymPtr,'String',num2str(vtot,4));
   
      InitialPtr = findobj(gcbf,'Tag','InitialO2');
      set(InitialPtr,'String',num2str(inito2,4));

      DiffusionPtr = findobj(gcbf,'Tag','DiffusionsKoefficient');
      set(DiffusionPtr,'String',num2str(dmol,4));
   
      if Model == 1
         ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
         set(ForgreningPtr,'String','');
   
         DifferensPtr = findobj(gcbf,'Tag','Differens');
         set(DifferensPtr,'String','');
      else 
         ForgreningPtr = findobj(gcbf,'Tag','Forgrening');
         set(ForgreningPtr,'String',num2str(gb,4));
   
         DifferensPtr = findobj(gcbf,'Tag','Differens');
         set(DifferensPtr,'String',num2str(diff,4));
      end
  
      if Model == 1      
         Enkelplott;
      end
   
      if Model == 2
         Dubbelplott;
      end

   
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



