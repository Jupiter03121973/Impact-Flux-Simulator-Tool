function status = odeprog(t,~,flag,varargin)
%status = odebarplot(t,y,flag,varargin)
%   ODE progress display function with interrupt control
%   Displays a vertical bar plot that fills as the simulation
%   nears its completion.  Also displays time ellapsed and estimates
%   time remaining in the simulation.  To avoid computation burden
%   refreshes are limited to every 0.5 seconds.
%
% Tim Franklin
% Virginia Tech
% Jesse Norris
% Wake Forrest
% May 2006

global odeprogglobvar

if nargin < 3 || isempty(flag) 
    if(etime(clock,odeprogglobvar(8:13))>0.5)  
        tfin=odeprogglobvar(1)/86400;
        sstrt=odeprogglobvar(2:7);
        perc=t(end)/(tfin*86400);
        area(findobj('Tag','Pbar'),[t(end)/86400 tfin-t(end)/86400;t(end)/86400, tfin-t(end)/86400]);
        title(findobj('Tag','Pbar'),[num2str(perc*100) '%']);
        set(findobj('Tag','eltime'),'String',etimev(clock,sstrt));
        set(findobj('Tag','esttime'),'String',etimev(etime(clock,sstrt)/perc*(1-perc)));
        odeprogglobvar(8:13)=clock;
    end
else
    switch(flag)
    case 'init'  
        odeprogglobvar=zeros(1,13);
        odeprogglobvar(1)=t(end);
        odeprogglobvar(2:7)=clock;
        odeprogglobvar(8:13)=clock;
        tfin=odeprogglobvar(1)/86400;
        sstrt=odeprogglobvar(2:7);
        fig = figure(95);
        set(fig, 'Name', 'Simulation Progress','NumberTitle','off')
        set(fig,'Position',[300,300,500,500]);
        ax = axes(figure(95),'Position',[0.5,0.25,0.25,0.6],'Tag','Pbar');
        axis(ax,[1,2,0,tfin]);
        set(ax,'XTickLabel',[],'NextPlot','replacechildren');
        ylabel(ax,'Simulation Progress - Time (d)');
        title(ax,'0%');
        area(ax,[0 tfin;0 tfin]);
        uicontrol(fig,'Style', 'pushbutton', 'String', 'Abort','Position', [7 460 90 30], 'Callback', @abortButtonCallback)
        uicontrol(fig,'Style', 'text', 'String', 'Ellapsed Time','Position', [7 100 90 15])
        uicontrol(fig,'Style', 'text', 'Tag', 'eltime', 'String', etimev(clock,sstrt),'Position', [7 80 90 15])
        uicontrol(fig,'Style', 'text', 'String', 'Time Remaining','Position', [7 60 90 15])
        uicontrol(fig,'Style', 'text', 'Tag', 'esttime', 'String', num2str(inf),'Position', [7 40 90 15])
        pause(0.1);

    case 'done'    
        if(ishandle(95))
            close(95);
        end
    end
end
status = 0;
drawnow;
end 

function [S] = etimev(t1,t0)
%ETIMEV  Verbose Elapsed time.
%   ETIMEV(T1,T0) returns string of the days, hours, minutes, seconds that have elapsed 
%   between vectors T1 and T0.  The two vectors must be six elements long, in
%   the format returned by CLOCK:
%
%       T = [Year Month Day Hour Minute Second]
%   OR
%   ETIMEV(t), t in seconds
if(exist('t1')&&exist('t0')&&length(t1)>2&&length(t0)>2)
    t=etime(t1,t0);     %Time in seconds
    if(t<0)
        t=-t;
    end
elseif(length(t1)==1)
    t=t1;
else
    t=0;
end
days=floor(t/(24*60*60));
t=t-days*24*60*60;
hours=floor(t/(60*60));
t=t-hours*60*60;
mins=floor(t/60);
t=floor(t-mins*60);
if(days>0)
    S=[num2str(days) 'd ' num2str(hours) 'h ' num2str(mins) 'm ' num2str(t) 's'];
elseif(hours>0)
    S=[num2str(hours) 'h ' num2str(mins) 'm ' num2str(t) 's'];
elseif(mins>0)
    S=[num2str(mins) 'm ' num2str(t) 's'];
else
    S=[num2str(t) 's'];
end
end

function abortButtonCallback(~, ~)    
    global abortSimulation; % Add this line

    abortSimulation = true; % Set the flag to stop the simulation
end

