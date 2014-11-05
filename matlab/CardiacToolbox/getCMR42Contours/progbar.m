function [handle]=progbar(varargin)

% progbar.m: general purpose progress bar for MATLAB
% type 'progbar demo' to get started
%
% pb=progbar(info)         initialises a new progbar
%    progbar(pb,progress)  updates an existing progbar
%    progbar(pb,-1)        closes an existing progbar
%    close(pb)             alternative way to close
%
% info, if supplied, should be a structure with some
% or all of the following fields
%
%   title       title, default='Progress'
%
%   msg         message, default=no message
%
%   size        a small integer dictating the size of the bar,
%				default=1
%
%   period      a time in seconds dictating the minimum update
%               rate of the progress bar, default=0
%
%	pos			any of 'topleft,top,top right,centreleft,centre,
%				centreright,bottomleft,bottom,bottomright'. note
%				american spelling accepted, but frowned upon.
%
%	color		RGB triple in MATLAB format, [R G B]. if you have
%				the color tool 'rgb.m' on your path, you may use
%				any format that rgb accepts.
%
%	clearance	number between 0 and 1, affects appearance of the
%				bar. default is 0.

%	VERSION:	04/10/2001
%	AUTHOR:		stupidkid
%	CONTACT:	stupidkid@benmitch.co.uk
%	WWW:		www.benmitch.co.uk\matlab
%	LOCATION:	interaction\guis\

if nargin==0 info.dummy='dummy'; end
if nargin==1 info=varargin{1}; end

if nargin<2
    % check type
	if isa(info,'char') & strcmp(info,'demo')
		progbar_demo
		return
	end
	
    if ~isstruct(info)
        error('Single argument should be an info structure or ''demo''');
    end
    
    % extract
    if ~isfield(info,'title') info.title='Progress'; end
    if ~isfield(info,'msg') info.msg=''; end
    if ~isfield(info,'size') info.size=1; end
    if ~isfield(info,'period') info.period=0; end
	if ~isfield(info,'pos') info.pos='centre'; end
	if ~isfield(info,'color') info.color=[0 0 1]; end
	if ~isfield(info,'clearance') info.clearance=0; end

    % check values
    if floor(info.size)~=info.size | info.size<1 | info.size>5
        error('size should be a small integer (1-5)');
    end
    if info.period<0
        error('period should not be negative');
    end
	if info.clearance<0 | info.clearance>1
		error('clearance should be between 0 and 1');
	end
    
    spars=get(0,'screensize');
	sl=spars(1);
	sb=spars(2);
	sw=spars(3);
	sh=spars(4);

	pw=300*info.size;			% progress bar width
	ph=20*info.size;			% progress bar height
	mh=~isempty(info.msg)*30;	% message bar height
	th=ph+mh;					% total height

	% check position
	switch info.pos
	case {'center','centre'}
		x=sl+sw/2;
		y=sb+sh/2;
	case {'centerleft','centreleft'}
		x=sl+sw/4;
		y=sb+sh/2;
	case {'centerright','centreright'}
		x=sl+3*sw/4;
		y=sb+sh/2;
	case {'top'}
		x=sl+sw/2;
		y=sb+3*sh/4;
	case {'topleft'}
		x=sl+sw/4;
		y=sb+3*sh/4;
	case {'topright'}
		x=sl+3*sw/4;
		y=sb+3*sh/4;
	case {'bottom'}
		x=sl+sw/2;
		y=sb+sh/4;
	case {'bottomleft'}
		x=sl+sw/4;
		y=sb+sh/4;
	case {'bottomright'}
		x=sl+3*sw/4;
		y=sb+sh/4;
	otherwise
		error(['pos ''' info.pos ''' was not recognised']);
	end
	
	olddfp=get(0,'defaultfigureposition');
	set(0,'defaultfigureposition',[x-pw/2 y-th/2 pw th])
	info.lastclock=clock;
	handle=figure(...
	'MenuBar','none',...
	'resize','off',...
	'numbertitle','off',...
	'userdata',info,...
	'name',info.title);
	set(0,'defaultfigureposition',olddfp)

	% get color
	if ~isempty(which('rgb.m'))
		color=rgb(info.color);
	else
		color=info.color;
	end
	
	% check color
	if ~isa(color,'double') | size(color,1)~=1 | size(color,2)~=3
		close(handle);
		error('Unrecognised color');
	end
	
    rectangle('position',[0 0 0.001 1],'facecolor',color,'edgecolor',color,'visible','off')
    set(gca,'position',[0 0 1 ph/(ph+mh)])
    axis([0 100 -info.clearance 1+info.clearance])
    set(gca,'Xtick',0:10:100);
    set(gca,'Ytick',[-1 2]);
    set(gca,'box','on');

    title(info.msg);
    drawnow
elseif nargin==2
	pb=varargin{1};
	progress=varargin{2};
	
	% check types
	if ~isnumeric(pb)
		error('pb should be a progress bar handle');
	end
	if ~isnumeric(progress)
		error('progress should be a number');
	end
	
	% close
    if(progress==-1)
		try close(pb); end
		return;
	end
    
    info=get(pb,'userdata');
    if (info.period>0)
        if etime(clock,info.lastclock)<info.period return; end
        info.lastclock=clock;
        set(pb,'userdata',info);
    end
    
    if(progress<=0) progress=0.001; end
	if(progress>100) progress=100; end
	
    BarHandle=get(get(pb,'Children'),'Children');
    set(BarHandle,'position',[0 0 progress 1],'visible','on');
    progress=ceil(progress-0.5);
    set(pb,'name',[int2str(progress) '% - ' info.title]);
    drawnow
else
    error('Wrong number of input args');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function progbar_demo

disp([10 'progbar.m demonstration: press a key to move on' 10]);
disp('Simple usage:');
disp([10 'p=progbar;'])
p=progbar;
pause
disp(['progbar(p,20)']);
progbar(p,20);
pause
disp(['progbar(p,70)']);
progbar(p,70);
pause
disp(['progbar(p,100)']);
progbar(p,100);
pause
disp(['close(p)']);
close(p)
pause
disp([10 'Set some parameters:']);
disp([10 'info.color=[1 0 0]']);
disp(['info.title=''LMS''']);
disp(['info.size=2']);
disp([10 'p=progbar(info)']);
disp(['progbar(p,50)']);
info.color=[1 0 0];
info.title='LMS';
info.size=2;
p=progbar(info);
progbar(p,50);
pause
disp('close(p)')
close(p)
pause
disp([10 'Set some more parameters:']);
disp([10 'info.clearance=0.5']);
disp('info.msg=''Processing Least Mean Square algorithm''');
disp('p=progbar(info);');
disp('progbar(p,20)');
info.clearance=0.5;
info.msg='Processing Least Mean Square algorithm';
p=progbar(info);
progbar(p,20);
pause
disp('close(p)')
close(p)
pause
disp([10 'For other uses of color and period, and detailed help, type help progbar' 10]);


