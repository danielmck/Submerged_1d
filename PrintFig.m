function PrintFig(filename,varargin)
    %  PrintFig(figno, filename) - Matlab 2014b and later version
    %  
	%  Chris Johnson (chris.johnson@manchester.ac.uk), 2010--2015
	%
    % filename argument should have no extension
    %
    % Options:
    %   'trim'      -       Trims a figure to the size of its contents
    %   'debug'     -       Outputs a pink border around the figure. Useful
    %                       to check that LaTeX substitutions aren't
    %                       overlapping the border, which makes the figure
    %                       a different size from that expected
    %
    %   'font-computermodern' - LaTeX Computer Modern font (default)
    %   'font-times' -      LaTeX Times font
    %   'font-helvetica' -  LaTeX Helvetica font
    %   'font-minion' -     LaTeX Minion font (requires MinionPro package)
    %   'font-palatino' -   LaTeX Palatino font (requires mathpazo package)
    %
    %   'small'     -       Sets fonts to LaTeX 'small' size
    %   'jfm'       -       Sets 'small' and 'font-times' flags; a fairly
    %                       good representation of what the JFM typesetters
    %                       use, although not identical since they use a 
    %                       non-free version of the Times font.
    %
    %   'pdf'       -       output a PDF file (default)
    %   'eps'       -       output an EPS file. When used without 'pdf',
    %                       only an EPS is output
    %
    %
    %  NOTE: Legends are likely to be of the wrong size, since the PSFrag
    %  substitution text used internally is used by MATLAB to set the
    %  legend size. The 'UserData' field of a legend is used to add 'X'
    %  characters to the PSFrag text, and this can be used to widen the
    %  legend. Usage is, for example:
    %
    %  legend({'Plot 1', 'Plot 2'}, 'UserData', 6);
    %
    %  which adds six characters of space to the legend. Some trial and
    %  error is likely needed. Some kind of manual setting of legend size
    %  is unavoidable, since MATLAB can't know the size of a LaTeX-typeset
    %  string.
    %  (Actually, we could do this by running LaTeX beforehand on all the
    %  legend strings and reading back box sizes from the LaTeX log or
    %  command line output.
    %  The appropriate LaTeX (returning width, height above baseline, and
    %  depth below baseline) is something like:
    %
    % \documentclass{article}
    % \usepackage[utf8x]{inputenc}
    % \newsavebox{\mybox}
    % \newcommand{\GetDimensions}[1]{\sbox{\mybox}{#1}%
    % \message{Dimensions: \the\wd\mybox, \the\ht\mybox, \the\dp\mybox}%
    % }
    % \begin{document}
    % \GetDimensions{Legend line 1}
    % \GetDimensions{Legend line 2}
    % ....
    % \end{document}
    %
    % )
	%

    newGraphics=false;
    v = sscanf (version, '%d.%d.%d');
    if ((v(1)==8 && v(2)>=4) || (v(1)>8)) % new graphics on R2014b and later
        newGraphics = true;
    end
    debug=0; % Turn debug on to enable pink bounding box
    noborder=0;
    font='computermodern';
    smallfont=0;
    generateEPS = 0;
    generatePDF = 1;
    
    generateFlagSet = 0;

    optargin = size(varargin,2);

    for i=1:optargin
        if (strcmp(varargin{i},'trim')==1)
            noborder=1;
        elseif (strcmp(varargin{i},'debug')==1)
            debug=1;
        elseif (strcmp(varargin{i},'font-times')==1)
            font='times';
        elseif (strcmp(varargin{i},'font-computermodern')==1)
            font='computermodern';
        elseif (strcmp(varargin{i},'font-minion')==1)
            font='minion';
        elseif (strcmp(varargin{i},'font-helvetica')==1)
            font='helvetica';
        elseif (strcmp(varargin{i},'font-palatino')==1)
            font='palatino';
        elseif (strcmp(varargin{i},'small')==1)
            smallfont=1;
        elseif (strcmp(varargin{i},'jfm')==1)
            smallfont=1;
            font='times';
        % generate EPS only if eps flag set, or both file types if both
        % eps and pdf flags set
        elseif (strcmp(varargin{i},'eps')==1)
            generateEPS=1;
            if (generateFlagSet == 0)
                generatePDF = 0;
            end
            generateFlagSet=1;
        elseif (strcmp(varargin{i},'pdf')==1)
            generatePDF=1;
            if (generateFlagSet == 0)
                generateEPS = 0;
            end
            generateFlagSet=1;
        end
    end

	% Stop MATLAB from trying to interpret LaTeX itself
    set(groot, 'DefaultTextInterpreter', 'none');
    set(groot, 'DefaultLegendInterpreter', 'none')

    pause(.5); % avoid stupid MATLAB race conditions
    v = sscanf (version, '%d.%d.%d') ;
    if (newGraphics)
        figno=get(gcf,'Number');
    else
        figno=gcf;
    end
    unis = get(gcf,'units');
    set(gcf,'Units','centimeters');
    posvec=get(gcf,'Position');
    width=posvec(3);
    set(gcf,'Units',unis);


    % Set fixed-width font for PSfrag substitutions so that 
    % number-width changes don't shift insertion points for LaTeX text.
    % This should probably be done in laprint
    hax = findobj(gcf,'Type','axes');
    set(hax,'FontName','FixedWidth');
    % Also using tabular figures in math mode, and math-mode tick lables,
    % to ensure that axis labels that say '1' look properly aligned.

    % Add corner annotations to stop bounding box shrinkage
    if (noborder==0)
        a=annotation('line',[.0001 .001],[.001 .0001],'Color',[.995 .995 .995]);
        set(a,'LineWidth',0.001);
        annotation('line',[.9999 .999],[.999 .9999],'Color',[.995 .995 .995]);
        set(a,'LineWidth',0.001);
    end

    % Visible bounding box for debugging
    if (debug==1)
        annotation('rectangle',[0 0 1 1],'Color',[1 .5 .5]);
    end
    
	% Run customised version of laprint   
    laprintChris(figno,filename, 'asonscreen', 'on', 'mathticklabels','on','head', 'off',...
				 'width', width, 'scalefonts','off','factor', 1.0, 'figcopy','off',...
				 'printcmd',sprintf('print -f%d -depsc2 -loose -r4000 -painters %s.eps',figno,filename));
	% 'keepfontprops', 'on',

	% Open up the exported EPS and tweak line caps, dash styles and encoding.
    fid = fopen([filename '.eps'], 'r');
    str = fread(fid);
    fclose(fid);
    str = char(str');
    str = strrep(str, 'ISOLatin1Encoding', 'StandardEncoding'); % This for old MATLAB only
    str = strrep(str, '0 cap', '1 cap');% This for old MATLAB only
    str = strrep(str, '2 cap', '1 cap');% This for old MATLAB only
    
    str = strrep(str, '0 setlinecap', '1 setlinecap');
    str = strrep(str, '2 setlinecap', '1 setlinecap');
    
    str = strrep(str, '2 LJ', '1 LJ 1 setlinecap'); % Added for modern MATLAB to get rounded dashes etc.
    
    % This one for old MATLAB only
    str = strrep(str, '/DA { [6 dpi2point mul] 0 setdash } bdef', '/DA { [3 dpi2point mul] 0 setdash } bdef');
    fid = fopen([filename '.eps'], 'w');
    fprintf(fid, '%s', str);
    fclose(fid);

	% Move psfrag subs file
    system(sprintf('mv %s.tex %s-PrintFig-Subs.tex',filename,filename));

    % Set up some filenames...
    [~,baseName]=system(['basename ' filename]);
    baseName=strtrim(baseName);
    [~,dirName]=system(['dirname ' filename]);
    dirName=strtrim(dirName);
    firstStepName=[filename '-PrintFig'];
    firstStepBaseName=[baseName '-PrintFig'];

	% Construct a LaTeX file that includes the figure and psfrag subs
    fid = fopen(sprintf('%s.tex',firstStepName), 'w');
    fprintf(fid, '\\documentclass[a2paper]{article}\n');
    fprintf(fid, '\\pagestyle{empty}\n');
    fprintf(fid, '\\usepackage[margin=0.2cm]{geometry}\n');
    fprintf(fid, '\\usepackage{amsmath,amssymb}\n');
    fprintf(fid, '\\usepackage{graphicx}\n');
    fprintf(fid, '\\usepackage{psfrag}\n');
    fprintf(fid, '\\usepackage{xcolor}\n');

    if (strcmp(font,'times'))
        fprintf(fid, '\\usepackage{mathptmx}\n');
        fprintf(fid, '\\usepackage[T1]{fontenc}\n');
    elseif (strcmp(font,'minion'))
        fprintf(fid, '\\usepackage[mathlf,mathtabular,swash,footnotefigures]{MinionPro}\n');
        fprintf(fid, '\\usepackage[T1]{fontenc}\n');
    elseif (strcmp(font,'palatino'))
        fprintf(fid, '\\usepackage{mathpazo}\n');
    elseif (strcmp(font,'helvetica'))
        fprintf(fid, '\\usepackage[scaled]{helvet}\n');
        fprintf(fid, '\\renewcommand*\\familydefault{\\sfdefault}\n');
        fprintf(fid, '\\usepackage[T1]{fontenc}\n');
        fprintf(fid, '\\usepackage{sfmath}\n');
    end

    fprintf(fid, '\\newcommand{\\df}{\\phantom{0}}\n');
    [status, result]=system(sprintf('dirname %s',filename));
    fprintf(fid, '\\graphicspath{{%s/}}\n',result(1:(numel(result)-1)));
    fprintf(fid, '\\begin{document}\n');
    if (smallfont==1)
        fprintf(fid,'\\small\n');
    end
    fprintf(fid, '\\input{%s-PrintFig-Subs}\n',filename);
    fprintf(fid, '\\end{document}\n');

    fclose(fid);

	% Compile this LaTeX file
    commandLine=['latex --interaction=nonstopmode -output-directory ' dirName ' ' firstStepName '.tex'];
    [s,retval]=system(commandLine);
    if (s)
        disp(['LaTeX failed: ' retval]);
        close(gcf);
        return;
    end

	% Run dvips
    commandLine=['cd ' dirName '; dvips -q -t a2 -o ' firstStepBaseName '.ps ' firstStepBaseName '.dvi'] ;
    [s,retval]=system(commandLine);
    if (s)
        fprintf(['dvips failed: ' retval]);
        close(gcf);
        return;
    end

	% If generating an EPS, convert the resulting .ps to an eps, erasing the matlab-exported EPS
    if (generateEPS)
         commandLine=['ps2eps < ' firstStepName '.ps > ' filename '.eps '];
         [s,retval]=system(commandLine);
         if (s)
             fprintf(['ps2eps failed: ' retval]);
             close(gcf);
             return;
         end

    end
    
	% If generating a PDF, convert the resulting .ps to a PDF and crop it
    if (generatePDF)
        commandLine=['ps2pdf -dPDFSETTINGS=/prepress -sPAPERSIZE=a2 ' firstStepName '.ps ' firstStepName '.pdf'];
        [s,retval]=system(commandLine);
        if (s)
            fprintf(['ps2pdf failed: ' retval]);
            close(gcf);
            return;
        end

        commandLine=['pdfcrop --hires ' firstStepName '.pdf  ' filename '.pdf']; % --margins 1
        %commandLine=['mv ' firstStepName '.pdf  ' filename '.pdf']; 
        [s,retval]=system(commandLine);
        if (s)
            fprintf(['pdfcrop failed: ' retval]);
            close(gcf);
            return;
        end
    end

	% Remove temporary files
    commandLine=['rm ' firstStepName '.tex ' filename '-PrintFig-Subs.tex '  firstStepName '.aux ' firstStepName '.log ' firstStepName '.dvi ' firstStepName '.ps ' firstStepName '.pfg'];
    if (generatePDF) % uncropped pstopdf output
         commandLine = [commandLine ' ' firstStepName '.pdf'];
    end
    if (~generateEPS) % if EPS generation is off, this file will still exist (it is the original MATLAB exported EPS) 
         commandLine = [commandLine ' ' filename '.eps']; 
    end
    [s,retval]=system(commandLine);
    if (s)
       fprintf(['rm failed: ' retval]);
        close(gcf);
        return;
    end

	% Close the figure; avoids the remanants of substituted text polluting the next figure if hold is on
    close(gcf);

end

function laprintChris(figno,filename,varargin)
%LAPRINTChris prints a figure for inclusion in LaTeX documents.
% 
% cgj: This is a custom version of Laprint designed to be used with
%      PrintFig, in Matlab 2014b and above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Initialize
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

laprintident = '3.16 cgj (22.v.2015)';   
vers = version;
vers = eval(vers(1:3));
if vers < 6.1
  error('Sorry. Matlab 6.1 or above is required.')
end

% no output
if nargout
   error('No output argument, please.')
end  

inter=get(groot,'defaulttextinterpreter');
if ~strcmp(inter,'none')
  warning('LaPrint:general',['It is recommended to switch off the '...
          'text interpreter\nbefore creating a figure to be saved '...
          'with LaPrint. Use the command\n',...
          '   >> set(groot,''defaulttextinterpreter'',''none'').'])
end

inter=get(groot,'defaultlegendinterpreter');
if ~strcmp(inter,'none')
  warning('LaPrint:general',['It is recommended to switch off the '...
          'text interpreter\nbefore creating a figure to be saved '...
          'with LaPrint. Use the command\n',...
          '   >> set(groot,''defaultlegendinterpreter'',''none'').'])
end


% nargin >=1 and ~isa(figno,'char')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% PART 1 of advanced usage:
%%%% Check inputs and initialize
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get settings 
LAPRINTOPT = prefsettings;
if isempty(LAPRINTOPT)
 LAPRINTOPT = factorysettings;
end

% modify prefs
if ~isa(figno,'double') 
figno
error('This is not a figure handle.') 
end
if ~any(get(0,'children')==figno)
figno
error('This is not a figure handle.') 
end
LAPRINTOPT.figno = figno;

if nargin>1
if ~isa(filename,'char')  
  filename
  error('This is not a file name.') 
end
LAPRINTOPT.filename=filename; 
end    


% read and check command line options  

try  % try old usage (Version 2.03)
   if nargin <=2
      error('lets take new usage')
   end
   % 2.03 defaults
   width           = 12;
   factor          = 0.8;
   scalefonts      = 1;
   keepfontprops   = 0;
   asonscreen      = 0;
   keepticklabels  = 0;
   mathticklabels  = 0;
   head            = 1;
   comment         = '';
   caption         = '';
   extrapicture    = 1;
   nzeros          = 5;
   verbose         = 0;
   figcopy         = 1;
   printcmd        = ['print(''-f<figurenumber>'',' ...
                      '''-deps'',''<filename.eps>'')'];
   package         = 'epsfig';
   color           = 0;
   createview      = 0;
   viewfilename    = [filename '_'];
   processview     = 0;
   cmd1            = '';
   cmd2            = '';
   cmd3            = '';
   cmd4            = '';
   cmd5            = '';
   cmd6            = '';
   cmd7            = '';
   cmd8            = '';
   for i=1:nargin-2
    if ~isa(varargin{i},'char')
      error('Options must be character arrays.')
    end  
    oriopt=varargin{i}(:)';
    opt=[ lower(strrep(oriopt,' ','')) '                   ' ];
    if strcmp(opt(1:7),'verbose')
      verbose=1;
    elseif strcmp(opt(1:10),'asonscreen')
      asonscreen=1;
    elseif strcmp(opt(1:14),'keepticklabels')
      keepticklabels=1;
    elseif strcmp(opt(1:14),'mathticklabels')
      mathticklabels=1;
    elseif strcmp(opt(1:13),'keepfontprops')
      keepfontprops=1;
    elseif strcmp(opt(1:14),'noextrapicture')
      extrapicture=0;
    elseif strcmp(opt(1:14),'noextrapicture')
      extrapicture=0;
    elseif strcmp(opt(1:5),'loose')
      printcmd = ['print(''-f<figurenumber>'',' ...
                      '''-deps'',''-loose'',''<filename.eps>'')'];
    elseif strcmp(opt(1:9),'nofigcopy')
      figcopy=0;
    elseif strcmp(opt(1:12),'noscalefonts')
      scalefonts=0;
    elseif strcmp(opt(1:6),'nohead')
      head=0;
    elseif strcmp(opt(1:7),'caption')
      eqpos=findstr(oriopt,'=');
      if isempty(eqpos)
	    caption='Matlab Figure';
      else	
	    caption=oriopt(eqpos+1:length(oriopt));
      end	
    elseif strcmp(opt(1:8),'comment=')
      eqpos=findstr(oriopt,'=');
      comment=oriopt(eqpos(1)+1:length(oriopt));
    elseif strcmp(opt(1:9),'viewfile=')
      createview=1;
      eqpos=findstr(oriopt,'=');
      viewfilename=oriopt(eqpos(1)+1:length(oriopt));
    elseif strcmp(opt(1:6),'width=')
      eval([ opt ';' ]);
    elseif strcmp(opt(1:7),'factor=')
      eval([ opt ';' ]);
    else
      error([ 'Option ' varargin{i} ' not recognized.'])
    end   
  end
 
  warning('LaPrint:general',['You are using the old LaPrint '...
          'syntax. This syntax might not be supported in '...
          'future releases of LaPrint.'])

catch % old usage doesn't work, take new one
  
  % restore preferences / factory defaults
  width           = LAPRINTOPT.width;
  factor          = LAPRINTOPT.factor;
  scalefonts      = LAPRINTOPT.scalefonts;
  keepfontprops   = LAPRINTOPT.keepfontprops;
  asonscreen      = LAPRINTOPT.asonscreen;
  keepticklabels  = LAPRINTOPT.keepticklabels;
  mathticklabels  = LAPRINTOPT.mathticklabels;
  head            = LAPRINTOPT.head;
  comment         = LAPRINTOPT.comment;
  caption         = LAPRINTOPT.caption;
  extrapicture    = LAPRINTOPT.extrapicture;
  nzeros          = LAPRINTOPT.nzeros;
  verbose         = LAPRINTOPT.verbose;
  figcopy         = LAPRINTOPT.figcopy;
  printcmd        = LAPRINTOPT.printcmd;
  package         = LAPRINTOPT.package;
  color           = LAPRINTOPT.color;
  createview      = LAPRINTOPT.createview;
  viewfilename    = LAPRINTOPT.viewfilename;
  processview     = LAPRINTOPT.processview;
  cmd1            = LAPRINTOPT.cmd1;
  cmd2            = LAPRINTOPT.cmd2;
  cmd3            = LAPRINTOPT.cmd3;
  cmd4            = LAPRINTOPT.cmd4;
  cmd5            = LAPRINTOPT.cmd5;
  cmd6            = LAPRINTOPT.cmd6;
  cmd7            = LAPRINTOPT.cmd7;
  cmd8            = LAPRINTOPT.cmd8;

  if nargin > 2
    if rem(nargin,2)
      error('Option names/values must appear in pairs.')
    end    
    for i=1:2:nargin-2
      if ~isa(varargin{i},'char')
         error('Option name must be a character array.')
      end  
      opt = lower(strrep(varargin{i}(:)',' ',''));
      val = varargin{i+1}(:)';
      switch opt
        case 'options'      
          if isa(val,'char')
            if strcmp(val,'factory')
              val = factorysettings;
            else
              load(val)
              val = LAPRINTOPT;
            end
          end
          if ~isa(val,'struct')
            error('Value of options must be a structure array.')
          end  
          % no error checking here!
          width           = val.width;
          factor          = val.factor;
          scalefonts      = val.scalefonts;
          keepfontprops   = val.keepfontprops;
          asonscreen      = val.asonscreen;
          keepticklabels  = val.keepticklabels;
          mathticklabels  = val.mathticklabels;
          head            = val.head;
          comment         = val.comment;
          caption         = val.caption;
          extrapicture    = val.extrapicture;
          nzeros          = val.nzeros;
          verbose         = val.verbose;
          figcopy         = val.figcopy;
          printcmd        = val.printcmd;
          package         = val.package;
          color           = val.color;
          createview      = val.createview;
          viewfilename    = val.viewfilename;
          processview     = val.processview;
          cmd1            = val.cmd1;
          cmd2            = val.cmd2;
          cmd3            = val.cmd3;
          cmd4            = val.cmd4;
          cmd5            = val.cmd5;
          cmd6            = val.cmd6;
          cmd7            = val.cmd7;
          cmd8            = val.cmd8;
        case 'width'     
          if ~isa(val,'double')  
            error('Value of width must be a double.')
          end  
          width = val;  
        case 'factor'     
          if ~isa(val,'double')  
            error('Value of factor must be a double.')
          end  
          factor=val;  
        case 'scalefonts'
          scalefonts = value01(val,opt); 
        case 'keepfontprops'
          keepfontprops = value01(val,opt); 
        case 'asonscreen'     
          asonscreen = value01(val,opt);
        case 'keepticklabels'
          keepticklabels = value01(val,opt); 
        case 'mathticklabels'
          mathticklabels = value01(val,opt) ;
        case 'head'
          head = value01(val,opt); 
        case 'comment'
          if ~isa(val,'char')
            error('Value of comment must be a character array.')
          end
          comment = val;
        case 'caption'
          if ~isa(val,'char')
            error('Value of caption must be a character array.')
          end
          caption = val;
        case 'extrapicture'
          extrapicture = value01(val,opt); 
        case 'nzeros'     
          if ~isa(val,'double')  
            error('Value of nzeros must be a double.')
          end  
          nzeros = val;
        case 'verbose'
          verbose = value01(val,opt);
        case 'figcopy'
          figcopy = value01(val,opt); 
        case 'printcmd'
          if ~isa(val,'char')
            error('Value of printcmd must be a character array.')
          end
          printcmd = val;
        case 'package'
          if ~isa(val,'char')
            error('Value of package must be a character array.')
          end
          val = lower(strrep(val,' ',''));
          switch val
            case {'graphicx','epsfig'}
              % fine
            otherwise
              error('Value of package is unknown.')
          end  
          package = val;
        case 'color'
          color = value01(val,opt); 
        case 'createview'
          createview = value01(val,opt);
        case 'viewfilename'
          if ~isa(val,'char')
            error('Value of viewfilename must be a character array.')
          end
          viewfilename = val;
        case 'processview'
          processview = value01(val,opt); 
        case 'cmd1'
          if ~isa(val,'char')
            error('Value of cmd1 must be a character array.')
          end
          cmd1 = val;
        case 'cmd2'
          if ~isa(val,'char')
            error('Value of cmd2 must be a character array.')
          end
          cmd2 = val;
        case 'cmd3'
          if ~isa(val,'char')
            error('Value of cmd3 must be a character array.')
          end
          cmd3 = val;
        case 'cmd4'
          if ~isa(val,'char')
            error('Value of cmd4 must be a character array.')
          end
          cmd4 = val;
        case 'cmd5'
          if ~isa(val,'char')
            error('Value of cmd5 must be a character array.')
          end
          cmd5 = val;
        case 'cmd6'
          if ~isa(val,'char')
            error('Value of cmd6 must be a character array.')
          end
          cmd6 = val;
        case 'cmd7'
          if ~isa(val,'char')
            error('Value of cmd7 must be a character array.')
          end
          cmd7 = val;
        case 'cmd8'
          if ~isa(val,'char')
            error('Value of cmd8 must be a character array.')
          end
          cmd8 = val;
        otherwise
          error(['Option ''' opt ''' unknown'])
      end % switch opt
    end % for i=3:2:nargin 
  end % if nargin > 2
end % try / catch    

if verbose, 
  disp([ 'This is LaPrint, version ' laprintident '.' ]); 
end  

comment   = strrep(strrep(comment,'\','\\'),'%','%%');
caption   = strrep(strrep(caption,'\','\\'),'%','%%');
iscaption = logical(length(caption));

if nzeros < 3
  warning('LaPrint:general',...
          'The value of nzero should be >=3. I will use nzeros=3.')
  nzeros=3;  
end

if processview 
  createview=1;
end

if mathticklabels
  Do='$';
else  
  Do='';
end  

% eps- and tex- filenames
[epsfullnameext,epsbasenameext,epsbasename,epsdirname] = ...
                       getfilenames(filename,'eps',verbose);
[texfullnameext,texbasenameext,texbasename,texdirname] = ...
                       getfilenames(filename,'tex',verbose);
if ~strcmp(texdirname,epsdirname)
   warning('LaPrint:files',['The eps-file and tex-file are '...
          'placed in different directories.']);
end

if createview | processview
  [viewfullnameext,viewbasenameext,viewbasename,viewdirname] = ...
                       getfilenames(viewfilename,'tex',verbose);
  if strcmp(texfullnameext,viewfullnameext)
    viewfilename=[ viewfilename '_'];
    warning('LaPrint:files',['The tex- and view-file coincide. '... 
           'I''ll use '' ' viewfilename ' ''. Hope that''s ok.' ])
  end  
  [viewfullnameext,viewbasenameext,viewbasename,viewdirname]= ...
                       getfilenames(viewfilename,'tex',verbose);
  if ~strcmp(texdirname,viewdirname)
    warning('LaPrint:files',['The eps-file and view-file are '...
	   'placed in different directories.' ])
  end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% PART 2 of advanced usage:
%%%% Create new figure, insert tags, and bookkeep original text
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show all
shh = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');

% preparing check for copyobj bug
figno_ori = figno;
number_children_ori = length(get(figno_ori,'children'));

% open new figure (if required) and set properties
if figcopy
  figno = copyobj(figno,0);
  set(figno,'visible','off')
  set(figno,'Numbertitle','off')
  set(figno,'MenuBar','none')
  pause(0.5)
end  

if asonscreen  
  xlimmodeauto       = findobj(figno,'xlimmode','auto');
  xtickmodeauto      = findobj(figno,'xtickmode','auto');
  xticklabelmodeauto = findobj(figno,'xticklabelmode','auto');
  ylimmodeauto       = findobj(figno,'ylimmode','auto');
  ytickmodeauto      = findobj(figno,'ytickmode','auto');
  yticklabelmodeauto = findobj(figno,'yticklabelmode','auto');
  zlimmodeauto       = findobj(figno,'zlimmode','auto');
  ztickmodeauto      = findobj(figno,'ztickmode','auto');
  zticklabelmodeauto = findobj(figno,'zticklabelmode','auto');
  set(xlimmodeauto,'xlimmode','manual')
  set(xtickmodeauto,'xtickmode','manual')
  set(xticklabelmodeauto,'xticklabelmode','manual')
  set(ylimmodeauto,'ylimmode','manual')
  set(ytickmodeauto,'ytickmode','manual')
  set(yticklabelmodeauto,'yticklabelmode','manual')
  set(zlimmodeauto,'ylimmode','manual')
  set(ztickmodeauto,'ytickmode','manual')
  set(zticklabelmodeauto,'yticklabelmode','manual')
end  
set(figno,'paperunits','centimeters');
set(figno,'units','centimeters');
orip = get(figno,'Position');

% determine width and height
if factor <= 0
  factor = width/orip(3);
end 
latexwidth = width;
epswidth   = latexwidth/factor;
epsheight  = epswidth*orip(4)/orip(3);

set(figno,'PaperPosition',[0 0 epswidth epsheight ])
set(figno,'papersize',[epswidth epsheight])
%set(figno,'Position',[orip(1)+0.5 orip(2)-0.5 epswidth epsheight ]) % cgj commented out
%set(figno,'Name',[ 'To be printed; size: ' num2str(factor,3) ...
%  ' x (' num2str(epswidth,3) 'cm x ' num2str(epsheight,3) 'cm)' ])

% some warnings
if verbose
  if (epswidth<13) | (epsheight<13*0.75)
    warning('LaPrint:size',['The size of the eps-figure is quite '...
       'small. The text objects might not be properly set. '...
       'Reducing ''factor'' might help.'])
  end
  if latexwidth/epswidth<0.5
    warning('LaPrint:size',['The size of the eps-figure is large ' ...
           'compared to the latex figure. '...
           'The text size might be too small. '...
           'Increasing ''factor'' might help.'])
  end  
  if (orip(3)-epswidth)/orip(3) > 0.1
    warning('LaPrint:size',['The size of the eps-figure is much '...
            'smaller than the original '...
            'figure on screen. Matlab might save different ticks '...
            'and ticklabels than in the original figure. '...
            'See option ''asonscreen''.'])
  end
  disp('Strike any key to continue.');
  pause
end  

%
% TEXT OBJECTS: modify new figure 
%

% find all text objects
hxl = get(findobj(figno,'type','axes'),'xlabel');
hyl = get(findobj(figno,'type','axes'),'ylabel');
hzl = get(findobj(figno,'type','axes'),'zlabel');
hti = get(findobj(figno,'type','axes'),'title');
hte = findobj(figno,'type','text');
%cgj
htcb = get(findobj(figno,'type','ColorBar'),'Label');

% array of all text handles
htext = unique([ celltoarray(hxl) celltoarray(hyl) celltoarray(hzl) ...
      celltoarray(hti) celltoarray(hte) celltoarray(htcb) ]);
nt = length(htext);

% set(celltoarray(hxl),'VerticalAlignment','top');
% get alignments
hora  = get(htext,'HorizontalAlignment');
vera  = get(htext,'VerticalAlignment');
align = cell(nt,1);
for i=1:nt
  align{i} = hora{i}(1);
  switch vera{i}
  case 'top'
    align{i} = [align{i} 't'];
  case 'cap'
%     if ~isempty(get(htext(i),'string'))
%       warning('LaPrint:text',['Using vertical ' ...
%             'alignment ''top'' instead of ''cap''.'])
%     end  
    align{i} = [align{i} 't'];
  case 'middle'
    align{i} = [align{i} 'c'];
  case 'baseline'
    align{i} = [align{i} 'B'];
  case 'bottom'
    align{i} = [align{i} 'b'];
  otherwise
    warning('LaPrint:text',['Vertical alignment ' vera{i} ...
            ' unknown. Using ''c''.'])
    align{i} = [align{i} 'c'];
  end
end  

% generate new strings and store old ones
oldstr   = get(htext,'string');
newstr   = cell(nt,1);
basestr  = ['s' char(48*ones(1,nzeros-1))];
extrastr = 0;
for i=1:nt
  osi = oldstr{i};
  oldstr{i} = ['\setlength{\tabcolsep}{0pt}\begin{tabular}{' ...
          align{i}(1) '}'];
  isnonempty_osi = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Start of legend processing
  if strcmp(get(get(htext(i),'parent'),'tag'),'legend')
    newstr1 = [];  
   
    extraChars = 0;
    udval=get(get(htext(i),'parent'),'UserData')
    if (numel(udval)==1)
       extraChars = udval;
    end
    
    if isa(osi,'cell')
      % Legend/cell : Don't use tabular, employ extra strings 
      nlines = length(osi);
      if nlines > 1
        newstr{nt+extrastr+nlines-1} = [];
        oldstr{nt+extrastr+nlines-1} = [];
        htext((nt+extrastr+1):(nt+extrastr+nlines-1))=htext(i);
        for line=1:nlines-1  
          oldstr{nt+extrastr+line} = strrep(strrep(osi{line},'\','\\'),'%','%%'); 
          paddingStr = char(ones(1,extraChars)*88);
          newstr{nt+extrastr+line} = [overwritetail(basestr,nt+extrastr+line) paddingStr];
          newstr1 = [newstr1; overwritetail(basestr,nt+extrastr+line) paddingStr]; 
        end    
        extrastr = extrastr+nlines-1;  
      end    
      if nlines > 0
        oldstr{i} = strrep(strrep(osi{nlines},'\','\\'),'%','%%'); 
        newstr{i} = overwritetail(basestr,i);
        newstr1   = [newstr1; overwritetail(basestr,i)]; 
      end  
      % replace strings in figure
      set(htext(i),'string',cellstr(newstr1));
    else
      % Legend/matrix : Don't use tabular, employ extra strings 
      nlines=size(osi,1);
      if nlines > 1
        newstr{nt+extrastr+nlines-1} = [];
        oldstr{nt+extrastr+nlines-1} = [];
        htext((nt+extrastr+1):(nt+extrastr+nlines-1))=htext(i);
        for line=1:nlines-1  
          oldstr{nt+extrastr+line} = strrep(strrep(osi(line,:),'\','\\'),'%','%%'); 
          paddingStr = char(ones(1,extraChars)*88);
          newstr{nt+extrastr+line} = [overwritetail(basestr,nt+extrastr+line) paddingStr];
          newstr1 = [newstr1; overwritetail(basestr,nt+extrastr+line) paddingStr]; 
        end    
        extrastr = extrastr+nlines-1;  
      end    
      if nlines > 0
        oldstr{i} = strrep(strrep(osi(nlines,:),'\','\\'),'%','%%'); 
        newstr{i} = overwritetail(basestr,i);
        newstr1   = [newstr1; overwritetail(basestr,i)]; 
      end  
      % replace strings in figure
      set(htext(i),'string',newstr1);
    end
  else %%% end of legend processing
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % text, not a legend  
    if isa(osi,'cell')
      nlines = length(osi);
      if nlines > 1
        for line=1:nlines-1  
          oldstr{i}=[oldstr{i} osi{line} '\\'];
          isnonempty_osi = isnonempty_osi+length(osi{line});  
        end    
        if align{i}(2) == 'B'
          warning('LaPrint:text',['Vertical Alignment ''baseline'' '...
                  'in text with multiple rows might not match.'])
          align{i}(2) = 't';
        end  
      end    
      if nlines > 0
        oldstr{i} = [oldstr{i} osi{nlines} '\end{tabular}'];
        isnonempty_osi = isnonempty_osi+length(osi{nlines});
      end  
      oldstr{i} = strrep(strrep(oldstr{i},'\','\\'),'%','%%');
      if isnonempty_osi
        newstr{i} = overwritetail(basestr,i);
      else  
        newstr{i} = '';    
      end
      % replace strings in figure
      set(htext(i),'string',newstr{i}); 
    else
      nlines=size(osi,1);
      if nlines > 1
        for line=1:nlines-1  
          oldstr{i} = [oldstr{i} osi(line,:) '\\'];
          isnonempty_osi = isnonempty_osi+length(osi(line,:));  
        end    
        if align{i}(2) == 'B'
          warning('LaPrint:text',['Vertical Alignment ''baseline'' '...
                  'in text with multiple rows might not match.'])
          align{i}(2) = 't';
        end  
      end
      if nlines > 0
        oldstr{i} = [oldstr{i} osi(nlines,:) '\end{tabular}'];
        isnonempty_osi = isnonempty_osi+length(osi(nlines,:));
      end  
      oldstr{i} = strrep(strrep(oldstr{i},'\','\\'),'%','%%');
   
      if isnonempty_osi
        newstr{i} = overwritetail(basestr,i);
      else  
        newstr{i} = '';    
      end
      % replace string in figure
      set(htext(i),'string',newstr{i});
    end % isa cell  
  end % isa legend  
end % for

% cgj_begin
%%% New style (Matlab 2014b and beyond) legend processing
%%% Legends are no longer axes and the text is not text items
%%% Instead need to scan legend objects and pull strings out of their 'String' property
%%% Text styling within legends not yet implemented

legendObjs = findobj(figno,'Type','legend');
nLegends = numel(legendObjs);
for i=1:nLegends
    legendObjs(i).Interpreter='none';
    thisLegend = legendObjs(i);
    extraChars = 0;
    if (numel(thisLegend.UserData)==1)
       extraChars = thisLegend.UserData; 
    end
    nStringsInLegend = numel(thisLegend.String);
  
    newstr{nt+extrastr+nStringsInLegend} = [];
    oldstr{nt+extrastr+nStringsInLegend} = [];
    for j=1:nStringsInLegend
        thisString = thisLegend.String{j};
        oldstr{nt+extrastr+j} = strrep(strrep(thisString,'\','\\'),'%','%%'); 
        
        paddingStr = char(ones(1,extraChars)*88);
        newstr{nt+extrastr+j} = [overwritetail(basestr,nt+extrastr+j) paddingStr];
        thisLegend.String{j}=newstr{nt+extrastr+j};
    end

    extrastr = extrastr+nStringsInLegend;
end

colorbarObjs = findobj(figno,'Type','colorbar');
nColorBars = numel(colorbarObjs);
for i=1:nColorBars
   
    thisColorBar = colorbarObjs(i);
    nTicksInColorbar = numel(thisColorBar.TickLabels);
  
    newstr{nt+extrastr+nTicksInColorbar} = [];
    oldstr{nt+extrastr+nTicksInColorbar} = [];
    for j=1:nTicksInColorbar
        thisTickLabel = thisColorBar.TickLabels{j};
        % Use math mode for color bar ticks if math mode used for axis ticks
        % 'Do' is set to '$' when math axis labels are on 
        oldstr{nt+extrastr+j} = [Do strrep(strrep(thisTickLabel,'\','\\'),'%','%%') Do]; 
        
        newstr{nt+extrastr+j} = overwritetail(basestr,nt+extrastr+j);
        thisColorBar.TickLabels{j}=newstr{nt+extrastr+j};
    end

    extrastr = extrastr+nTicksInColorbar;
end
ntp = nt+extrastr;
% cgj_end


% Alignment of Legends
if extrastr > 0
  align{ntp} = [];
  [align{nt+1:ntp}] = deal('lc');
end

% get font properties and create commands
if ntp > 0
  [fontsizecmd{1:ntp}]   = deal('');
  [fontanglecmd{1:ntp}]  = deal('');
  [fontweightcmd{1:ntp}] = deal('');
  [colorcmd{1:ntp}]      = deal('');
  [colorclose{1:ntp}]    = deal('');
end
selectfontcmd = '';

if keepfontprops

  % fontsize
  set(htext,'fontunits','points');
  fontsize = get(htext,'fontsize');
  for i=1:ntp
    fontsizecmd{i} = [ '\\fontsize{' num2str(fontsize{i}) '}{' ...
	  num2str(fontsize{i}*1.5) '}'  ];
  end
    
  % fontweight
  fontweight = get(htext,'fontweight');
  for i=1:ntp
    switch fontweight{i}
    case 'light'
      fontweightcmd{i} = [ '\\fontseries{l}\\mathversion{normal}' ];
    case 'normal'
      fontweightcmd{i} = [ '\\fontseries{m}\\mathversion{normal}' ];
    case 'demi'
      fontweightcmd{i} = [ '\\fontseries{sb}\\mathversion{bold}' ];
    case 'bold'
      fontweightcmd{i} = [ '\\fontseries{bx}\\mathversion{bold}' ];
    otherwise
      warning('LaPrint:text',['Unknown fontweight: ' fontweight{i} ])
      fontweightcmd{i} = [ '\\fontseries{m}\\mathversion{normal}' ];
    end
  end  

  % fontangle
  fontangle = get(htext,'fontangle');
  for i=1:ntp
    switch fontangle{i}
    case 'normal'
      fontanglecmd{i} = [ '\\fontshape{n}' ];
    case 'italic'
      fontanglecmd{i} = [ '\\fontshape{it}' ];
    case 'oblique'
      fontanglecmd{i} = [ '\\fontshape{it}' ];
    otherwise
      warning('LaPrint:text',['unknown fontangle: ' fontangle{i} ])
      fontanglecmd{i} = [ '\\fontshape{n}' ];
    end
  end  
  selectfontcmd = '\\selectfont ';
   
end

if color & nt>0 % cgj - was ntp (number of text + legend entries, now just number of text)
  col   = get(htext,'color');
  bgcol = get(htext,'BackgroundColor');
  ecol  = get(htext,'EdgeColor');
  
  for i=1:nt % cgj - was ntp
    col0           = get(get(htext(i),'parent'),'color'); 
    [coli,isc]     = char2rgb(col{i},[0 0 0]);
    [bgcoli,isbgc] = char2rgb(bgcol{i},col0);
    [ecoli,isec]   = char2rgb(ecol{i},col0);
    if isbgc | isec
      set(htext(i),'BackgroundColor','none')
      set(htext(i),'EdgeColor','none')
      colorcmd{i} = ['\\setlength{\\fboxsep}{2pt}\\fcolorbox[rgb]{' ...
        num2str(ecoli(1)) ',' num2str(ecoli(2)) ',' num2str(ecoli(3)) '}{' ...
        num2str(bgcoli(1)) ',' num2str(bgcoli(2)) ',' ...
        num2str(bgcoli(3)) '}{\\color[rgb]{' ...
        num2str(coli(1)) ',' num2str(coli(2)) ',' num2str(coli(3)) '}' ];  
      colorclose{i} = '}';   
    else  
      colorcmd{i} = ['\\color[rgb]{' ...
        num2str(coli(1)) ',' num2str(coli(2)) ',' num2str(coli(3)) '}' ];
    end  
  end  
end

%
% LABELS: modify new figure
%

if ~keepticklabels

  % all axes
  hax = celltoarray(findobj(figno,'type','axes'));
  na  = length(hax);

%   % try to figure out if we have 3D axes an warn
%   issuewarning = 0;
%   for i=1:na
%     issuewarning = max(issuewarning,is3d(hax(i)));
%   end
%   if issuewarning
%     warning('LaPrint:label',['This seems to be a 3D plot. '...
%             'The LaTeX labels are possibly incorrect. '...
%             'The option  ''keepticklabels'' might help. '...
%             'Setting ''figcopy'' to ''off'' might be wise, too.'])
%   end

  % try to figure out if we linear scale with extra factor 
  % and determine powers of 10
  powers = NaN*zeros(na,3);  % matrix with powers of 10 
  for i=1:na                    % all axes
    allxyz = { 'x', 'y', 'z' };
    for ixyz=1:3                % x,y,z
      xyz = allxyz{ixyz};
      ticklabelmode = get(hax(i),[ xyz 'ticklabelmode']);
      if strcmp(ticklabelmode,'auto')
        tick      = get(hax(i),[ xyz 'tick']);
        ticklabel = get(hax(i),[ xyz 'ticklabel']);	      
	    nticklabels    = size(ticklabel,1);
	    nticks    = length(tick);
	    if nticks==0,
          powers(i,ixyz)=0;
          nticklabels=0;
	    end  
	    if nticklabels==0,
          powers(i,ixyz)=0;
  	    end  
        for k=1:nticklabels    % all ticks
	      label = str2num(ticklabel{k});
	      if length(label)==0, 
	        powers(i,ixyz) = 0;
	        break; 
	      end  
	      if ( label==0 ) & ( abs(tick(k))>1e-10 )
	        powers(i,ixyz) = 0;
	        break; 
          end	      
	      if label~=0    
            expon  = log10(tick(k)/label);
	        rexpon = round(expon);
	        if abs(rexpon-expon)>1e-10
              powers(i,ixyz) = 0;
	          break; 
            end	
            if isnan(powers(i,ixyz))
	          powers(i,ixyz) = rexpon;
	        else 	
	          if powers(i,ixyz)~=rexpon
        	    powers(i,ixyz) = 0;
	            break; 
              end		
	        end 
          end  	    
	    end % k	    
      else % if 'auto'
        powers(i,ixyz) = 0;
      end % if 'auto'
    end % ixyz
  end % i
  
  % place text to be replaced by powers on y-axis
  for i=1:na             
    allxyz = { 'x', 'y', 'z' };
    ixyz=2;                % x,y,z
    xyz = allxyz{ixyz};
    leftright=get(hax(i),'yaxislocation');
    if powers(i,ixyz) & ~is3d(hax(i)) & isequal(leftright,'left')
        powertext = ['ypower' int2str(i)];
        xlimit    = get(hax(i),'xlim');
        ylimit    = get(hax(i),'ylim');
        htext     = text(xlimit(1),ylimit(2)+...
                  0.01*(ylimit(2)-ylimit(1)),...
                  powertext);
        set(htext,'VerticalAlignment','Baseline');
    end
  end % i

  % replace all ticklabels and bookkeep
  nxlabel = zeros(1,na);
  nylabel = zeros(1,na);
  nzlabel = zeros(1,na);
  allxyz={ 'x', 'y', 'z' }; 
  for ixyz=1:3
    xyz = allxyz{ixyz};
    k=1;
    if strcmp(xyz,'y') 
      basestr = [ 'v' char(48*ones(1,nzeros-1))];
    else
      basestr = [ xyz char(48*ones(1,nzeros-1))];
    end  
    oldtl  = cell(na,1);
    newtl  = cell(na,1);
    nlabel = zeros(1,na);
    for i=1:na
      % set(hax(i),[ xyz 'tickmode' ],'manual')
      % set(hax(i),[ xyz 'ticklabelmode' ],'manual')
      oldtl{i}  = chartocell(get(hax(i),[ xyz 'ticklabel' ]));
      nlabel(i) = length(oldtl{i});
      newtl{i}  = cell(1,nlabel(i));
      for j=1:nlabel(i)
        newtl{i}{j} = overwritetail(basestr,k);
        k = k+1;
        oldtl{i}{j} = deblank(strrep(strrep(oldtl{i}{j},'\','\\'),...
                             '%','%%'));
      end
      set(hax(i),[ xyz 'ticklabel' ],newtl{i});
    end  
    eval([ 'old' xyz 'tl=oldtl;' ]);
    eval([ 'new' xyz 'tl=newtl;' ]);
    eval([ 'n' xyz 'label=nlabel;' ]);
  end

  % determine latex commands for font properties
  
  if keepfontprops

    % ticklabel font size
    afsize = zeros(na,1);
    for i=1:na
      afsize(i) = get(hax(i),'fontsize');
    end          
    if (any(afsize ~= afsize(1) ))
      warning('LaPrint:text',['Different font sizes for axes not '...
              'supported. All axes will have font size ' ...
	           num2str(afsize(1)) '.' ] )
    end      
    afsizecmd = [ '\\fontsize{' num2str(afsize(1)) '}{' ...
	  num2str(afsize(1)*1.5) '}'  ];

    % ticklabel font weight
    afweight = cell(na,1);
    for i=1:na
      afweight{i} = get(hax(i),'fontweight');
    end
    switch afweight{1}
    case 'light'
      afweightcmd = [ '\\fontseries{l}\\mathversion{normal}' ];
    case 'normal'
      afweightcmd = [ '\\fontseries{m}\\mathversion{normal}' ];
    case 'demi'
      afweightcmd = [ '\\fontseries{sb}\\mathversion{bold}' ];
    case 'bold'
      afweightcmd = [ '\\fontseries{bx}\\mathversion{bold}' ];
    otherwise
      warning('LaPrint:text',['unknown fontweight: ' afweight{1} ])
      afweightcmd = [ '\\fontseries{m}\\mathversion{normal}' ];
    end
    for i=1:na
      if ~strcmp(afweight{i},afweight{1})
        warning('LaPrint:text',['Different font weights for axes '...
                'are not supported. All axes will have font weight ' ...
                afweightcmd '.'])
      end      
    end      

    % ticklabel font angle
    afangle = cell(na,1);
    for i=1:na
      afangle{i} = get(hax(i),'fontangle');
    end
    switch afangle{1}
    case 'normal'
      afanglecmd = [ '\\fontshape{n}' ];
    case 'italic'
      afanglecmd = [ '\\fontshape{it}' ];
    case 'oblique'
      afanglecmd = [ '\\fontshape{it}' ];
    otherwise
      warning('LaPrint:text',['unknown fontangle: ' afangle{1} ])
      afanglecmd=[ '\\fontshape{n}' ];
    end
    for i=1:na
      if ~strcmp(afangle{i},afangle{1})
        warning('LaPrint:text',['Different font angles for axes not '...
                'supported. All axes will have font angle ' ...
                afanglecmd '.'] )
      end      
    end      
  
  end

  % ticklabel color
  acolcmd='';
  if color
    acol=[];
    allxyz={ 'x', 'y', 'z' }; 
    acolwarn = 0;
    for i=1:na
      for ixyz=1:3
        xyzcolor = [allxyz{ixyz} 'color'];
        if  ~isempty(get(hax(i),[allxyz{ixyz} 'ticklabel']))
          if isempty(acol)
            acol = char2rgb(get(hax(i),xyzcolor));  
          else
            if any(char2rgb(get(hax(i),xyzcolor))~=acol)
              acolwarn = 1;
            end 
          end
        end   
      end 
    end
    if acolwarn
      warning('LaPrint:label',['Different colors for axes not ' ...
            'supported. All ticklabels will have color [ ' ...
             num2str(acol) ' ].' ] )
    end
    if ~isempty(acol)
      if any(acol~=[0 0 0])
        acolcmd = [ '\\color[rgb]{' num2str(acol(1)) ',' ...
              num2str(acol(2)) ',' num2str(acol(3)) '}' ];
      end	
    end 
  end

  % ticklabel alignment
    xyzalign = char([116*ones(na,1) 114*ones(na,1) 114*ones(na,1)]); 
    for i=1:na
      switch get(hax(i),'XAxisLocation')
      case 'top'
        xyzalign(i,1)='B';
      end
      switch get(hax(i),'YAxisLocation')
      case 'right'
        xyzalign(i,2)='l';
      end
    end

end

%
% extra picture environment
%

if extrapicture
  unitlength = zeros(na,1);
  ybound     = zeros(na,1);
  for i=na:-1:1   % reverse order, to keep axes in original order
    if ~is3d(hax(i))
      xlim = get(hax(i),'xlim');
      ylim = get(hax(i),'ylim');
      axes(hax(i));
      hori = text(ylim(1),ylim(1),[ 'origin' int2str(i) ]);
      set(hori,'VerticalAlignment','bottom');
      set(hori,'Fontsize',2);
      set(hax(i),'Units','normalized')
      pos = get(hax(i),'Position');
      unitlength(i) = pos(3)*epswidth;
      ybound(i) = (pos(4)*epsheight)/(pos(3)*epswidth);
    else
      warning('LaPrint:extrapic',['Option ''extrapicture'' for 3D ' ...
                  'axes not supported.'])
    end
  end 
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% PART 3 of advanced usage:
%%%% save eps and tex files
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prevent matlab print command to modify lims and ticks 
% (empty, if asonscreen=1)
if ~keepticklabels
  xlimmodeauto       = findobj(figno,'xlimmode','auto');
  xtickmodeauto      = findobj(figno,'xtickmode','auto');
  xticklabelmodeauto = findobj(figno,'xticklabelmode','auto');
  ylimmodeauto       = findobj(figno,'ylimmode','auto');
  ytickmodeauto      = findobj(figno,'ytickmode','auto');
  yticklabelmodeauto = findobj(figno,'yticklabelmode','auto');
  zlimmodeauto       = findobj(figno,'zlimmode','auto');
  ztickmodeauto      = findobj(figno,'ztickmode','auto');
  zticklabelmodeauto = findobj(figno,'zticklabelmode','auto');
  set(xlimmodeauto,'xlimmode','manual')
  set(xtickmodeauto,'xtickmode','manual')
  set(xticklabelmodeauto,'xticklabelmode','manual')
  set(ylimmodeauto,'ylimmode','manual')
  set(ytickmodeauto,'ytickmode','manual')
  set(yticklabelmodeauto,'yticklabelmode','manual')
  set(zlimmodeauto,'ylimmode','manual')
  set(ztickmodeauto,'ytickmode','manual')
  set(zticklabelmodeauto,'yticklabelmode','manual')
end

% create eps file
cmd = strrep(printcmd,'<filename.eps>',epsfullnameext);
cmd = strrep(cmd,'<filename>',filename);
cmd = strrep(cmd,'<figurenumber>',int2str(figno));
  
if verbose
  disp([ 'executing: '' ' cmd ' ''' ]);
end
eval(cmd);

%
% create latex file
%
if verbose
  disp([ 'writing to: '' ' texfullnameext ' ''' ])
end
fid = fopen(texfullnameext,'w');

% head
if head
  fprintf(fid,[ '%% This file is generated by the MATLAB m-file' ...
       ' laprint.m. It can be included\n']);
  fprintf(fid,[ '%% into LaTeX documents using the packages ']);
  fprintf(fid,package);
  if color
      fprintf(fid,', color');
  end    
  fprintf(fid,[ ' and psfrag.\n' ]);
  fprintf(fid,  ['%% It is accompanied by a postscript file. ',... 
     'A sample LaTeX file is:\n']);
  fprintf(fid, '%%    \\documentclass{article}\\usepackage{');
  fprintf(fid,package);
  if color
     fprintf(fid,',color');
  end    
  fprintf(fid, ',psfrag}\n');
  fprintf(fid,[ '%%    \\begin{document}\\input{' ...
	texbasename '}\\end{document}\n' ]);
    fprintf(fid, [ '%% See http://www.mathworks.de/matlabcentral'...
      '/fileexchange/loadFile.do?objectId=4638\n']);
    fprintf(fid, [ '%% for recent versions of laprint.m.\n' ]);
  fprintf(fid,  '%%\n');
  fprintf(fid,[ '%% created by:           ' 'LaPrint version ' ...
	laprintident '\n' ]);
  fprintf(fid,[ '%% created on:           ' datestr(now) '\n' ]);
  fprintf(fid,[ '%% eps bounding box:     ' num2str(epswidth) ...
      ' cm x ' num2str(epsheight) ' cm\n' ]);
  fprintf(fid,[ '%% comment:              ' comment '\n' ]);
  fprintf(fid,'%%\n');
else 
  fprintf(fid,[ '%% generated by laprint.m\n' ]);
  fprintf(fid,'%%\n');
end

% go on
fprintf(fid,'\\begin{psfrags}%%\n');
%fprintf(fid,'\\fontsize{10}{12}\\selectfont%%\n');
fprintf(fid,'\\psfragscanon%%\n');

% text strings

numbertext=0;
for i=1:ntp % cgj - this nt changed to ntp, for figures with no text strings but some extra text strings
            % this never happens in pre-2014b
  numbertext = numbertext+length(newstr{i});
end
if numbertext>0,
  fprintf(fid,'%%\n');
  fprintf(fid,'%% text strings:\n');
  for i=1:ntp
    if length(newstr{i})
      alig = strrep(align{i},'c','');
      fprintf(fid,[ '\\psfrag{' newstr{i} '}[' alig '][' alig ']{' ...
        fontsizecmd{i} fontweightcmd{i} fontanglecmd{i}  ...
        selectfontcmd colorcmd{i} oldstr{i} colorclose{i} '}%%\n' ]);
    end
  end
end

% labels

if ~keepticklabels
  if ~isempty(acolcmd)
     fprintf(fid,'%%\n');
     fprintf(fid,'%% axes ticklabel color:\n');
     fprintf(fid,[ acolcmd '%%\n' ]);
  end    
  if keepfontprops
    fprintf(fid,'%%\n');
    fprintf(fid,'%% axes font properties:\n');
    fprintf(fid,[ afsizecmd afweightcmd '%%\n' ]);
    fprintf(fid,[ afanglecmd '\\selectfont%%\n' ]);
  end  
  nxlabel = zeros(1,na);
  nylabel = zeros(1,na);
  nzlabel = zeros(1,na);
  for i=1:na
    nxlabel(i) = length(newxtl{i});
    nylabel(i) = length(newytl{i});
    nzlabel(i) = length(newztl{i});
  end    
      
  allxyz = { 'x', 'y', 'z' };
  for ixyz=1:3
    xyz = allxyz{ixyz};
    eval([ 'oldtl=old' xyz 'tl;' ]);
    eval([ 'newtl=new' xyz 'tl;' ]);
    eval([ 'nlabel=n' xyz 'label;' ]);
    if sum(nlabel) > 0
      fprintf(fid,'%%\n');
      fprintf(fid,[ '%% ' xyz 'ticklabels:\n']);
      for i=1:na
        poss = ['[' xyzalign(i,ixyz) '][' xyzalign(i,ixyz) ']']; 
        if nlabel(i)
          if strcmp(get(hax(i),[ xyz 'scale']),'linear')
	        % lin scale
            rexpon = powers(i,ixyz);
            if ~rexpon 
              % no powers
              for j=1:nlabel(i)
                fprintf(fid,[ '\\psfrag{' newtl{i}{j} '}' poss '{' ...
		                    Do oldtl{i}{j} Do '}%%\n' ]);
              end 
            else
              % powers
              if ixyz==2 
                leftright=get(hax(i),'yaxislocation');
                if ~is3d(hax(i)) & isequal(leftright,'left')
                  for j=1:nlabel(i)
                    fprintf(fid,[ '\\psfrag{' newtl{i}{j} '}' poss '{' ...
	 	                    Do oldtl{i}{j} Do '}%%\n' ]);
                  end 
                  fprintf(fid,[ '\\psfrag{ypower' int2str(i) ...
                     '}[Bl][Bl]{$\\times 10^{' ...
 		             int2str(rexpon) '}$}%%\n' ]);
                else
                  for j=1:nlabel(i)-1
                    fprintf(fid,[ '\\psfrag{' newtl{i}{j} '}' poss '{' ...
		                      Do oldtl{i}{j} Do '}%%\n' ]);
                  end 
                  if ~is3d(hax(i))
	                fprintf(fid,[ '\\psfrag{' newtl{i}{nlabel(i)} ...
                     '}' poss '{' ... 
                     Do oldtl{i}{nlabel(i)} Do '$\\times 10^{'...
		             int2str(rexpon) '}$}%%\n' ]);
                  else
 	                fprintf(fid,[ '\\psfrag{' newtl{i}{nlabel(i)} ...
                     '}' poss '{\\shortstack{' ... 
                     Do oldtl{i}{nlabel(i)} Do '\\\\$\\times 10^{'...
		             int2str(rexpon) '}\\ $}}%%\n' ]);

                  end 
                end  
              elseif ixyz==1
                for j=1:nlabel(i)-1
                  fprintf(fid,[ '\\psfrag{' newtl{i}{j} '}' poss '{' ...
		                    Do oldtl{i}{j} Do '}%%\n' ]);
                end 
                leftright=get(hax(i),'xaxislocation');
                if isequal(leftright,'bottom')
	              fprintf(fid,[ '\\psfrag{' newtl{i}{nlabel(i)} ...
                     '}' poss '{\\shortstack{' ... 
                     Do oldtl{i}{nlabel(i)} Do '\\\\$\\times 10^{'...
		             int2str(rexpon) '}\\ $}}%%\n' ]);
	            else
                  fprintf(fid,[ '\\psfrag{' newtl{i}{nlabel(i)} ...
                     '}' poss '{\\shortstack{$\\times 10^{' ...
                     int2str(rexpon) '}\\ $\\\\' ...
                     Do oldtl{i}{nlabel(i)} Do '}}%%\n' ]);
                end
              else
                for j=1:nlabel(i)-1
                  fprintf(fid,[ '\\psfrag{' newtl{i}{j} '}' poss '{' ...
		                    Do oldtl{i}{j} Do '}%%\n' ]);
                end 
	            fprintf(fid,[ '\\psfrag{' newtl{i}{nlabel(i)} ...
		           '}' poss '{' Do oldtl{i}{nlabel(i)} Do ...
                   '\\setlength{\\unitlength}{1ex}' ...
		           '\\begin{picture}(0,0)\\put(0.5,1.5){$\\times 10^{' ...
		           int2str(rexpon) '}$}\\end{picture}}%%\n' ]);
              end 
            end % rexpon 
          else
            % log scale
            for j=1:nlabel(i)
              fprintf(fid,[ '\\psfrag{' newtl{i}{j} '}' poss '{${' ...
                oldtl{i}{j} '}$}%%\n' ]);
            end % for (log)
          end % if linear
        end  % if nlabel(i) 
      end
    end
  end
end  

% extra picture
if extrapicture
  fprintf(fid,'%%\n');
  fprintf(fid,'%% extra picture(s):\n');
  for i=1:na
    fprintf(fid,[ '\\psfrag{origin' int2str(i) '}[lb][lb]{' ...
                  '\\setlength{\\unitlength}{' ...
		  num2str(unitlength(i),'%5.5f') 'cm}%%\n' ]);
    fprintf(fid,[ '\\begin{picture}(1,' ...
		  num2str(ybound(i),'%5.5f') ')%%\n' ]);
    %fprintf(fid,'\\put(0,0){}%% lower left corner\n');
    %fprintf(fid,[ '\\put(1,' num2str(ybound(i),'%5.5f') ...
    %	          '){}%% upper right corner\n' ]);
    fprintf(fid,'\\end{picture}%%\n');
    fprintf(fid,'}%%\n');
  end
end  

% figure
fprintf(fid,'%%\n');
fprintf(fid,'%% Figure:\n');
if iscaption
  fprintf(fid,[ '\\parbox{' num2str(latexwidth) 'cm}{\\centering%%\n' ]);
end  
if ~scalefonts
  switch package
  case 'epsfig'   
     fprintf(fid,[ '\\epsfig{file=' epsbasenameext ',width=' ...
	 num2str(latexwidth) 'cm}%%\n' ]);
  case 'graphicx'   
  %   fprintf(fid,[ '\\includegraphics[width=' num2str(latexwidth) ...
    %         'cm]{' epsbasenameext '}%%\n' ]);
     fprintf(fid,[ '\\includegraphics{' epsbasenameext '}%%\n' ]);
  otherwise  
    warning('LaPrint:general',['Package ''' package ''' not known. '...
            'I hope you know what you are doing...'])    
  end
else
  switch package
  case 'epsfig'  
     fprintf(fid,[ '\\resizebox{' num2str(latexwidth) 'cm}{!}' ...
      '{\\epsfig{file=' epsbasenameext '}}%%\n' ]);
  case 'graphicx' 
     fprintf(fid,[ '\\resizebox{' num2str(latexwidth) 'cm}{!}' ...
      '{\\includegraphics{' epsbasenameext '}}%%\n' ]);
  otherwise
    warning('LaPrint:general',['Package ''' package ''' not known. '...
            'I hope you know what you are doing...'])    
  end
end
if iscaption
  fprintf(fid,[ '\\caption{' caption '}%%\n' ]);
  fprintf(fid,[ '\\label{fig:' texbasename '}%%\n' ]);
  fprintf(fid,[ '}%%\n' ]);
end  
fprintf(fid,'\\end{psfrags}%%\n');
fprintf(fid,'%%\n');
fprintf(fid,[ '%% End ' texbasenameext '\n' ]);
fclose(fid);

set(figno,'Name','Printed by LaPrint')
if figcopy
  if verbose
    disp('Strike any key to continue.');
    pause
  end
  hlegend = findobj(figno,'Tag','legend');
  set(hlegend,'DeleteFcn','')
  close(figno)
end

% check for copyobj-bug --> should be ok now!
number_children_new = length(get(figno_ori,'children'));
if number_children_new < number_children_ori
   if figcopy 
      warning(['LaPrint:general','Objects in the figure have been '...
              'deleted! This is due to a bug in matlabs '...
              '''copyopj''. You might want to try to set '...
              '''figcopy'' to ''off''.']) 
   else    
      warning('LaPrint:general',['Objects in the figure have been '...
              'deleted!'])
   end
 end
    
%
% create view file
%

if createview | processview
  if verbose
    disp([ 'writing to: '' ' viewfullnameext ' ''' ])
  end
  fid = fopen(viewfullnameext,'w');

  if head
    fprintf(fid,[ '%% This file is generated by laprint.m.\n' ]);
    fprintf(fid,[ '%% It calls ' texbasenameext ...
		  ', which in turn  calls ' epsbasenameext '.\n' ]);
    fprintf(fid,[ '%% Process this file using, e.g.,\n' ]);
    fprintf(fid,[ '%%   latex ' viewbasenameext '\n' ]);
    fprintf(fid,[ '%%   dvips -o' viewbasename '.ps ' viewbasename ...
            '.dvi\n']);
    fprintf(fid,[ '%%   ghostview ' viewbasename '.ps&\n' ]);
  else 
    fprintf(fid,[ '%% generated by laprint.m\n' ]);
  end

  fprintf(fid,[ '\\documentclass{article}\n' ]);
  fprintf(fid,[ '\\usepackage{' ]);
  fprintf(fid,package);
  if color
      fprintf(fid,',color');
  end    
  fprintf(fid,[ ',psfrag,a4}\n' ]);
  fprintf(fid,[ '\\usepackage[latin1]{inputenc}\n' ]);
  if ~strcmp(epsdirname,viewdirname)
    fprintf(fid,[ '\\graphicspath{{' epsdirname '}}\n' ]);
  end  
  fprintf(fid,[ '\\begin{document}\n' ]);
  fprintf(fid,[ '\\pagestyle{empty}\n' ]);
  if strcmp(texdirname,viewdirname) 
    fprintf(fid,[ '    \\input{' texbasenameext '}\n' ]);
  else
    fprintf(fid,[ '    \\input{' texdirname texbasenameext '}\n' ]);
  end
  fprintf(fid,[ '\\end{document}\n' ]);
  fclose(fid);
end

% process view file
    
if processview
    conti=1;
    
    for i=1:8
      eval(['cmdi=cmd' int2str(i) ';'])
      if ~isempty(cmdi) & conti
        cmd = strrep(cmdi,'<viewfile>',viewbasename);
        cmd = strrep(cmd,'<filename>',filename);
        disp([ 'executing: '' ' cmd ' ''' ]);
        [stat,resu]=system(cmd);
        if stat %| isempty(strfind(resu,'Output written on'))
          disp(resu)
          conti=0;
        end
      end 
    end

    if ~conti
       disp('An error occured in the latex/friends sequence.')
    end

end

set(0,'ShowHiddenHandles',shh);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% functions used
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fullnameext,basenameext,basename,dirname] = getfilenames(...
    filename,extension,verbose)
% appends an extension to a filename (as '/home/tom/tt') and determines  
% fullnameext: filename with extension and dirname, as '/home/tom/tt.tex'
% basenameext: filename with extension without dirname, as 'tt.tex'
% basename   : filename without extension without dirname, as 'tt'
% dirname    : dirname without filename, as '/home/tom/'
% In verbose mode, it asks if to overwrite or to modify.
%
[dirname, basename] = splitfilename(filename);
fullnameext = [ dirname basename '.' extension ];
basenameext = [ basename '.' extension ];
if verbose
  quest = (exist(fullnameext)==2);
  while quest
    yn = input([ strrep(strrep(fullnameext,'\','\\'),'%','%%') ...
            ' exists. Overwrite? (y/n) '],'s');
    if strcmp(yn,'y') 
      quest = 0;
    else
      filename = input( ...
	             [ 'Please enter new filename (without extension .' ...
	             extension '): ' ],'s');
      [dirname, basename] = splitfilename(filename);
      fullnameext = [ dirname basename '.' extension ];
      basenameext = [ basename '.' extension ];
      quest = (exist(fullnameext)==2);
    end
  end
end
if ( exist(dirname)~=7 & ~strcmp(dirname,[ '.' filesep ]) ...
      & ~strcmp(dirname,filesep) )
  error([ 'Directory ' dirname ' does not exist.' ] )
end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dirname,basename] = splitfilename(filename)
% splits filename into dir and base
slashpos  = findstr(filename,filesep);
nslash    = length(slashpos);
nfilename = length(filename);
if nslash
  dirname  = filename(1:slashpos(nslash));
  basename = filename(slashpos(nslash)+1:nfilename);
else
  dirname = pwd;
  nn=length(dirname);
  if ~strcmp(dirname(nn),filesep)
    dirname = [ dirname filesep ];
  end   
  basename = filename;
end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yesno = is3d(haxes)
% tries to figure out if axes is 3D
yesno = 0;
CameraPosition = get(haxes,'CameraPosition');
CameraTarget = get(haxes,'CameraTarget');
CameraUpVector = get(haxes,'CameraUpVector');
if CameraPosition(1)~=CameraTarget(1)
  yesno = 1;
end  
if CameraPosition(2)~=CameraTarget(2)
  yesno = 1;
end  
if any(CameraUpVector~=[0 1 0])
  yesno = 1;
end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = celltoarray(a)
% converts a cell of doubles to an array
if iscell(a),
  b = [];
  for i=1:length(a),
    b = [b a{i}]; 
  end  
else
  b = a(:)';
end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = chartocell(a)
% converts a character array into a cell array of characters

% convert to cell 
if isa(a,'char')
  n = size(a,1);
  b = cell(1,n);
  for j=1:n
    b{j}=a(j,:); 
  end  
else
  b = a;
end  
% convert to char
n=length(b);
for j=1:n
  if isa(b{j},'double')
    b{j} = num2str(b{j});
  end  
end	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = overwritetail(a,k)
% overwrites tail of a by k
% a,b: strings
% k: integer
ks = int2str(k);
b = [ a(1:(length(a)-length(ks))) ks ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rgb,isc] = char2rgb(c,c0)
% convert color definitions from character to rgb-vector

isc = 1;
if ~ischar(c)
  rgb = c; 
else  
  switch c
  case {'y','yellow'}
    rgb = [1 1 0];
  case {'m','magenta'}
    rgb = [1 0 1];
  case {'c','cyan'}
    rgb = [0 1 1];
  case {'r','red'}
    rgb = [1 0 0];
  case {'g','green'}
    rgb = [0 1 0];
  case {'b','blue'}
    rgb = [0 0 1];
  case {'w','white'}
    rgb = [1 1 1];
  case {'k','black'}
    rgb = [0 0 0];
  case 'none' 
    if nargin==2
      rgb = char2rgb(c0);
    else
      rgb = [1 1 1];    
    end    
    isc = 0;
  otherwise
    warning('LaPrint:general',['Unknown Color: ''' c ...
            '''. Taking black.'])
    rgb = [0 0 0];
  end    
end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = value01(val,opt)
% convert off/on to 0/1

if nargin==2
   txt = ['Value of ' opt ' must be ''on'' or ''off'''];
else
   txt = ['Value must be ''on'' or ''off'''];
end

if  ~isa(val,'char')  
  error(txt)
end
val = lower(strrep(val,' ',''));
switch val
  case 'on'
    v = 1;
  case 'off'
    v = 0;
  otherwise
    error(txt)
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function opt = prefsettings
    if ispref('LaPrint','LAPRINTOPT')  
      opt = getpref('LaPrint','LAPRINTOPT');
    else
      opt = [];  
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function opt=factorysettings

% try to find LaTeX and friends 
if ispc
  try 
    latexpath = winqueryreg('HKEY_LOCAL_MACHINE',...
	  'SOFTWARE\MiK\MiKTeX\CurrentVersion\MiKTeX','Install Root');
    latexcmd = [latexpath '\miktex\bin\latex.exe -halt-on-error '...
               '-interaction nonstopmode <viewfile>.tex'];
    dvipscmd = [latexpath '\miktex\bin\dvips.exe -D600 -E* '...
               '-o<viewfile>.eps <viewfile>.dvi'];
  catch    % hoping the path variable is properly set
    latexcmd = ['latex.exe -halt-on-error '...
               '-interaction nonstopmode <viewfile>.tex'];
    dvipscmd = ['dvips.exe -D600 -E* '...
               '-o<viewfile>.eps <viewfile>.dvi'];
  end
  epstoolcmd = ['C:\Ghostgum\epstool\bin\epstool.exe '...
               '--bbox --copy --output '...
               '<filename>_final.eps <viewfile>.eps'];
  delcmd =     ['del <viewfile>.eps <viewfile>.dvi ',...
               '<viewfile>.aux <viewfile>.log ',...
               '<viewfile>.pfg'];
  gsviewcmd = 'C:\Ghostgum\gsview\gsview32.exe <filename>_final.eps&';
else % hoping the path variable is properly set
  latexcmd =   ['latex -halt-on-error '...
               '-interaction nonstopmode <viewfile>.tex'];
  dvipscmd =   ['dvips -D600 -E* '...
               '-o<viewfile>.eps <viewfile>.dvi'];
  epstoolcmd = ['epstool --bbox --copy --output '...
               '<filename>_final.eps <viewfile>.eps'];
  delcmd =     ['rm <viewfile>.eps <viewfile>.dvi ',...
               '<viewfile>.aux <viewfile>.log ',...
               '<viewfile>.pfg'];
  gsviewcmd =  'ghostview <filename>_final.eps&';
end

vers = version;
vers = eval(vers(1:3));
if vers < 6.5
  colorvalue=0;
else
  colorvalue=1;
end

   opt = struct(...
      'figno',{1},...
      'filename','unnamed',...
      'width',12,...
      'factor',0.8,...
      'scalefonts',1,...
      'keepfontprops',0,...
      'asonscreen',0,...
      'keepticklabels',0,...
      'mathticklabels',0,...
      'head',1,...
      'comment','',...
      'caption','',...
      'extrapicture',0,...
      'nzeros',3,...
      'verbose',0,...
      'figcopy',1,...
      'printcmd',['print(''-f<figurenumber>'',' ...
                      '''-depsc'',''-painters'','...
                      '''<filename.eps>'')'],...
      'package','graphicx',...
      'color',colorvalue,...
      'createview',0,...
      'viewfilename','unnamed_',...
      'processview',0,...
      'cmd1',latexcmd,...
      'cmd2',dvipscmd,...
      'cmd3',epstoolcmd,...
      'cmd4',delcmd,...
      'cmd5',gsviewcmd,...
      'cmd6','',...
      'cmd7','',...
      'cmd8','');
end
