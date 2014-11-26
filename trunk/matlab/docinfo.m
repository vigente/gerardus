function docinfo
% This is a placeholder matlab function to provide info about creating
% documentation in gerardus.
%
% HTML files created with the matlab publisher go in the /doc folder
% For them to appear as text/html mimetype on a server such as googlecode, the
% appropriate svn mimetypes need to be set. An example is shown below in this
% function's file.
%
% also see:
% http://manjeetdahiya.com/2010/09/29/serving-html-documentation-from-google-code-svn/
% https://code.google.com/p/flexlib/wiki/HowToBuild
% http://uk.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html
%
% Click <a href="matlab:web('../doc/docinfo/doc.html')">here</a> for Documentation with examples: 


%%% The following is example lines from a configuration file as found in
%%% ~/.subversion.config on posix systems, to set appropriate mimetypes for
%%% common file-extensions.

% ### Section for configuring miscelleneous Subversion options.
% [miscellany]
% ### Set enable-auto-props to 'yes' to enable automatic properties
% ### for 'svn add' and 'svn import', it defaults to 'no'.
% ### Automatic properties are defined in the section 'auto-props'.
% enable-auto-props = yes
% 
% ### Section for configuring automatic properties.
% [auto-props]
% ### The format of the entries is:
% ###   file-name-pattern = propname[=value][;propname[=value]...]
% ### The file-name-pattern can contain wildcards (such as '*' and
% ### '?').  All entries which match (case-insensitively) will be
% ### applied to the file.  Note that auto-props functionality
% ### must be enabled, which is typically done by setting the
% ### 'enable-auto-props' option.
% 
% ### Mime-types for common file types (adapted from: https://code.google.com/p/flexlib/wiki/HowToBuild)
% ### Also see: http://manjeetdahiya.com/2010/09/29/serving-html-documentation-from-google-code-svn/
% 
% # Scriptish formats
% *.bat        = svn:eol-style=native; svn:keywords="Rev Date Author"; svn-mine-type=text/plain
% *.bsh        = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/x-beanshell
% *.cgi        = svn:eol-style=native; svn:keywords="Rev Date Author"; svn-mine-type=text/plain
% *.cmd        = svn:eol-style=native; svn:keywords="Rev Date Author"; svn-mine-type=text/plain
% *.js         = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/javascript
% *.php        = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/x-php
% *.pl         = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/x-perl; svn:executable
% *.pm         = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/x-perl
% *.py         = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/x-python; svn:executable
% *.sh         = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/x-sh; svn:executable
% 
% # Image formats
% *.bmp        = svn:mime-type=image/bmp
% *.gif        = svn:mime-type=image/gif
% *.ico        = svn:mime-type=image/ico
% *.jpeg       = svn:mime-type=image/jpeg
% *.jpg        = svn:mime-type=image/jpeg
% *.png        = svn:mime-type=image/png
% *.tif        = svn:mime-type=image/tiff
% *.tiff       = svn:mime-type=image/tiff
% 
% # Data formats
% *.pdf        = svn:mime-type=application/pdf
% *.avi        = svn:mime-type=video/avi
% *.eps        = svn:mime-type=application/postscript
% *.gz         = svn:mime-type=application/gzip
% *.mov        = svn:mime-type=video/quicktime
% *.mp3        = svn:mime-type=audio/mpeg
% *.ps         = svn:mime-type=application/postscript
% *.psd        = svn:mime-type=application/photoshop
% *.rtf        = svn:mime-type=text/rtf
% *.swf        = svn:mime-type=application/x-shockwave-flash
% *.tgz        = svn:mime-type=application/gzip
% *.wav        = svn:mime-type=audio/wav
% *.zip        = svn:mime-type=application/zip
% 
% # Microsoft-specific formats (see: http://filext.com/faq/office_mime_types.php)
% *.doc        = svn:mime-type=application/msword
% *.docx       = svn:mime-type=application/vnd.openxmlformats-officedocument.wordprocessingml.document
% *.xls        = svn:mime-type=application/vnd.ms-excel
% *.xlsx       = svn:mime-type=application/vnd.openxmlformats-officedocument.spreadsheetml.sheet
% *.ppt        = svn:mime-type=application/vnd.ms-powerpoint
% *.pptx       = svn:mime-type=application/vnd.openxmlformats-officedocument.presentationml.presentation
% 
% # Matlab-specific binary formats (adapted from: http://www.mathworks.co.uk/help/matlab/matlab_prog/set-up-svn-source-control.html)
% *.mdl        = svn:mime-type=application/octet-stream
% *.mat        = svn:mime-type=application/octet-stream 
% *.slx        = svn:mime-type= application/octet-stream
% *.mdlp       = svn:mime-type=application/octet-stream
% *.slxp       = svn:mime-type=application/octet-stream
% *.sldd       = svn:mime-type=application/octet-stream
% *.p          = svn:mime-type=application/octet-stream
% *.mexa64     = svn:mime-type=application/octet-stream
% *.mexw32     = svn:mime-type=application/octet-stream
% *.mexw64     = svn:mime-type=application/octet-stream
% *.mexmaci64  = svn:mime-type=application/octet-stream
% 
% # Text formats
% .htaccess    = svn:mime-type=text/plain
% *.css        = svn:mime-type=text/css
% *.dtd        = svn:mime-type=text/xml
% *.html       = svn:mime-type=text/html
% *.ini        = svn:mime-type=text/plain
% *.sql        = svn:mime-type=text/x-sql
% *.txt        = svn:mime-type=text/plain
% *.xhtml      = svn:mime-type=text/xhtml+xml
% *.xml        = svn:mime-type=text/xml
% *.xsd        = svn:mime-type=text/xml
% *.xsl        = svn:mime-type=text/xml
% *.xslt       = svn:mime-type=text/xml
% *.xul        = svn:mime-type=text/xul
% *.yml        = svn:mime-type=text/plain
% CHANGES      = svn:mime-type=text/plain
% COPYING      = svn:mime-type=text/plain
% INSTALL      = svn:mime-type=text/plain
% Makefile*    = svn:mime-type=text/plain
% README       = svn:mime-type=text/plain
% TODO         = svn:mime-type=text/plain
% 
% # Code formats
% *.c          = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/plain
% *.cpp        = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/plain
% *.h          = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/plain
% *.java       = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/plain
% *.as         = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/plain
% *.mxml       = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/plain
% *.m          = svn:eol-style=native; svn:keywords="Rev Date Author"; svn:mime-type=text/plain