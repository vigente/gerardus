

# Contributing to the project #

  1. Any communications should be carried on the [mailing list](http://groups.google.com/group/gerardus-users/) (gerardus-users@googlegroups.com)
  1. There are two ways to contribute code:
    * Sending patches to the mailing list. Then one of the developers will apply the patch, test it, and if there are no problems, commit the changes to the repository
    * If you are in the [committers list](http://code.google.com/p/gerardus/people/list), then you can commit your code directly to the repository

## Adding new Matlab funtions to the repository ##

  1. Add the new function to the desired folder in the repository
  1. Add a test script to the directory "gerardus\matlab\test"
  1. Open the SVN properties window and import "common.svnprops" from the gerardus directory - do this for both files
  1. Update the ChangeLog with details of the commit
  1. Select "SVN Commit" of the gerardus folder

## Sending patches ##

If you spot a bug in the code, or want to improve a function, but you are not in the [committers list](http://code.google.com/p/gerardus/people/list), this is the simplest way to contribute.

  1. Check out a read-only version of Gerardus
    * In Linux
```
$ svn co http://gerardus.googlecode.com/svn/trunk/ gerardus
```
    * In Windows (e.g. using TortoiseSVN) select the URL "`http://gerardus.googlecode.com/svn/trunk/`" and the destination folder "`gerardus`"
  1. Or if you already have checked out the code in the past, make sure that you have the latest version
    * In Linux
```
$ svn up
```
    * In Windows, select right-click and select "Update" from the TortoiseSVN menu
  1. Edit the files you want to improve
  1. Create a patch with your changes
    * In Linux
```
$ cd gerardus
$ svn diff > file.patch
```
    * In Windows, right-click on the "`gerardus`" folder and select "Diff" from the TortoiseSVN menu. Save the result as e.g. `file.patch`
  1. Email `file.patch` to the mailing list with an explanation of what you have done

# Notes about reading vectors from Matlab into different types of C++ vector-like objects #

(This section was spawned from this [email sent to the gerardus-user mailing list](https://groups.google.com/d/msg/gerardus-users/pLH0iR0H74o/rNC2Qr20_iMJ))

Our `MatlabImportFilter` can read vectors as input parameters directly from Matlab, and pass them to ITK without intermediate steps. For example, in the `itk::VotingBinaryIterativeHoleFillingImageFilter`

https://code.google.com/p/gerardus/source/browse/trunk/matlab/ItkToolbox/ItkImFilter.cpp?r=1113#506

```
// The variance for the discrete Gaussian kernel. Sets the
// variance independently for each dimension. The default is 0.0
// in each dimension (ITK)
typename FilterType::ArrayType defVariance;
defVariance.Fill(0.0);
filter->SetVariance(matlabImport->
                     GetRowVectorArgument<typename FilterType::ArrayType::ValueType, 
                                             typename FilterType::ArrayType,
                                             VImageDimension>(2, "VAR", defVariance));
```

The method `GetRowVectorArgument()` is declared here

https://code.google.com/p/gerardus/source/browse/trunk/matlab/MatlabImportFilter.h?r=1113#194

Because we may want to read Matlab vectors into very different classes, e.g. `itk::Size<Dimension>` vs. `CGAL::Direction_3<CGAL::Simple_cartesian<double> >`, the definition of this method relies heaviy on the idea of a vector wrapper,

https://code.google.com/p/gerardus/source/browse/trunk/matlab/MatlabImportFilter.hxx?r=1113#326

```
VectorWrapper<VectorValueType, VectorType, mxLogical, VectorSize> paramWrap;
```

That is, we have a class `VectorWrapper` that provides a common interface from Matlab to different sorts of vectors. For the moment, by default we assume that the output vector is

https://code.google.com/p/gerardus/source/browse/trunk/matlab/VectorWrapper.h?r=1113#75

```
std::vector
```

but we have partial specialisations to deal with

https://code.google.com/p/gerardus/source/browse/trunk/matlab/VectorWrapper.h?r=1113#121

```
itk::Size<Dimension>
itk::FixedArray<ValueType, Dimension>
```

https://code.google.com/p/gerardus/source/browse/trunk/matlab/VectorWrapper.h?r=1113#205

```
CGAL::Point_3<CGAL::Simple_cartesian<ValueType> >
CGAL::Direction_3<CGAL::Simple_cartesian<ValueType> >
```

For ITK this already covers quite a bit of possibilities, because many types get reduced to the two above. For example, the `itk::BoxFilterType<...,...>::RadiusType` used in the Median filter is actually an `itk::FixedArray<double, Dimension>`.

https://code.google.com/p/gerardus/source/browse/trunk/matlab/ItkToolbox/ItkImFilter.cpp?r=1113#752


# Matlab code guidelines #

  * All Matlab (.m) and MEX (.cpp, .hpp) functions live in directory [`gerardus\matlab`](http://code.google.com/p/gerardus/source/browse/#svn%2Ftrunk%2Fmatlab), under 6 toolbox directories:
    * [CardiacToolbox](http://code.google.com/p/gerardus/source/browse/#svn%2Ftrunk%2Fmatlab%2FCardiacToolbox): functions that are only relevant for heart modelling and image processing
    * [FileFormatToolbox](http://code.google.com/p/gerardus/source/browse/#svn%2Ftrunk%2Fmatlab%2FFileFormatToolbox): functions to load and save data from and to files, and convert between file formats
    * [FiltersToolbox](http://code.google.com/p/gerardus/source/browse/#svn%2Ftrunk%2Fmatlab%2FFiltersToolbox): functions to run filters or general image processing on 2D/3D images
    * [ItkToolbox](http://code.google.com/p/gerardus/source/browse/#svn%2Ftrunk%2Fmatlab%2FItkToolbox): MEX functions to run ITK filters and transforms from Matlab
    * [PointsToolbox](http://code.google.com/p/gerardus/source/browse/#svn%2Ftrunk%2Fmatlab%2FPointsToolbox): functions to process sets of points, landmarks, configurations...
    * [ThirdParty](http://code.google.com/p/gerardus/source/browse/#svn%2Ftrunk%2Fmatlab%2FThirdPartyToolbox): functions written by people from other projects, but that Gerardus uses
  * All function names are lower-case, e.g. `scinrrd_intersect_plane()`
  * All variable names, input and output arguments are lower-case. An exception to this is when the variable or argument is a boolean flag, when sometimes it's OK to make it all upper-case, e.g. `FOUND = true;`
  * All functions must have a help header. Have a look at the current functions to get an idea of the format, e.g.
```
% SCINRRD_INTERSECT_PLANE  Intersection of a plane with an image volume
%
% [IM, GX, GY, GZ, MIDX] = SCINRRD_INTERSECT_PLANE(NRRD, M, V)
%
%   IM is an image that displays the intersection of the plane with the
%   image volume in SCI NRRD format. Voxels that fall outside the image
%   volume are returned as NaN.
%
%   GX, GY, GZ are matrices of the same size as IM, and contain the
%   Cartesian coordinates of the voxels in IM. You can visualize the
%   resulting plane using
%
%     >> surf(gx, gy, gz, im, 'EdgeColor', 'none')
[...]
```
  * All non-third party functions must have a copyright notice, that when first created will look something like this
```
% Authors: John Doe <johndoe@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.0
% $Rev$
% $Date$
% 
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
```
  * If it's a new function, you need to add to the file the subversion property "`svn:eol-style=native`" (so that end-of-line characters will be correct in any operating system), and the keywords "`Rev Date`" (so that the help header gets automatically updated with the revision number and date after each commit)
    * In Linux
```
$ svn propset svn:keywords "Rev Date" newfile.m
$ svn propset svn:eol-style native newfile.m
```
  * The code must have comments. Start comments with lower-case and no period at the end, e.g.
```
% keep track of where the rotation point is using a boolean vector for the
% row position, and another one for the column position. The rotation point
% is at the center of the grid
```
  * The number of input/output arguments is checked in the Matlab file with
```
% check arguments
error(nargchk(3, 4, nargin, 'struct'));
error(nargoutchk(0, 5, nargout, 'struct'));
```
  * Values for missing or empty arguments are set using, e.g.
```
% defaults
if (nargin < 4 || isempty(interp))
    interp = 'nn';
end
```
  * No white spaces before/after brackets, e.g.
```
lmax = ceil(sqrt(sum(([nrrd.axis.size] - 1).^2)));
```
  * White space after commas and operations, e.g.
```
[gc, gr] = meshgrid(-lmax:lmax, -lmax:lmax);
grcs(1, :) = grcs(1, :) + idxm(1);
```

# Integrating the Boost C++ libraries with Gerardus #

These are notes of what we did to get the Boost libraries code into Gerardus.

  1. Go to the [Boost C++ Libraries](http://www.boost.org/) and click on the current release link (1.53.0 at the time of this writing)
  1. Download the tarball with the whole project, e.g. [boost\_1\_53\_0.tar.bz2](https://sourceforge.net/projects/boost/files/boost/1.53.0/boost_1_53_0.tar.bz2/download)
  1. Extract `boost_1_53_0` to `gerardus/cpp/src/third-party`. This takes up 456M, which is too much to distribute with Gerardus
  1. Change to the boost source code directory
```
cd gerardus/cpp/src/third-party/boost_1_53_0
```
  1. Delete documentation, tests and examples to reduce the size to 144M
```
find . -name test | xargs rm -rf
find . -name doc | xargs rm -rf
find . -name examples | xargs rm -rf
find . -name example | xargs rm -rf
```
  1. Add and commit the boost source code to the gerardus repository
```
cd ..
svn add boost_1_53_0
svn ci boost_1_53_0
```

# Debugging problems with MEX functions in Windows #

MEX files in Windows are DLLs that can link to other DLLs. Sometimes a MEX file may not run giving the error
```
Invalid MEX-file <mexfilename>:
The specified module could not be found.
```
As explained in Mathworks support page ["Invalid MEX-File Error"](http://www.mathworks.co.uk/help/matlab/matlab_external/invalid-mex-file-error.html), in Windows one can install the program Dependency Walker, and run it within Matlab to figure out whether all DLLs the MEX file depends on are visible, e.g.
```
>> !"C:\Program Files\depends22_x64\depends.exe" ..\ItkToolbox\itk_imfilter.mexw64
```
To make a DLL visible, its has to be in the same directory as the MEX file, or it's directory has to be in the system path (which is different from the Matlab toolbox path). In Gerardus, we solve this problem in the script [add\_gerardus\_paths.m](https://code.google.com/p/gerardus/source/browse/trunk/matlab/add_gerardus_paths.m) with the code
```
%% System paths to DLLs for Windows
if ~isempty(strfind(getenv('OS'), 'Windows'))
    
    % full paths to directories with DLLs
    pathsToDlls = {cd(cd('..\lib')), cd(cd('..\lib\bin')), cd(cd('..\cpp\src\third-party\CGAL-4.2\auxiliary\gmp\lib'))};
    
    % get system paths
    systemPath = getenv('PATH');
    
    % add paths to directories with DLLs unless they are already in the
    % system path
    for I = 1:length(pathsToDlls)
        if isempty(strfind(lower(systemPath), lower(pathsToDlls{I})))
            disp(['Adding ' pathsToDlls{I} ' to system path'])
            systemPath = [pathsToDlls{I} ';' systemPath];
        end
    end
    setenv('PATH', systemPath);
end
```
This code searches the system path, and adds any DLL directory that is missing.

# Usage of TortoiseSVN in Windows #

  1. [Download and install TortoiseSVN](http://tortoisesvn.net/downloads.html) for Windows. Reboot Windows after installation is successful.
  1. Check out Gerardus by right-clicking inside a folder where you want to have it, and selecting "SVN Checkout" in the contextual menu with the following parameters:
    * URL of repository: https://gerardus.googlecode.com/svn/trunk
    * Checkout directory: C:\path\to\folder\where\you\want\gerardus
    * Checkout Depth: Fully recursive.
    * Revision: HEAD revision.
  1. Make your changes to the repository.
  1. Select the files you want to commit. Right click and choose "SVN Commit..." on the contextual menu. If this is the first time you try to commit to Gerardus, a dialogue will pop up asking for your usename/password:
```
<https://gerardus.googlecode.com:443> Google Code Subversion Repository

Requests a username and password
Username:
Password:
□ Save authentication
```
    * Your username is your Gmail address (e.g. john.doe@gmail.com).
    * Go to the Gerardus Googlecode website to obtain your [Gerardus password](https://code.google.com/hosting/settings). This password is different from your Gmail password, and is automatically generated by Google Code.
    * Tick the "Save authentication" box so that in future your commits won't require entering your username/password again.