

# Quickstart #

If you just need pure Matlab functions, it's as easy as downloading the code with Subversion (instructions [here](https://code.google.com/p/gerardus/source/checkout)) and then running from Matlab the script `gerardus/matlab/add_gerardus_paths.m` to add the Gerardus toolboxes to your Matlab path.

If you want to use also the Matlab MEX functions and C++ programs, then you need to also build and/or install some dependencies and run Gerardus' CMake build process (instructions below).

[Gerardus release 0.9.0](https://code.google.com/p/gerardus/source/browse/?r=1176#svn%2Ftags%2Frelease-0.9.0) builds on Linux 64 bit and Windows 64 bit (Visual Studio 2010). Mac OSX does not build out-of-the-box, although Tom Doel is looking into it.

# Build and/or install dependencies #

## Summary ##

  * Subversion client
  * g++ compiler
  * CMake (currently tested with v2.8)
  * m4
  * OpenGL library
  * Insight Toolkit v4.3.1 (Windows only)
  * Matlab (currently tested with R2012a, R2012b, R2013b, R2014a)

  * One-liner for linux that installs everything except ITK and Matlab. ITK will be downloaded and built by the Gerardus CMakeLists.txt script, and Matlab has to be installed separately
```
sudo apt-get install subversion cmake g++ g++-4.7 m4 libgl1-mesa-dev
```

## In detail ##

  * [Subversion client](http://subversion.apache.org/)
    1. Download and install the subversion package appropriate for your operating system
    * Note for Linux users: Subversion is probably packaged for your distribution. For example, Ubuntu users can install it doing
```
$ sudo apt-get install subversion
```
    * Note for Windows users: A popular choice for Windows users is the [TortoiseSVN client](http://tortoisesvn.net/downloads.html)
    * Note for Mac users: You have a choice of installing a graphical front-end such as [Svnx](http://code.google.com/p/svnx/) or using subversion from the command-line. The command-line tools are installed with Xcode and can be found in the following directory:
```
/Applications/Xcode.app/Contents/Developer/usr/bin
```
    * Note for other users: See the [subversion package download page](http://subversion.apache.org/packages.html) for more options
  * g++ Compiler
    * Linux users must install both g++ and the version required by Matlab, e.g. g++-4.4 for Matlab2012a/Matlab2012b or g++-4.7 for Matlab2013b. E.g. Ubuntu/Debian users can install it with
```
$ sudo apt-get install g++ g++-4.4

or

$ sudo apt-get install g++ g++-4.7
```
    * Mac users can install [Xcode](https://developer.apple.com/xcode/) free from Apple's App Store. This includes gcc and Subversion. After installing Xcode you need to install the command-line tools - to do this open the preferences window in Xcode (start Xcode, open the Xcode menu and choose Preferences). Then click the Downloads tab and click the install button next to the `Command Line Tools` in the Components list.
    * For windows users, [Visual C++ Express](http://www.microsoft.com/visualstudio/en-us/products/2010-editions/visual-cpp-express) is available free from Microsoft and includes C++ compilers.
  * [CMake](http://www.cmake.org/) (Gerardus requires version >= 2.8.8)
    1. Download and install CMake for your system (Linux, Windows, MacOS X, etc).
    * Note for Linux users: CMake is probably packaged for your distribution. Ubuntu packages CMake 2.8.6 in Precise Pangolin (12.04), so it's necessary to install cmake 2.8.9 from the upcoming Quantal Quetzal (12.10). This can be done by replacing `precise` by `quantal` in `/etc/apt/sources.list`, then running
```
$ sudo apt-get update
$ sudo apt-get install cmake
```
> > > Then undo the edit and run again
```
$ sudo apt-get update
```
    * Note for Windows users: Download and run the Win32 Installer from the [CMake downloads page](http://www.cmake.org/cmake/resources/software.html)
    * Note for Mac users: Download and run the Mac OSX 64/32-bit Universal Installer from the [CMake downloads page](http://www.cmake.org/cmake/resources/software.html)
    * Note for other users: See the [CMake resources page](http://www.cmake.org/cmake/resources/software.html) for more options
  * [m4](http://www.gnu.org/software/m4/) (GNU M4 is an implementation of the traditional Unix macro processor)
    1. It's a small program, simply download and install it
      * Note for Linux users: In Ubuntu, it can be easily installed doing
```
$ sudo apt-get install m4
```
  * [OpenGL library](https://www.opengl.org/)
    * Gerardus has this dependency because `MatlabImportFilter` uses the `_image` type from `CGAL_ImageIO`, which depends on OpenGL to build.
    * Linux users: You need the development package of the OpenGL library. E.g. in Ubuntu
```
$ sudo apt-get install libgl1-mesa-dev
```
  * [Insight Toolkit](http://www.itk.org/itkindex.html) (Gerardus tested with version 4.3.1)
    * Linux users: Starting from `gerardus/CMakeLists.txt` v0.9.0, ITK is automatically downloaded and built as part of the Gerardus building process, so it doesn't need to be installed separately. This may also work for MacOS X users, but I have not tested it.
    * Windows users: Download the [[Insight Toolkit (ITK) v4.3.1 source code](http://www.itk.org/ITK/resources/software.html). Follow the instructions to build and install ITK using CMake (see the [ITK Software Guide](http://www.itk.org/ItkSoftwareGuide.pdf) for details)
    * MacOS X users: In case the Gerardus build system fails to automatically download and build ITK, then [this site](http://worldwidepenguin.com/2010/05/how-to-install-itk-on-a-mac-cmake-macports/) has instructions for building ITK separately (although they forget to mention that you click "Generate" on the CMake GUI at the end)
    * If you build ITK separately, you need to set the following variables in your CMake build of ITK
```
CMAKE_BUILD_TYPE   Release (this one ignored by Windows)
BUILD_SHARED_LIBS  ON
ITK_USE_REVIEW     ON
ITK_LEGACY_REMOVE  ON
```
      * Windows users can do this checking the "Advanced" box in the CMake interface, and then checking the `ITK_USE_REVIEW` and `ITK_LEGACY_REMOVE`. If using Visual Studio, then the `CMAKE_BUILD_TYPE` entry will be ignored anyway, because Visual Studio is a multi-configuration sytem
    * Note for everyone: You can also untick the following variables to save time, as Gerardus does not need the ITK examples, documentation or tests
```
BUILD_DOCUMENTATION   OFF
BUILD_EXAMPLES        OFF
BUILD_TESTING         OFF
```
    * Note for Linux and MacOSX: ITK needs to be compiled with the same version of the compiler used for Matlab. This can be ensured setting the following CMake variables (e.g. [Matlab2012a/Matlab2012b/Matlab2013a for Linux require gcc 4.4](http://www.mathworks.com/support/compilers/R2012a), and [Matlab2013b/Matlab2014a](http://www.mathworks.com/support/compilers/R2013b/?sec=glnxa64) require gcc 4.7)
```
CMAKE_C_COMPILER gcc-4.4
CMAKE_CXX_COMPILER g++-4.4

or

CMAKE_C_COMPILER gcc-4.7
CMAKE_CXX_COMPILER g++-4.7
```
    * Linux and MacOSX users can save time running cmake like this:
```
$ cmake -DBUILD_SHARED_LIBS=ON -DITK_USE_REVIEW=ON \
-DCMAKE_BUILD_TYPE=Release -DITK_LEGACY_REMOVE=ON \
-DCMAKE_C_COMPILER=/usr/bin/gcc-4.4 \
-DCMAKE_CXX_COMPILER=/usr/bin/g++-4.4  ..

or

$ cmake -DBUILD_SHARED_LIBS=ON -DITK_USE_REVIEW=ON \
-DCMAKE_BUILD_TYPE=Release -DITK_LEGACY_REMOVE=ON \
-DCMAKE_C_COMPILER=/usr/bin/gcc-4.7 \
-DCMAKE_CXX_COMPILER=/usr/bin/g++-4.7  ..
```
    * Notes for Windows users:
      * Start Visual Studio as **Administrator**. This is required to install ITK
      * The ITK folder must be close to the root directory C:/>, because sometimes ITK gives compilation problems in Windows if the paths to files are very long strings
      * It's important to select "Release" instead of "Debug" from the pull-down configuration menu above the Solution Explorer. This will enable to build and install the Release version of ITK, which will then be needed by the Release build of Gerardus
      * Click on "Start" -> "My computer" -> "View system information", and in the "Advanced" tab click on "Environment Variables". Select the System variable `Path`, and click "Edit". Add the path to the `ITKCOMMON.DLL` library, e.g. "C:\Program Files\ITK\bin".
  * [Matlab](http://www.mathworks.com/products/matlab/)
    1. Matlab has its own installer for each operating system (Linux, Windows, MacOS X, etc)
    * Note for Linux / MacOS X users: Ensure the matlab binary is in the system path. For example, you can add a symlink from a directory already in the path. Assuming your `matlab` binary is in `/opt/matlab/bin/`, then
```
$ cd /usr/local/sbin
$ sudo link -s /opt/matlab/bin/matlab
```
    * Note for Windows users: Similarly to the previous step, add to the `Path` environment variable the path to your Matlab `LIBMEX.DLL` and `LIBMX.DLL` libraries, e.g. "C:\Program Files\MATLAB\R2010b\bin\win32" (32 bit) or "C:\Program Files\MATLAB\R2010b\bin\win64" (64 bit)
  * [Gerardus](http://code.google.com/p/gerardus/)
    1. Check out the Gerardus project code
    * Note for Linux and Mac users: If you are not going to commit any code back, you can check out the latest version of the code from the command line doing
```
$ svn co http://gerardus.googlecode.com/svn/trunk/ gerardus
```
> > > If you are going to commit code, then you will need to use `https` instead of `http`
```
$ svn co https://gerardus.googlecode.com/svn/trunk/ gerardus
```
    * Note for Windows and Mac users: If you are using a graphic interface like TortoiseSVN, choose "Checkout" from the menu options, and enter `https://gerardus.googlecode.com/svn/trunk/` (commiters) or `http://gerardus.googlecode.com/svn/trunk/` (read only) as the repository's URL. Name the destination folder `gerardus`

# Build and install Gerardus #

## Instructions for Linux and Mac users ##

  1. Run from the command line
```
# Create a build directory for Gerardus
$ cd gerardus
$ mkdir bin
$ cd bin
# First configuration run (this will download and build ITK)
$ cmake ..
# Second configuration run (this will configure all Gerardus)
$ cmake ..
# Build and install Gerardus code (add e.g. flag -j6 to use 6 processors)
$ make install
```
  1. If the `matlab` command is not in the path, Gerardus will not know where to find Matlab. In that case, it's necessary to provide to the cmake lines above a flag with the path to the Matlab root directory, e.g. `-DMATLAB_ROOT=/usr/local/MATLAB/R2012a/`
  1. Open Matlab from the `matlab` directory and run script `add_gerardus_paths.m` to add the Gerardus toolboxes to the Matlab path
```
# from the linux command line
$ cd ../matlab
$ ./matlab &
# from the Matlab command line
>> add_gerardus_paths.m
```
  1. If you get an error message similar to this when running a MEX file
```
Invalid MEX-file
'/home/jdoe/workspace/gerardus/matlab/CgalToolbox/cgal_meshseg.mexa64':
/usr/local/MATLAB/R2012b/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6:
version `GLIBCXX_3.4.15' not found (required by
/home/jdoe/workspace/gerardus/matlab/CgalToolbox/cgal_meshseg.mexa64)
```

> the problem is that your system `libstdc++.so` (linked to by the ITK libraries) is [slightly different from Matlab `libstdc++.so`](https://groups.google.com/forum/#!topic/gerardus-users/48vhO1TWhkg) (linked to by the `libmex.so` library). Even if both `libstdc++.so` correspond to e.g. g++ v4.4, they have different minor versions, and that's were the error comes from. The solution is to delete the soft link in Matlab that points to Matlab's version of the library, and replace it by a soft link to the system's library. Library locations and names will vary between Linux distributions and Matlab versions, but for example, for Matlab R2012b with a default install:
```
cd /usr/local/MATLAB/R2012b/sys/os/glnxa64/
sudo rm libstdc++.so.6
sudo ln -s /usr/lib/gcc/x86_64-linux-gnu/4.4/libstdc++.so
```

## Instructions for Windows users ##

  1. Open file `gerardus\CMakeLists.txt` with your CMake program
  1. Click on the "Configure" button. Specify the generator for this project (the one we have tested is "Visual Studio 10" with "Use default native compilers").
  1. Click on the "Generate" button
  1. Open file `gerardus\bin\GERARDUS.sln` with your development environment (e.g. Visual Studio 10)
  1. Select "Release" from the configuration pull-down menu above the Solution Explorer, instead of "Debug"
  1. Select "Build -> Build Solution" from the menu. This will build the libraries and MEX files
  1. In the "Solution Explorer", right click on "INSTALL" and select "Build". This will install the libraries and MEX files
  1. Click on "Start" -> "My computer" -> "View system information", and in the "Advanced" tab click on "Environment Variables". Select the System variable `Path`, and click "Edit". Add the path to the `ITKCOMMON-*.DLL` library, e.g. "C:\Program Files\ITK\bin". Make sure that this variable also contains the path to where the Matlab libmex.dll file is, e.g. "`C:\Program Files\MATLAB\R2012\bin\win64`", or add them if necessary. After editing the Path variable you'll need to reboot Matlab.
  1. Launch Matlab, change directory to `gerardus\test`, and run script `add_gerardus_paths.m` to add the Gerardus toolboxes to the Matlab path, and the directories with DLLs to the system path
```
% from the Matlab command line
>> add_gerardus_paths.m
```


# Notes on MEX functions provided by Gerardus #

Gerardus contains Matlab function that can be run directly without building any code. There are also C++ programs and C++ MEX files that need to be compiled and installed before they can be run. These instructions explain how to that.

After building and installing Gerardus, you should have the command line programs and Matlab MEX functions listed below.

Command line programs (file names given following the linux convention. In Windows, executables have the extension `.exe`):

  * `gerardus/programs/`
    * `extractVoxelCoordinatesFromSegmentationMask`
    * `padSegmentationMaskWithVoxels`
    * `resize3DImage`
    * `rigidRegistration2D`
    * `rotate3DImage`
    * `skeletonize3DSegmentation`
    * `vesselness3DImage`

Matlab MEX functions (file names given following the linux 64 bit convention. In linux 32 bit, Windows 32/64 bit and Mac 32/64 bit the extension is different):

  * `gerardus/matlab/`
    * `PointsToolbox/sparse_breakdown.mexa64`
    * `PointsToolbox/mba_surface_interpolation.mexa64`
    * `ItkToolbox/itk_imfilter.mexa64`
    * `ItkToolbox/itk_registration.mexa64`
    * `ItkToolbox/itk_pstransform.mexa64`
    * `ThirdPartyToolbox/dijkstra.mexa64`
    * `FiltersToolbox/bwregiongrow.mexa64`
    * `FiltersToolbox/im2dmatrix.mexa64`
    * `CgalToolbox/cgal_insurftri.mexa64`
    * `CgalToolbox/cgal_closest_trifacet.mexa64`