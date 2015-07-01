Documentation for the Elastix Toolbox.
========================================

Install
----------------

The Elastix Toolbox is part of project Gerardus at the Institute of
Biomedical Engineering (IBME), University of Oxford. The most recent
development copy can be found at

https://github.com/rcasero/gerardus/tree/master/matlab/ElastixToolbox

The toolbox needs to be in Matlab's path.

In addition, you need to install in your system the binaries from
project elastix developed by Stefan Klein & Marius Staring

http://elastix.isi.uu.nl/index.php

To make the binaries available to Matlab, e.g. in linux you can
install them to a directory called /opt/elastix, and then link to the
binaries from /usr/local/bin/.

```bash
ln -s  /opt/elastix/elastix /usr/local/bin/elastix
ln -s  /opt/elastix/transformix /usr/local/bin/transformix
```

Quick guide
-----------

There are two main binaries in project elastix:

* 'elastix:' compute a transform to register two images
* 'transformix:' apply a transform to an image

This toolbox provides, correspondingly, functions elastix.m and
transformix.m.

The elastix project binaries rely on text files to describe input
registration parameters or an output transform. Our Matlab functions
use structs.

We provide two functions to convert between files and structs:

```matlab
>> t = elastix_read_file2param('t.txt');
>> file = elastix_write_param2file(t);
```

The images to be registered can have any number of channels, and can
be provided as arrays or file names. For example, as arrays

```matlab
>> imf = imread('fixed.png');
>> imm = imread('moving.png');
>> regParam = elastix_read_file2param('rigid.txt');
>> [t, immreg] = elastix(regParam, imf, imm);
```

The second example is as filenames

```matlab
>> [t, immreg] = elastix('rigid.txt', 'fixed.png', 'moving.png');
```

Transforms can be applied to an image also in array format

```matlab
>> t = elastix_read_file2param('t.txt');
>> imout = transformix(t, im);
```

or in file name format

```matlab
>> imfile = transformix('t.txt', 'im.png');
```

Devil in the details of composition of transformations
------------------------------------------------------

Let IMF be the fixed image, and IMM0 be a moving image.

Let T0 be an initial transform for IMM0. There are two ways we can
apply the initial transform, and it produces different results.

'Case 1:' We pass the initial transform with the -t0 parameter, and
use IMM0 for the registration.

```matlab
>> OPT.t0 = T0;
>> T1 = elastix(PARAM, IMF, IMM0, OPT);
```

In this case, T1.InitialTransformParametersFileName = T0. That is, T0
is the initial transform followed by T1. The registered image is

```matlab
>> IMMREG = transformix(T1, IMM0);
```

'Case 2:' We apply the transform explicitly to IMM0 and then
apply the registration to the result

```matlab
>> IMM = transformix(T0, IMM0);
>> T2 = elastix(PARAM, IMF, IMM);
```

In this case, T2.InitialTransformParametersFileName =
'NoInitialTransform'. If we want to concatenate both transforms, we
need to do T0.InitialTransformParametersFileName = T2. Note how the
order has swapped. Now, T2 is the initial transform, followed by T0.

The same registered image can be obtained with these two syntaxes:

```matlab
>> IMMREG = transformix(T2, IMM);
```

or

```matlab
>> T0.InitialTransformParametersFileName = T2;
>> IMMREG = transformix(T0, IMM0);
```

Note also that T1.TransformParameters is different from
T2.TransformParameters.
