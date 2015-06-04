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
