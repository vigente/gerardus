

# Summary #

The [Filters Toolbox](https://code.google.com/p/gerardus/source/browse/trunk/matlab/FiltersToolbox) contains functions that are applied to images (e.g. edge detection, segmentation, image enhancement).

# List of functions #

## skeleton\_label ##

Code: [skeleton\_label.m](https://code.google.com/p/gerardus/source/browse/trunk/matlab/FiltersToolbox/skeleton_label.m)

```
% skeleton_label  Give each branch of a skeleton a different label, and
% sort the voxels within each branch.
```

Script: [test\_skeleton\_label.m](https://code.google.com/p/gerardus/source/browse/trunk/matlab/test/test_skeleton_label.m).

Load binary image:

```
im = imread('circles.png');
imagesc(im)
axis equal off
colormap(gray)
```

<img width='50%' src='https://gerardus.googlecode.com/svn/wiki/figures/test_skeleton_label-circles.png'>

Compute skeleton:<br>
<br>
<pre><code>im = bwmorph(im, 'skel', Inf);<br>
imagesc(im)<br>
axis equal off<br>
colormap(gray)<br>
</code></pre>

<img width='50%' src='https://gerardus.googlecode.com/svn/wiki/figures/test_skeleton_label-skel.png'>

Split the skeleton into branches, and assign a different label to each branch:<br>
<br>
<pre><code>im = skeleton_label(im);<br>
imagesc(im)<br>
axis equal off<br>
colormap('default')<br>
</code></pre>

<img width='50%' src='https://gerardus.googlecode.com/svn/wiki/figures/test_skeleton_label-lab.png'>