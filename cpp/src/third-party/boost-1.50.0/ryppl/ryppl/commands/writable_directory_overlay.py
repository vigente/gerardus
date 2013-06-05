# Copyright Dave Abrahams 2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

from ryppl.support._argparse import existing_directory
import os
import shutil
import stat

def command_line_interface(cli):
    '''Overlay the contents of the source directory on the contents of the
    destination directory, leaving all new files writable'''

    cli.add_argument(
        'srcdir'
      , type=existing_directory
      , help='A directory to be overlaid on DSTDIR')

    cli.add_argument(
        'dstdir'
      , type=existing_directory
      , help='A directory on which to overlay SRCDIR')

def overlay(srcdir, dstdir):

    def copy_files(_, src_subdir, fnames):
        dst_subdir = os.path.join(dstdir, os.path.relpath(src_subdir, srcdir))

        if not os.path.isdir(dst_subdir):
            os.makedirs(dst_subdir)

        for fname in fnames:
            src_file = os.path.join(src_subdir,fname)
            if os.path.isdir(src_file):
                continue

            dst_file = os.path.join(dst_subdir,fname)

            shutil.copy2(src_file,dst_file)
            st = os.stat(dst_file)
            if not (st.st_mode & stat.S_IWRITE):
                os.chmod(dst_file, st.st_mode | stat.S_IWRITE)

    os.path.walk(srcdir, copy_files, None)

def run(args):
    overlay(args.srcdir, args.dstdir)

def test():
    # A (very) simple test of efficacy
    import tempfile
    root = tempfile.mkdtemp()
    srcdir = os.path.join(root, 'src')
    dstdir = os.path.join(root, 'dst')
    os.makedirs(os.path.join(srcdir,'y'))
    open(os.path.join(srcdir,'y','z'), 'w')
    os.makedirs(dstdir)
    src_file = os.path.join(srcdir, 'x')
    open(src_file, 'w').write('foo')
    os.chmod(src_file, stat.S_IREAD)
    assert not os.stat(src_file).st_mode & stat.S_IWRITE
    overlay(srcdir, dstdir)
    assert os.path.isfile(os.path.join(dstdir,'y','z'))
    assert os.stat(os.path.join(dstdir,'x')).st_mode & stat.S_IWRITE

if __name__ == '__main__':
    test()
