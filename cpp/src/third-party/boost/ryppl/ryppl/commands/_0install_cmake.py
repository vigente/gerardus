# Copyright Dave Abrahams 2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
import sys
from ryppl.support import _argparse
from ryppl.support.cmake import cmake
from os import environ as env
import os
from logging import info
from multiprocessing import cpu_count

def command_line_interface(cli):
    '''Used in 0install feeds for building CMake-based Ryppl projects from source'''

    cli.add_argument(
        '--overlay'
      , type=_argparse.existing_directory
      , help='A directory to be overlaid on the composite source directory before building')

    cli.add_argument(
        '--add-subdirectory', nargs=2, action='append', metavar=('SOURCEDIR','SUBDIR')
      , help='''Copy SOURCEDIR into SUBDIR of the composite source
 directory and generate a top-level CMakeLists.txt file that includes
 it.  This is sometimes necessary to handle evil dependency
 clusters.'''
        )

    cli.add_argument(
        '--build-type', default='Release'
      , help='Passed to CMake with -DCMAKE_BUILD_TYPE=%(metavar)s'
        )

    cli.add_argument(
        'components'
      , nargs='+'
      , help='''The CMake components to be built.  If subdirectories have
been added with --add-subdirectory, each component should be
prefixed by its subdirectory and a slash, e.g. "regex/dev".''')

def run(args):
    if args.overlay or args.add_subdirectory:
        SRCDIR = os.path.join(os.getcwd(),'ryppl-composite-src')
        info('preparing source directory...')
    else:
        SRCDIR = env['SRCDIR']

    DISTDIR = env['DISTDIR']

    if args.overlay:
        cmake(['-E', 'copy_directory', env['SRCDIR'], SRCDIR])
        cmake(['-E', 'copy_directory', args.overlay, SRCDIR])

    if args.add_subdirectory:
        with open(os.path.join(SRCDIR,'CMakeLists.txt'), 'w') as root_cmakelists_txt:
            root_cmakelists_txt.write('cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)\n')
            for src,dst in args.add_subdirectory:
                cmake(['-E', 'copy_directory', src, os.path.join(SRCDIR, dst)])
                root_cmakelists_txt.write('add_subdirectory(%s)\n' % dst)

    build_docs = 'doc' in (x.split('/')[-1] for x in args.components)

    info('configuring...')
    for i in range(2):
        cmake([
            '-DRYPPL_INITIAL_PASS=%d' % (i == 0)
          , '-DCMAKE_MODULE_PATH=' + env['RYPPL_CMAKE_MODULE_PATH']
          , '-DCMAKE_BUILD_TYPE='+ args.build_type
          , '-DRYPPL_DISABLE_TESTS=1'
          , '-DRYPPL_DISABLE_EXAMPLES=1'
          , '-DRYPPL_DISABLE_DOCS=%d' % int(not build_docs)
          , '-DBUILD_SHARED_LIBS=1'
          , SRCDIR
          ])

    info('building...')
    build_args = ['--build', '.']
    if build_docs:
        build_args += ['--target', 'documentation']
    if os.path.exists('Makefile'):
        build_args+= ['--', '-j%d' % (cpu_count() + 1)]
    cmake(build_args)
    
    for c in args.components:
        subdir_component = c.split('/')
        info('installing %s...', c)
        if len(args.components) > 1:
            prefix = os.path.join(DISTDIR, *subdir_component)
        else:
            prefix = DISTDIR
        component = subdir_component[-1]
        
        cmake(
            ([] if component == 'preinstall' else ['-DCOMPONENT='+component])
            + [ '-DCMAKE_INSTALL_PREFIX='+prefix, '-P', 'cmake_install.cmake' ]
            )


