#!/usr/bin/env python
# Copyright (C) 2012 Dave Abrahams <dave@boostpro.com>
#
# Distributed under the Boost Software License, Version 1.0.  See
# accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt

import sys_path_setup
import os
import sys
import shutil
from tempdir import TempDir
from subprocess import *
from path import Path

home = Path(os.environ['HOME'])
feeds = home / 'src' / 'ryppl' / 'feeds'
if '--src' in sys.argv[1:]:
    dump_dir = feeds / 'src_dumps'
else:
    dump_dir = feeds / 'dumps'

boost_zero = home / 'src' / 'ryppl' / 'boost-zero'

if __name__ == '__main__':
    with TempDir() as build_dir:

        if os.path.exists(dump_dir):
            print '### deleting old dumps...'
            shutil.rmtree(dump_dir)
        os.makedirs(dump_dir)

        try:
            os.chdir(build_dir)
            print '### build_dir =',build_dir
            print '### CMake first pass...'
            args = ['cmake', '-DRYPPL_PROJECT_DUMP_DIRECTORY='+dump_dir]
            if '--src' not in sys.argv[1:]:
                args += ['-DRYPPL_DISABLE_TESTS=1', '-DRYPPL_DISABLE_DOCS=1' ]
            args += [boost_zero]

            env = os.environ.copy()
            env['PATH'] += os.pathsep + '/opt/local/lib/openmpi/bin'
            p = Popen(args, env=env, stdout=open(os.devnull,'w'), stderr=PIPE)
            stderr = p.communicate()[1]

            assert p.returncode != 0, 'Expected to exit with an error after initial configuration pass'
            if not 'Initial pass successfully completed, now run again!' in stderr:
                print stderr
                raise CalledProcessError(p.returncode, 'cmake', output=stderr)

            print '### CMake second pass...'
            check_call(args, env=env)
        finally:
            print '### Cleaning source directories...'
            os.chdir(boost_zero)
            check_call(['git', 'submodule', 'foreach', 
                        'if [[ "$(pwd)" != %r ]] ; then git clean -f -d ; fi' % str(boost_zero/'ryppl')]
                       , stdout=open(os.devnull,'w')
                       )
        
        print '### Done.'
