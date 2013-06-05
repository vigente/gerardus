# Copyright Dave Abrahams 2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
from logging import info
import subprocess
import os

command = ('cmake',)

def cmake(args, **kw):
    info('cmake %s', args)
    return subprocess.check_call(command + tuple(args), **kw)

def configure_for_circular_dependencies(*args, **kw):
    info('CMake first pass')
    cmd = command + tuple(args)
    p = subprocess.Popen(
        cmd,
        stdout=open(os.devnull,'w'), 
        stderr=subprocess.PIPE, 
        **kw)
    stderr = p.communicate()[1]
    if p.wait() == 0 or 'Initial pass successfully completed, now run again!' not in stderr:
        raise subprocess.CalledProcessError(
            p.returncode, cmd,
            output=stderr if p.returncode  else
            'Expected to exit with an error after initial configuration pass\n' + stderr)

    info('CMake second pass')
    cmake(args, **kw)

def generators():
    help = subprocess.check_output(command + ('--help',))

    generator_block = help.split('The following generators are available on this platform:')[1]
    ret = []
    for line in generator_block.replace('\r\n','\n').split('\n'):
        g = line.split('=')[0].strip()
        if g:
            ret.append(g)
    return ret

if __name__ == '__main__':
    print generators()
    configure_for_circular_dependencies('-G', 'Xcode', '../src')

