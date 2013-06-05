import os
from path import *
from subprocess import check_call
from tempdir import *

__all__ = ['sign_feed']

# From http://stackoverflow.com/a/379535/125349
# by Suraj Barkale
def which(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate

    return None

def sign_feed(feed_path):
    # hook gpg to suppress annoying output
    os_path = [Path(x) for x in os.environ['PATH'].split(Path.pathsep)]
    with TempDir() as hook_dir:
        wrapper = hook_dir/'gpg'
        open(wrapper, 'w').write(
'''#!/usr/bin/env bash
%r --batch "$@"
''' % str(which('gpg')))
        os.chmod(wrapper, 0744)

        env = os.environ.copy()
        env['PATH'] = os.pathsep.join([hook_dir] + os_path)
        check_call(['0publish', '--xmlsign', feed_path], env=env)
