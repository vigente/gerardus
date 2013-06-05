# Copyright Dave Abrahams 2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
from zeroinstall.injector import cli
from subprocess import check_call
import sys

command = (
    sys.executable
    , '-c', 'import sys\n'
    'from zeroinstall.injector import cli\n'
    'cli.main(sys.argv[1:])\n'
)

def launch(args, **kw):
    '''Effectively calls 0launch on with the given command-line
    arguments.  Keyword arguments should either be noreturn=True,
    indicating that this process should be *replaced* by whatever
    0launch invokes, or will be forwarded to subprocess.Popen and
    interpreted accordingly.'''

    if kw.pop('noreturn', None):
        assert len(kw) == 0
        # The caller has told us he doesn't need to regain control, so
        # launch the command directly.
        cli.main(args)
    else:
        # Run 0launch in a subprocess so we get control back
        # afterwards.  Otherwise the 0launch process *replaces* itself
        # with the process being launched (on posix).
        check_call(command + tuple(args), **kw)

