# Copyright Dave Abrahams 2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
import os
from path import Path

def executable_path(name):
    for d in os.environ['PATH'].split(Path.pathsep):
        p = Path(d)/name
        if os.name == 'nt':
            for ext in os.environ['PATHEXT'].split(Path.pathsep):
                if os.path.isfile(p+ext):
                    return p+ext
        elif os.path.isfile(p):
            return p
