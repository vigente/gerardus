# Copyright Dave Abrahams 2012. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
"""Supplies the main entry point for the ryppl command"""

import argparse
import commands
import sys
import logging
from commands import *

def run():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-commands')
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument('--verbose', '-v', action='append_const', const=1)
    verbosity.add_argument('--quiet', '-q', action='append_const', const=1)

    for module_name in commands.__all__:
        cmd_module = getattr(commands, module_name)

        cmd_name = module_name.replace('_', '-')[int(module_name.startswith('_')):]

        kw = dict(description=getattr(cmd_module.command_line_interface, '__doc__', None))

        subparser = subparsers.add_parser(cmd_name, **kw)
        subparser.set_defaults(runner=cmd_module.run)
        cmd_module.command_line_interface(subparser)

    args = parser.parse_args(sys.argv[1:])

    log_level = logging.WARN
    if args.quiet:
        log_level += 10*len(args.quiet)
    if args.verbose:
        log_level -= 10*len(args.verbose)
    logging.getLogger().setLevel(log_level)

    args.runner(args)
    
