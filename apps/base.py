#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Basic support for running library as script
"""


from __future__ import print_function

import logging
import os
import os.path as op
import sys

from glob import glob

from TDGP.apps.font import *

from TDGP import __copyright__, __version__

TDGPHELP = "TDGP utility libraries v{} [{}]\n".format(__version__, __copyright__)


# logging output color format
# https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
COLORS = {
    'WARNING': yellow,
    'INFO': green,
    'DEBUG': green,
    'CRITICAL': yellow,
    'ERROR': red
}
class ColoredFormatter(logging.Formatter):

    def __init__(self, msg, use_color=True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color
    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            color_level = COLORS[levelname](levelname)
            record.levelname = color_level
        return logging.Formatter.format(self, record)

class ColoredLogger(logging.Logger):
    formats = magenta("%(asctime)s <%(module)s>") +"[%(levelname)s]"+ green(" %(message)s")
    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.formats)
        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)


def debug(level=logging.DEBUG):
    import logging
    logging.setLoggerClass(ColoredLogger)
    formats = magenta("%(asctime)s <%(module)s>")
    formats += yellow(" [%(levelname)s]")
    formats += green(" %(message)s")
    logging.basicConfig(level=level, format=formats, datefmt="%H:%M:%S")

'''
COLORS = {
    'WARNING': yellow,
    'INFO': green,
    'DEBUG': green,
    'CRITICAL': yellow,
    'ERROR': red
}
def debug(level=logging.DEBUG):
    """
    Basic config logging format
    """
    from TDGP.apps.font import magenta, green, yellow
    formats = magenta("%(asctime)s <%(module)s>")
    formats += yellow(" [%(levelname)s]")
    formats += green(" %(message)s")
    logging.basicConfig(level=level, format=formats, datefmt="%H:%M:%S")
'''

debug()

def main():
    pass


def dmain(mainfile, type="action"):
    """
    modify from https://github.com/tanghaibao/jcvi/blob/master/jcvi/apps/base.py
    """
    cwd = op.dirname(mainfile)
    pyscripts = [x for x in glob(op.join(cwd, "*", '__main__py'))] \
        if type == "module" \
        else glob(op.join(cwd, "*.py"))
    actions = []
    for ps in sorted(pyscripts):
        action = op.basename(op.dirname(ps)) \
            if type == 'module' \
            else op.basename(ps).replace(".py", "")
        if action[0] == "_":
            continue
        pd = get_module_docstring(ps)
        action_help = [x.rstrip(":.,\n") for x in pd.splitlines(True) \
                       if len(x.strip()) > 10 and x[0] != "%"][0] \
            if pd else "no docstring found"
        actions.append((action, action_help))

    a = ActionDispatcher(actions)
    a.print_help()


def splitall(path):
    """
    split all path and return a list
    >>> splitall("/code/TDGP/utils")
    ["code", "TDGP", "utils"]
    """

    allparts = []
    while True:
        path, p1 = op.split(path)
        if not p1:
            break
        allparts.append(p1)
    allparts = allparts[::-1]
    return allparts



class ActionDispatcher(object):
    """
    The action dispatch function.
    Copy from jcvi(https//:github.com/tanghaibao/jcvi)
    """

    def __init__(self, actions):

        self.actions = actions
        if not actions:
            actions = [(None, None)]
        self.valid_actions, self.action_helps = zip(*actions)

    def get_meta(self):
        args = splitall(sys.argv[0])[-3:]
        args[-1] = args[-1].replace(".py", "")
        if args[-2] == "bioway":
            meta = "MODULE"
        elif args[-1] == "__main__":
            meta = "SCRIPT"
        else:
            meta = "ACTION"
        return meta, args

    def print_help(self):
        meta, args = self.get_meta()
        if meta == "MODULE":
            del args[0]
            args[-1] = meta
        elif meta == "SCRIPT":
            args[-1] = meta
        else:
            args[-1] += " " + meta

        help = "Usage:\n    python -m {0}\n\n\n".format('.'.join(args))
        help += "Available {0}s:\n".format(meta)
        max_action_len = max(len(action) for action, ah in self.actions)
        for action, action_help in sorted(self.actions):
            action = action.rjust(max_action_len + 4)
            help += " | ".join((action, action_help.capitalize() + "\n"))

        help += "\n" + TDGPHELP

        sys.stderr.write(help)
        sys.exit(1)

    def dispatch(self, globals):
        from difflib import get_close_matches
        meta = "ACTION"
        if len(sys.argv) == 1:
            self.print_help()

        action = sys.argv[1]

        if not action in self.valid_actions:
            print("[error] {0} not a valid {1}\n".format(action, meta),
                  file=sys.stderr)
            alt = get_close_matches(action, self.valid_actions)
            print("Did you mean one of these?\n\t{0}\n".
                  format(", ".join(alt)), file=sys.stderr)
            self.print_help

        globals[action](sys.argv[2:])

'''
class OptionParser(OptionP):
    """
    OptionParser modify from https://githup.com/tanghaibao/jcvi.git
    """

    def __init__(self, doc):
        OptionP.__init__(self, doc, epilog=BIOWAYHELP)
        pass

    def parse_args(self, args=None):
        dests = set()
        ol = []
        for g in [self] + self.option_groups:
            ol += g.option_list
        for o in ol:
            if o.dest in dests:
                continue
            self.add_help_from_choices(o)
            dests.add(o.dest)

        return OptionP.parse_args(self, args)

    def add_help_from_choices(self, o):
        if o.help == SUPPRESS_HELP:
            return

        default_tag = "%default"
        assert o.help, "Option {0} do not have help string".format(o)
        help_pf = o.help.capitalize()
        if "[" in help_pf:
            help_pf = help_pf.rsplit("[", 1)[0]
        help_pf = help_pf.strip()

        if o.type == "choice":
            if o.default is None:
                default_tag = "guess"
            ctext = "|".join(sorted(str(x) for x in o.choices))
            if len(ctext) > 100:
                ctext = ctext[:100] + "..."
            choice_text = "must be one of {0}".format(ctext)
            o.help = "{0}, {1} [default: {2}]".format(help_pf,
                                                      choice_text, default_tag)
        else:
            o.help = help_pf
            if o.default is None:
                default_tag = "disabled"
            if o.get_opt_string() not in ("--help", "--version") \
                    and o.action != "store_false":
                o.help += " [default:{0}".format(default_tag)
'''

def get_module_docstring(filepath):
    """
    get module docstring
    """
    co = compile(open(filepath).read(), filepath, 'exec')
    if co.co_consts and isinstance(co.co_consts[0], six.string_types):
        docstring = co.co_consts[0]
    else:
        docstring = None
    return docstring


