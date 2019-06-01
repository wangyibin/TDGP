#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
The format of the text, such as font color.
"""

class bcolor(object):
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[91m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

FONTCOLORS = dict(zip(
    ('grey', 'red', 'green', 'yellow', 'blue', 'magenta','cyan', 'white',),
    ["0;%d" % x for x in range(30, 38)]))
END = '\033[0m'

class ColoredText:
    """
    Set font color.
    """
    def __init__(self, text, attr=None):
        self.text = text 
        self.attr = attr

    def __str__(self):
        """Colorize text.

        Available text color:
            red, green, yellow, blue, magenta, cyan, white.
        
        Example:
            ColoredText('Hello World!', 'red')
        """
        ctext = None

        fmt_str = '\033[%sm%s'

        if self.attr:
            ctext = fmt_str % (FONTCOLORS[self.attr],self.text)
        ctext += END
    
        return ctext or self.text
    
    __repr__ = __str__


grey = lambda s: str(ColoredText(s,'grey'))
red = lambda s: str(ColoredText(s, 'red'))
green = lambda s: str(ColoredText(s, 'green'))
yellow = lambda s: str(ColoredText(s, 'yellow'))
blue = lambda s: str(ColoredText(s, 'blue'))
magenta = lambda s: str(ColoredText(s,'magenta'))
cyan = lambda s: str(ColoredText(s, 'cyan'))
white = lambda s: str(ColoredText(s, 'white'))



def test():
    print(red('hello'))
    print(yellow('hello'))


if __name__ == "__main__":
    test()
