import logging
import warnings
#logging.basicConfig(level=logging.INFO)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('numpy').setLevel(logging.WARNING)
logging.getLogger('numexpr').setLevel(logging.WARNING)

if not sys.warnoptions:
    warnings.simplefilter("ignore")

warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


__author__ = ("Yibin Wang")
__copyright__ = "Copyright (c) 2019, Yibin Wang"
__email__ = "yibinwang96@outlook.com"
__license__ = "BSD"
__status__ = "Development"
__version__ = "0.0.3"
