"""``Thermo RAW`` file parser and converter to ISA-Tab.

``raw2isa`` is a parser/converter that translates Thermo Fisher Scientific
``.raw`` mass spectrometry files[1]_ to ISA-TAB format (Investigation-Study-Assay
with TAB separated values)[2]_.

The metadata of the ``.raw`` files is extracted with the ``epyraw`` library,
a Python wrapper around the Thermo Fisher ``RawFileReader`` .NET library.

References:
    [1] https://planetorbitrap.com/rawfilereader
    [2] http://isa-tools.org/format.html

License:
    GNU General Public License version 3.0 (GPLv3)
"""
__author__ = "raw2isa contributors"
__version__ = "0.1.0"
__license__ = "GPLv3"

__all__ = ["isa", "parsing", "raw", "usermeta"]

try:
    from . import isa, parsing, raw, usermeta
except ImportError:
    pass
