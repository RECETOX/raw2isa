"""Miscellaneous utilities for raw2isa.

Most of the code in this module was adapted from the ``mzml2isa`` package,
which is distributed under the GNU General Public License version 3.0.

License:
    GNU General Public License version 3.0 (GPLv3)
"""

import functools
import string


class PermissiveFormatter(string.Formatter):
    """A formatter that replaces wrong and missing keys with a blank."""

    def __init__(self, missing="", bad_fmt=""):
        self.missing = missing
        self.bad_fmt = bad_fmt

    def get_field(self, field_name, args, kwargs):
        # Handle a key not found
        try:
            val = super().get_field(field_name, args, kwargs)
        except (KeyError, AttributeError, IndexError, TypeError):
            val = None, field_name
        return val

    def format_field(self, value, spec):
        # handle an invalid format
        if value is None:
            return self.missing
        try:
            return super().format_field(value, spec)
        except ValueError:
            if self.bad_fmt is not None:
                return self.bad_fmt
            else:
                raise


def star_args(func):
    """Unpack arguments before calling the wrapped function."""

    @functools.wraps(func)
    def new_func(args):
        return func(*args)

    return new_func
