"""Higher-level interface to parser/converter functionalities.

This module exposes the basic API of raw2isa, either being called from the
command line interface with argument parsing via the **main** function, or
from another Python program via the **convert** function which works the
same.

Example:
    >>> from raw2isa.parsing import convert
    >>> convert('examples/raw', 'output_directory', 'STUDY1')

License:
    GNU General Public License version 3.0 (GPLv3)
"""

import argparse
import contextlib
import multiprocessing.pool
import os
import sys
import warnings

from . import __version__
from ._impl import tqdm
from .isa import ISA_Tab
from .raw import ThermoRawFile
from .usermeta import UserMetaLoader


def _raw_files_in(directory):
    """List the RAW files of a local directory, case-insensitively.

    Arguments:
        directory (str): the directory to scan.

    Returns:
        list: the sorted list of RAW file paths.
    """
    return sorted(
        (
            os.path.join(directory, entry.name)
            for entry in os.scandir(directory)
            if entry.is_file() and entry.name.lower().endswith(".raw")
        ),
        key=str.lower,
    )


def _parse_file(path):
    """Parse a single RAW file and return its metadata dictionary.

    Arguments:
        path (str): the local filesystem path to the RAW file.

    Returns:
        dict: a dictionary containing the extracted metadata.
    """
    return ThermoRawFile(path).metadata


def convert(
    in_path,
    out_path,
    study_identifier,
    usermeta=None,
    split=True,
    jobs=1,
    template_directory=None,
    verbose=True,
):
    """Parse a study from the given *in_path* and then create an ISA-Tab study.

    Arguments:
        in_path (str): path to the directory containing the Thermo RAW files
        out_path (str): path to the output directory
        study_identifier (str): study identifier (e.g. MTBLSxxx)

    Keyword Arguments:
        usermeta (str, optional): the path to a json file, a xlsx file or
            directly a json formatted string containing user-defined
            metadata [default: None]
        split (bool): split assay files based on the polarity of the scans
            [default: True]
        jobs (int): the number of jobs to use for parsing RAW files.
            Note that the .NET runtime used by ``epyraw`` does not always
            behave well with multithreading; keep 1 if in doubt.
            [default: 1]
        template_directory (str, optional): the path to a directory
            containing custom templates to use when writing ISA tab
            [default: None]
        verbose (bool): display more output [default: True]
    """
    # open user metadata file if any
    meta_loader = UserMetaLoader(usermeta)

    # check the input directory
    if not os.path.isdir(in_path):
        raise NotADirectoryError(
            f"{in_path} is not a directory: raw2isa needs a local folder containing "
            "the Thermo RAW files to parse."
        )

    # get all RAW files
    paths_iter = _raw_files_in(in_path)

    if not paths_iter:
        warnings.warn(f"No files were found in {in_path}.", UserWarning, stacklevel=2)
        return

    # wrap in a progress bar if needed
    if not verbose and tqdm is not None:
        paths_iter = tqdm.tqdm(paths_iter)

    # parse using threads if needed
    if jobs > 1:
        with contextlib.closing(multiprocessing.pool.ThreadPool(jobs)) as pool:
            metalist = pool.map(_parse_file, paths_iter)
    else:
        metalist = list(map(_parse_file, paths_iter))

    # write ISA-Tab files
    if metalist:
        if verbose:
            print("Parsing RAW meta information into ISA-Tab structure")
        isa_tab = ISA_Tab(
            out_path,
            study_identifier,
            usermeta=meta_loader.usermeta,
            template_directory=template_directory,
        )
        isa_tab.write(metalist, "raw", split=split)


def main(argv=None):
    """Run **raw2isa** from the command line.

    Arguments:
        argv (list, optional): the list of arguments to run raw2isa with
            (if None, then sys.argv is used) [default: None]
    """
    p = argparse.ArgumentParser(
        prog="raw2isa",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Extract meta information from Thermo RAW files and create ISA-tab structure""",
        usage="raw2isa -i IN_PATH -o OUT_PATH -s STUDY_ID [options]",
    )

    p.add_argument(
        "-i",
        dest="in_path",
        help="input folder containing Thermo RAW files",
        required=True,
    )

    p.add_argument(
        "-o",
        dest="out_path",
        help="out folder (new files will be created here)",
        required=True,
    )
    p.add_argument(
        "-s",
        dest="study_id",
        help="study identifier (e.g. MTBLSxxx)",
        required=True,
    )
    p.add_argument(
        "-m",
        dest="usermeta",
        help="additional user provided metadata (JSON format)",
        default=None,
        required=False,
    )
    p.add_argument(
        "-j",
        dest="jobs",
        help="launch different threads for parsing",
        action="store",
        required=False,
        default=1,
        type=int,
    )
    p.add_argument(
        "-n",
        dest="split",
        help="do NOT split assay files based on polarity",
        action="store_false",
        default=True,
    )
    p.add_argument(
        "-W",
        dest="wrng_ctrl",
        help="warning control (with python default behaviour)",
        action="store",
        default="once",
        required=False,
        choices=["ignore", "always", "error", "default", "module", "once"],
    )
    p.add_argument(
        "-t",
        dest="template_dir",
        help="directory containing default template files",
        action="store",
        default=None,
    )
    p.add_argument(
        "--version", action="version", version=f"raw2isa {__version__}"
    )
    p.add_argument(
        "-v",
        dest="verbose",
        help="show more output (default if tqdm is not installed)",
        action="store_true",
        default=False,
    )

    args = p.parse_args(argv or sys.argv[1:])

    if not tqdm:
        args.verbose = True

    if args.verbose:
        print(f"{os.linesep} input path: {args.in_path}")
        print(f"output path: {args.out_path}")
        print(f"Study identifier:{args.study_id}{os.linesep}")

    with warnings.catch_warnings():
        warnings.filterwarnings(args.wrng_ctrl)
        convert(
            args.in_path,
            args.out_path,
            args.study_id,
            usermeta=args.usermeta,
            split=args.split,
            verbose=args.verbose,
            jobs=args.jobs,
            template_directory=args.template_dir,
        )


if __name__ == "__main__":
    main()
