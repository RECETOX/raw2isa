"""Minimal ISA-Tab serializer for raw2isa.

This module contains a single class, ``ISA_Tab``, which is used to dump a
list of ``ThermoRawFile.metadata`` dictionaries to ISA-Tab files. It was
adapted from the ``mzml2isa`` package (GPLv3).

License:
    GNU General Public License version 3.0 (GPLv3)
"""

import copy
import csv
import functools
import os
import re
from collections import ChainMap

from . import __version__, templates
from ._impl import resource_files
from .utils import PermissiveFormatter


class ISA_Tab:
    """Export a list of RAW metadata dictionaries to ISA-Tab files.

    Attributes:
        usermeta (dict): a dictionary containing metadata defined by the
            user, such as "Study Publication" or "Submission Date".
            Defaults to None.
        isa_env (dict): a dictionary containing various environment
            variables associated with the current ISA_Tab object (such as
            'Study file name', 'Written assays', etc.)
    """

    def __init__(self, out_dir, name, usermeta=None, **kwargs):
        """Setup the environments and the directories.

        Arguments:
            out_dir (str): the path to the output directory.
            name (str): the name of the *omics study to generate.

        Keyword Arguments:
            usermeta (dict, optional): a dictionary containing metadata
                defined by the user, such as "Study Publication" or
                "Submission Date". [default: None]
            template_directory (str, optional): the path to a directory
                containing custom ISA-Tab templates. Not all templates are
                required, so if for instance only your "a_raw.txt" is
                non-standard, then you only have to provide a new "a_raw.txt"
                in your custom template directory. If None, the templates
                shipping with raw2isa are used. [default: None]
        """
        usermeta = kwargs.get("usermeta", usermeta)
        template_default = resource_files(templates.__name__)
        template_directory = kwargs.get("template_directory") or template_default

        self.usermeta = copy.deepcopy(usermeta) if usermeta is not None else {}
        self._set_default_usermeta(name)
        self.isa_env = {
            "out_dir": out_dir,
            "Study Identifier": name,
            "Study file name": f"s_{name}.txt",
            "Assay polar file name": f"a_{name}_{{}}_metabolite_profiling_mass_spectrometry.txt",
            "Assay file name": f"a_{name}_metabolite_profiling_mass_spectrometry.txt",
            "default_path": template_default,
            "template_path": template_directory,
            "Technology type": [],
            "Measurement type": [],
            "Written assays": [],
            "Platform": [],
            "raw measurement": {
                "name": "metabolite profiling",
                "accession": "http://purl.obolibrary.org/obo/OBI_0000366",
                "ref": "OBI",
            },
            "raw technology": {
                "name": "mass spectrometry",
                "accession": "http://purl.obolibrary.org/obo/OBI_0000470",
                "ref": "OBI",
            },
        }

        self.isa_env["Converter"] = "raw2isa"
        self.isa_env["Converter version"] = __version__

    def _set_default_usermeta(self, name):
        """Populate minimal study metadata required by ISA validators."""
        defaults = {
            "investigation": {
                "title": name,
                "description": "Generated from Thermo RAW metadata by raw2isa.",
            },
            "study": {
                "title": name,
                "description": "Generated from Thermo RAW metadata by raw2isa.",
            },
        }
        for section, values in defaults.items():
            target = self.usermeta.setdefault(section, {})
            for key, value in values.items():
                target.setdefault(key, value)

    def write(self, metalist, datatype="raw", **kwargs):
        """Generate and write the ISA files.

        Arguments:
            metalist (list): a list of RAW metadata dictionaries
            datatype (str): the datatype of the study (only 'raw' is
                supported for now)

        Keyword Arguments:
            split (bool, optional): a boolean stating if assay files should
                be split based on their polarities. [default: True]
        """
        split = kwargs.get("split", True)

        platforms = [
            meta["Instrument"]["name"] for meta in metalist if "Instrument" in meta
        ]

        if len(set(platforms)) > 1:
            print(
                "WARNING: The RAW files are derived from multiple instrument types, this can be problematic"
                'as the "platform" used in the ISAcreator typically uses 1 instrument type per assay,'
                'please check the Investigation output file to ensure the correct "platform" has been assigned'
                "to the correct assay file."
            )

        platform_str = ", ".join(sorted(set(platforms)))

        if split:
            self.isa_env["Platform"] = [
                platform_str
            ] * self._count_polarities(metalist)
        else:
            self.isa_env["Platform"] = [platform_str]

        if not os.path.exists(self.isa_env["out_dir"]):
            os.makedirs(self.isa_env["out_dir"])

        h, d = self.make_assay_template(metalist, datatype)

        self.create_assay(metalist, h, d, split=split)
        self.create_study(metalist, datatype)
        self.create_investigation(metalist, datatype)

    @staticmethod
    def _count_polarities(metalist):
        polarities = set()
        for meta in metalist:
            polarity = meta.get("Scan polarity", {}).get("name")
            if polarity is None:
                polarity = "nopolarity"
            polarities.add(polarity)
        return max(len(polarities), 1)

    def _template_path(self, prefix, datatype):
        """Resolve a template file, allowing custom templates to be partial."""
        path = os.path.join(self.isa_env["template_path"], f"{prefix}_{datatype}.txt")
        if not os.path.exists(path):
            path = os.path.join(self.isa_env["default_path"], f"{prefix}_{datatype}.txt")
        return path

    def make_assay_template(self, metalist, datatype):
        """Generate the assay template rows.

        Parameters:
            metalist (list): a list of RAW metadata dictionaries
            datatype (str): the datatype of the study

        Returns:
            tuple: a tuple containing the list of assay headers and the
                list of assay data fields based on the cardinality of
                elements in the metalist
        """
        template_a_path = self._template_path("a", datatype)

        with open(template_a_path) as a_in:
            headers, data = [x.strip().replace('"', "").split("\t") for x in a_in.readlines()]

        i = 0
        while i < len(headers):
            header, datum = headers[i], data[i]

            if "{{" in datum and "Term" not in header:
                entry_list = metalist[0][self.unparameter(header)]["entry_list"]
                hsec, dsec = (
                    (headers[i : i + 3], data[i : i + 3])
                    if headers[i + 1] == "Term Source REF"
                    else (headers[i : i + 1], data[i : i + 1])
                )

                # remove the section we are going to format and re-insert
                headers[:] = headers[:i] + headers[i + len(hsec) :]
                data[:] = data[:i] + data[i + len(dsec) :]

                for n in reversed(range(len(entry_list))):
                    for (h, d) in zip(reversed(hsec), reversed(dsec), strict=False):
                        headers.insert(i, h)
                        data.insert(i, d.format(n))

            i += 1

        return headers, data

    def create_assay(self, metalist, headers, data, **kwargs):
        """Write the assay file(s), possibly split by polarity.

        Arguments:
            metalist (list): a list of RAW metadata dictionaries
            headers (list): the list containing the assay headers row
            data (list): the list containing the assay data row

        Keyword Arguments:
            split (bool, optional): a boolean stating if assay files should
                be split based on their polarities. [default: True]
        """
        split = kwargs.get("split", True)
        fmt = PermissiveFormatter()

        if split:
            polarities = (
                {m["Scan polarity"]["name"] for m in metalist if "Scan polarity" in m}
                if any("Scan polarity" in m for m in metalist)
                else {"nopolarity"}
            )
        else:
            polarities = {"nosplit"}

        new_a_path = (
            os.path.join(self.isa_env["out_dir"], self.isa_env["Assay file name"])
            if len(polarities) == 1
            else os.path.join(self.isa_env["out_dir"], self.isa_env["Assay polar file name"])
        )

        for polarity in sorted(polarities):
            if polarity == "nosplit":
                rows = metalist
            elif polarity == "nopolarity":
                rows = [x for x in metalist if "Scan polarity" not in x]
            else:
                rows = [
                    x
                    for x in metalist
                    if x.get("Scan polarity", {}).get("name") == polarity
                ]
            if not rows:
                continue

            csv_wopen = functools.partial(open, mode="w", newline="")

            with csv_wopen(new_a_path.format(polarity[:3].upper())) as a_out:
                self.isa_env["Written assays"].append(
                    os.path.basename(new_a_path.format(polarity[:3].upper()))
                )
                self.isa_env["Technology type"].append(self.isa_env["raw technology"])
                self.isa_env["Measurement type"].append(self.isa_env["raw measurement"])

                writer = csv.writer(
                    a_out, quotechar='"', quoting=csv.QUOTE_ALL, delimiter="\t"
                )
                writer.writerow(headers)

                for meta in rows:
                    writer.writerow(
                        [
                            fmt.vformat(x, None, ChainMap(meta, self.usermeta))
                            for x in data
                        ]
                    )

    def create_study(self, metalist, datatype):
        """Write the study file.

        Arguments:
            metalist (list): a list of RAW metadata dictionaries
            datatype (str): the datatype of the study
        """
        template_s_path = self._template_path("s", datatype)

        new_s_path = os.path.join(self.isa_env["out_dir"], self.isa_env["Study file name"])

        fmt = PermissiveFormatter()

        with open(template_s_path) as s_in:
            headers, data = s_in.readlines()

        with open(new_s_path, "w") as s_out:
            s_out.write(headers)
            for meta in metalist:
                s_out.write(fmt.vformat(data, None, ChainMap(meta, self.usermeta)))

    def create_investigation(self, metalist, datatype):
        """Write the investigation file.

        Arguments:
            metalist (list): a list of RAW metadata dictionaries
            datatype (str): the datatype of the study
        """
        template_i_path = self._template_path("i", datatype)

        new_i_path = os.path.join(self.isa_env["out_dir"], "i_Investigation.txt")

        meta = metalist[0]
        fmt = PermissiveFormatter()

        chained = ChainMap(self.isa_env, meta, self.usermeta)

        with open(template_i_path) as i_in:
            with open(new_i_path, "w") as i_out:
                for line in i_in:
                    if "{{" in line:
                        line, value = line.strip().split("\t")
                        label = value[3:].split("[")[0]

                        if label in chained:
                            for k in range(len(chained[label])):
                                line = "\t".join([line, value.format(k)])
                            line += "\n"
                        else:
                            line = f'{line}\t""\n'

                    line = fmt.vformat(line, None, chained)
                    i_out.write(line)

    @staticmethod
    def unparameter(string):
        """Extract string 's' from 'Parameter Value[s]'."""
        match = re.match(r"Parameter Value\[(.*)\]", string)
        return match.group(1) if match is not None else string
