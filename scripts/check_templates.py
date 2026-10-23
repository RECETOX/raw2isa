"""Check that the assay template header and data rows have matching lengths."""
import sys

for path in [
    "raw2isa/templates/a_raw.txt",
    "raw2isa/templates/s_raw.txt",
]:
    with open(path) as f:
        lines = f.read().splitlines()
    lengths = [len(line.split("\t")) for line in lines]
    print(path, lengths)
    if len(set(lengths)) != 1:
        print("MISMATCH")
        sys.exit(1)
print("OK")
