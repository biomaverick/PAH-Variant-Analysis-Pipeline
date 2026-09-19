#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

def parse_aa_change(change: str):
    """Parse 'R261Q' -> ('R', 261, 'Q')"""
    m = re.match(r"([A-Za-z])(\d+)([A-Za-z])", change.strip())
    if not m:
        raise ValueError(f"Could not parse AA change: {change}")
    return m.group(1).upper(), int(m.group(2)), m.group(3).upper()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variant_csv", help="CSV with columns: gene, aa_change")
    ap.add_argument("--pdb", required=True, help="PDB ID, e.g. 1PAH")
    ap.add_argument("--chain", default="A", help="Chain ID in the PDB structure")
    ap.add_argument("--outdir", default="results/structural", help="Output directory")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variants = []
    with open(args.variant_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            wt, pos, mut = parse_aa_change(row["aa_change"])
            variants.append((wt, pos, mut))

    # --- PyMOL script: load structure, highlight each mutated residue ---
    pymol_lines = [
        f"fetch {args.pdb}, async=0",
        "hide everything",
        "show cartoon",
        "color grey80",
    ]
    for wt, pos, mut in variants:
        sel = f"resi_{pos}"
        pymol_lines += [
            f"select {sel}, chain {args.chain} and resi {pos}",
            f"show sticks, {sel}",
            f"color red, {sel}",
            f"label {sel} and name CA, \"{wt}{pos}{mut}\"",
        ]
    pymol_lines.append(f"png {outdir}/{args.pdb}_variants_overview.png, width=1200, height=900, dpi=150, ray=1")
    (outdir / "highlight_variants.pml").write_text("\n".join(pymol_lines) + "\n")
    print(f"Wrote PyMOL script: {outdir/'highlight_variants.pml'}")
    print(f"  Run with: pymol -cq {outdir/'highlight_variants.pml'}")

    # --- FoldX individual_list.txt ---
    # Format: WTaa+Chain+Position+MUTaa;  e.g. RA261Q;
    foldx_lines = []
    for wt, pos, mut in variants:
        foldx_lines.append(f"{wt}{args.chain}{pos}{mut};")
    (outdir / "individual_list.txt").write_text("\n".join(foldx_lines) + "\n")
    print(f"Wrote FoldX mutation list: {outdir/'individual_list.txt'}")
    print("  Next steps for FoldX (run manually, not scripted here):")
    print(f"    1. foldx --command=RepairPDB --pdb={args.pdb}.pdb")
    print(f"    2. foldx --command=BuildModel --pdb={args.pdb}_Repair.pdb "
          f"--mutant-file={outdir/'individual_list.txt'}")
    print("    3. Compare Average_BuildModel_*.fxout ddG values across variants")

if __name__ == "__main__":
    main()
