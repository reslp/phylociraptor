#!/usr/bin/env python3
"""
Config validator for phylociraptor.
Validates the data CSV file and the config.yaml file.
"""

import sys
import os
import re
import csv
import io

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML is required but not installed.")
    sys.exit(1)


# ─── ANSI colors ────────────────────────────────────────────────────────────────
RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
RESET = "\033[0m"


def print_ok(msg):
    print(f"  {GREEN}[OK]{RESET} {msg}")


def print_warn(msg):
    print(f"  {YELLOW}[WARNING]{RESET} {msg}")


def print_err(msg):
    print(f"  {RED}[ERROR]{RESET} {msg}")


# ─── CSV validation ─────────────────────────────────────────────────────────────

def validate_csv(csv_path):
    """Validate the data CSV file. Returns (errors, warnings) counts."""
    errors = 0
    warnings = 0

    print(f"\nValidating data CSV file: {csv_path}")
    print("-" * 60)

    if not os.path.isfile(csv_path):
        print_err(f"File not found: {csv_path}")
        return 1, 0

    # Read raw bytes to check line endings
    with open(csv_path, "rb") as f:
        raw = f.read()

    # 1) Check line endings
    has_crlf = b"\r\n" in raw
    has_cr_only = bool(re.search(b"\r(?!\n)", raw))
    has_lf = b"\n" in raw

    if has_cr_only:
        print_err("File uses old Mac-style line endings (\\r). Convert to Unix (\\n).")
        errors += 1
    elif has_crlf:
        print_err("File uses Windows-style line endings (\\r\\n). Convert to Unix (\\n).")
        errors += 1
    elif has_lf:
        print_ok("Line endings are Unix-style (\\n).")
    else:
        print_warn("File appears to be a single line (no line endings detected).")
        warnings += 1

    # Decode for further checks, stripping an optional UTF-8 BOM
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        print_err("File is not valid UTF-8.")
        return errors + 1, warnings

    lines = text.splitlines()
    if len(lines) == 0:
        print_err("File is empty.")
        return errors + 1, warnings

    # 2) Check delimiter – the header must be comma-separated
    header_line = lines[0]
    if "\t" in header_line and "," not in header_line:
        print_err("File appears to be tab-delimited instead of comma-delimited.")
        errors += 1
    elif ";" in header_line and "," not in header_line:
        print_err("File appears to be semicolon-delimited instead of comma-delimited.")
        errors += 1
    elif "," in header_line:
        print_ok("File is comma-delimited.")
    else:
        print_warn("Could not detect delimiter (header has no commas, tabs, or semicolons).")
        warnings += 1

    # Parse rows once with csv.reader so malformed row widths can be reported cleanly
    parsed_rows = list(csv.reader(io.StringIO(text)))
    if len(parsed_rows) == 0:
        print_err("No header line detected.")
        return errors + 1, warnings

    fieldnames = [field.strip() for field in parsed_rows[0]]

    # 3) Check header line
    if len(fieldnames) == 0:
        print_err("No header line detected.")
        return errors + 1, warnings

    print_ok(f"Header detected with {len(fieldnames)} columns: {', '.join(fieldnames)}")

    # Check required columns
    required_cols = ["species", "web_local"]
    for col in required_cols:
        if col not in fieldnames:
            print_err(f"Required column '{col}' not found in header.")
            errors += 1

    has_species = "species" in fieldnames
    has_web_local = "web_local" in fieldnames
    has_mode = "mode" in fieldnames

    if has_mode:
        print_ok("Optional column 'mode' is present.")

    # Check row widths and then build dict rows with normalized field names
    expected_columns = len(fieldnames)
    rows = []
    for i, row_values in enumerate(parsed_rows[1:], start=2):
        if len(row_values) != expected_columns:
            print_err(
                f"Line {i}: Expected {expected_columns} columns based on header, found {len(row_values)}."
            )
            errors += 1
        padded_row_values = row_values + [""] * max(0, expected_columns - len(row_values))
        row = dict(zip(fieldnames, padded_row_values[:expected_columns]))
        rows.append(row)

    # 4 & 5) Check data rows
    # Allowed characters: alphanumeric, underscore, space, dot, hyphen, equals, slash (for paths)
    special_char_pattern = re.compile(r"[^A-Za-z0-9_ .\-=\/]")
    valid_modes = {"genome", "transcriptome", "proteome", ""}

    if len(rows) == 0:
        print_err("File has a header but no data rows.")
        errors += 1
        return errors, warnings

    print_ok(f"File has {len(rows)} data rows.")

    for i, row in enumerate(rows, start=2):  # line 2 onwards (1-indexed, line 1 = header)
        # Check species column
        if has_species:
            species_val = (row.get("species", "") or "").strip()
            match = special_char_pattern.search(species_val)
            if match:
                print_err(f"Line {i}: Column 'species' contains special character "
                          f"'{match.group()}' in value: {species_val}")
                errors += 1

        # Check web_local column
        if has_web_local:
            web_local_val = (row.get("web_local", "") or "").strip()
            match = special_char_pattern.search(web_local_val)
            if match:
                print_err(f"Line {i}: Column 'web_local' contains special character "
                          f"'{match.group()}' in value: {web_local_val}")
                errors += 1

        # Check mode column (if present)
        if has_mode:
            mode_val = (row.get("mode", "") or "").strip()
            if mode_val not in valid_modes:
                print_err(f"Line {i}: Column 'mode' has invalid value '{mode_val}'. "
                          f"Allowed: genome, transcriptome, proteome, or empty.")
                errors += 1

    if errors == 0:
        print_ok("All data rows passed validation.")

    return errors, warnings


# ─── YAML config validation ─────────────────────────────────────────────────────

def validate_config_yaml(yaml_path):
    """Validate the config.yaml file. Returns (errors, warnings) counts."""
    errors = 0
    warnings = 0

    print(f"\nValidating config file: {yaml_path}")
    print("-" * 60)

    if not os.path.isfile(yaml_path):
        print_err(f"File not found: {yaml_path}")
        return 1, 0

    # Try to parse YAML
    try:
        with open(yaml_path) as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        print_err(f"YAML parsing error: {e}")
        return 1, 0

    if config is None:
        print_err("Config file is empty or not valid YAML.")
        return 1, 0

    if not isinstance(config, dict):
        print_err("Config file top-level structure must be a mapping (dict).")
        return 1, 0

    print_ok("YAML file is syntactically valid.")

    # --- Check top-level keys and types ---

    # species (path to CSV) – required, string
    if "species" not in config:
        print_err("Missing required key: 'species' (path to data CSV file).")
        errors += 1
    else:
        if not isinstance(config["species"], str):
            print_err("'species' must be a string (path to CSV file).")
            errors += 1
        else:
            print_ok(f"'species' is set to: {config['species']}")

    # email – optional, string
    if "email" in config:
        if not isinstance(config["email"], str):
            print_err("'email' must be a string.")
            errors += 1
        else:
            print_ok(f"'email' is set.")

    # concurrency – optional, integer
    if "concurrency" in config:
        if not isinstance(config["concurrency"], int):
            print_err("'concurrency' must be an integer.")
            errors += 1
        else:
            print_ok(f"'concurrency' = {config['concurrency']}")

    # seed – optional, integer or None
    if "seed" in config:
        if config["seed"] is not None and not isinstance(config["seed"], int):
            print_err("'seed' must be an integer or empty (null).")
            errors += 1
        else:
            print_ok(f"'seed' = {config['seed']}")

    # --- orthology section ---
    if "orthology" in config:
        orth = config["orthology"]
        if not isinstance(orth, dict):
            print_err("'orthology' must be a mapping.")
            errors += 1
        else:
            print_ok("'orthology' section present.")
            if "method" in orth and orth["method"] not in ("busco", "orthofinder"):
                print_warn(f"'orthology.method' has unusual value: {orth['method']}")
                warnings += 1
            if "threads" in orth and not isinstance(orth["threads"], int):
                print_err("'orthology.threads' must be an integer.")
                errors += 1
            if "busco_options" in orth:
                bo = orth["busco_options"]
                if not isinstance(bo, dict):
                    print_err("'orthology.busco_options' must be a mapping.")
                    errors += 1
                else:
                    if "set" in bo and not isinstance(bo["set"], str):
                        print_err("'orthology.busco_options.set' must be a string.")
                        errors += 1
                    if "version" in bo and not isinstance(bo["version"], str):
                        print_err("'orthology.busco_options.version' must be a string.")
                        errors += 1
                    if "mode" in bo:
                        valid_busco_modes = ("genome", "transcriptome", "proteins")
                        if bo["mode"] not in valid_busco_modes:
                            print_warn(f"'orthology.busco_options.mode' = '{bo['mode']}' "
                                       f"is not one of {valid_busco_modes}.")
                            warnings += 1

    # --- filtering section ---
    if "filtering" in config:
        filt = config["filtering"]
        if not isinstance(filt, dict):
            print_err("'filtering' must be a mapping.")
            errors += 1
        else:
            print_ok("'filtering' section present.")
            if "cutoff" in filt:
                if not isinstance(filt["cutoff"], (int, float)):
                    print_err("'filtering.cutoff' must be a number.")
                    errors += 1
                elif not (0 <= filt["cutoff"] <= 1):
                    print_warn(f"'filtering.cutoff' = {filt['cutoff']} is outside [0, 1].")
                    warnings += 1
            if "minsp" in filt and not isinstance(filt["minsp"], int):
                print_err("'filtering.minsp' must be an integer.")
                errors += 1
            if "seq_type" in filt:
                if filt["seq_type"] not in ("aa", "nu"):
                    print_warn(f"'filtering.seq_type' = '{filt['seq_type']}' – expected 'aa' or 'nu'.")
                    warnings += 1
            if "dupseq" in filt:
                if filt["dupseq"] not in ("persample", "perfile"):
                    print_warn(f"'filtering.dupseq' = '{filt['dupseq']}' – expected 'persample' or 'perfile'.")
                    warnings += 1

    # --- alignment section ---
    if "alignment" in config:
        aln = config["alignment"]
        if not isinstance(aln, dict):
            print_err("'alignment' must be a mapping.")
            errors += 1
        else:
            print_ok("'alignment' section present.")
            valid_aln_methods = {"mafft", "clustalo", "muscle", "tcoffee", "prank"}
            if "method" in aln:
                methods = aln["method"] if isinstance(aln["method"], list) else [aln["method"]]
                for m in methods:
                    if m not in valid_aln_methods:
                        print_warn(f"'alignment.method' contains unknown method: '{m}'")
                        warnings += 1
            if "threads" in aln and not isinstance(aln["threads"], int):
                print_err("'alignment.threads' must be an integer.")
                errors += 1

    # --- trimming section ---
    if "trimming" in config:
        trim = config["trimming"]
        if not isinstance(trim, dict):
            print_err("'trimming' must be a mapping.")
            errors += 1
        else:
            print_ok("'trimming' section present.")
            valid_trim_methods = {"trimal", "aliscore", "bmge", "clipkit"}
            if "method" in trim:
                methods = trim["method"] if isinstance(trim["method"], list) else [trim["method"]]
                for m in methods:
                    if m not in valid_trim_methods:
                        print_warn(f"'trimming.method' contains unknown method: '{m}'")
                        warnings += 1
            if "minsp" in trim and not isinstance(trim["minsp"], int):
                print_err("'trimming.minsp' must be an integer.")
                errors += 1
            if "min_parsimony_sites" in trim and not isinstance(trim["min_parsimony_sites"], int):
                print_err("'trimming.min_parsimony_sites' must be an integer.")
                errors += 1
            if "max_rcv_score" in trim:
                if not isinstance(trim["max_rcv_score"], (int, float)):
                    print_err("'trimming.max_rcv_score' must be a number.")
                    errors += 1
                elif not (0 <= trim["max_rcv_score"] <= 1):
                    print_warn(f"'trimming.max_rcv_score' = {trim['max_rcv_score']} is outside [0, 1].")
                    warnings += 1

    # --- modeltest section ---
    if "modeltest" in config:
        mt = config["modeltest"]
        if not isinstance(mt, dict):
            print_err("'modeltest' must be a mapping.")
            errors += 1
        else:
            print_ok("'modeltest' section present.")
            if "threads" in mt and not isinstance(mt["threads"], int):
                print_err("'modeltest.threads' must be an integer.")
                errors += 1
            if "bootstrap" in mt and not isinstance(mt["bootstrap"], int):
                print_err("'modeltest.bootstrap' must be an integer.")
                errors += 1

    # --- genetree_filtering section ---
    if "genetree_filtering" in config:
        gf = config["genetree_filtering"]
        if not isinstance(gf, dict):
            print_err("'genetree_filtering' must be a mapping.")
            errors += 1
        else:
            print_ok("'genetree_filtering' section present.")
            if "bootstrap_cutoff" in gf:
                cutoffs = gf["bootstrap_cutoff"]
                if isinstance(cutoffs, list):
                    for c in cutoffs:
                        if not isinstance(c, (int, float)):
                            print_err(f"'genetree_filtering.bootstrap_cutoff' contains non-numeric value: {c}")
                            errors += 1
                elif not isinstance(cutoffs, (int, float)):
                    print_err("'genetree_filtering.bootstrap_cutoff' must be a number or list of numbers.")
                    errors += 1

    # --- speciestree section ---
    if "speciestree" in config:
        st = config["speciestree"]
        if not isinstance(st, dict):
            print_err("'speciestree' must be a mapping.")
            errors += 1
        else:
            print_ok("'speciestree' section present.")

    # --- mltree section ---
    if "mltree" in config:
        ml = config["mltree"]
        if not isinstance(ml, dict):
            print_err("'mltree' must be a mapping.")
            errors += 1
        else:
            print_ok("'mltree' section present.")
            if "threads" in ml and not isinstance(ml["threads"], dict):
                print_err("'mltree.threads' must be a mapping (method: value).")
                errors += 1
            if "bootstrap" in ml and not isinstance(ml["bootstrap"], dict):
                print_err("'mltree.bootstrap' must be a mapping (method: value).")
                errors += 1

    # --- njtree section ---
    if "njtree" in config:
        nj = config["njtree"]
        if not isinstance(nj, dict):
            print_err("'njtree' must be a mapping.")
            errors += 1
        else:
            print_ok("'njtree' section present.")

    # --- bitree section ---
    if "bitree" in config:
        bi = config["bitree"]
        if not isinstance(bi, dict):
            print_err("'bitree' must be a mapping.")
            errors += 1
        else:
            print_ok("'bitree' section present.")

    return errors, warnings


# ─── Main entry point ───────────────────────────────────────────────────────────

def main(config_file="data/config.yaml"):
    """Run all validation checks."""
    print("=" * 60)
    print("phylociraptor config validator")
    print("=" * 60)

    total_errors = 0
    total_warnings = 0

    # First validate the config YAML to extract the CSV path
    yaml_errors, yaml_warnings = validate_config_yaml(config_file)
    total_errors += yaml_errors
    total_warnings += yaml_warnings

    # Try to get the CSV path from the config
    csv_path = None
    try:
        with open(config_file) as f:
            config = yaml.safe_load(f)
        if isinstance(config, dict) and "species" in config:
            csv_path = config["species"]
    except Exception:
        pass

    if csv_path:
        csv_errors, csv_warnings = validate_csv(csv_path)
        total_errors += csv_errors
        total_warnings += csv_warnings
    else:
        print(f"\n{RED}Skipping CSV validation: could not determine CSV path from config.{RESET}")

    # Summary
    print("\n" + "=" * 60)
    print("Validation summary")
    print("=" * 60)
    if total_errors == 0 and total_warnings == 0:
        print(f"  {GREEN}All checks passed!{RESET}")
    else:
        if total_warnings > 0:
            print(f"  {YELLOW}Warnings: {total_warnings}{RESET}")
        if total_errors > 0:
            print(f"  {RED}Errors:   {total_errors}{RESET}")
    print()

    return total_errors


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Validate phylociraptor config and data files.")
    parser.add_argument("--config-file", default="data/config.yaml",
                        help="Path to the config.yaml file (default: data/config.yaml)")
    pargs = parser.parse_args()
    sys.exit(0 if main(pargs.config_file) == 0 else 1)
