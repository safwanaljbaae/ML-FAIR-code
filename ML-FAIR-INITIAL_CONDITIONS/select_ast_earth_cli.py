#!/usr/bin/env python3
"""Filter asteroid candidates and generate REBOUND-style orbital elements.

This refactor preserves the behavior of the original ``Select_ast_earth.py`` while
making the workflow configurable from the command line.

Main steps
----------
1. Read an SBDB-style CSV file.
2. Filter objects by semi-major axis (default: a < 0.98 AU).
3. Write a two-column asteroid list (default: ``a`` and ``pdes``).
4. Query JPL Horizons for orbital elements at a chosen epoch.
5. Write orbital elements to ``Particles.el``.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Iterable, List

import astropy.table
import pandas as pd
from astroquery.jplhorizons import Horizons

REQUIRED_INPUT_COLUMNS = {"a", "pdes"}
DEFAULT_AST_COLUMNS = ["a", "pdes"]
PARTICLE_COLUMNS = ["a", "e", "incl", "Omega", "w", "M"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare initial-condition files from an SBDB query table."
    )
    parser.add_argument(
        "--input-csv",
        default="sbdb_query_results.csv",
        help="Input CSV file with asteroid data (default: %(default)s).",
    )
    parser.add_argument(
        "--a-max",
        type=float,
        default=0.98,
        help="Keep rows with semi-major axis a < a_max (default: %(default)s).",
    )
    parser.add_argument(
        "--sort-by",
        default="a",
        help="Column used to sort the filtered table (default: %(default)s).",
    )
    parser.add_argument(
        "--designation-column",
        default="pdes",
        help="Column passed to JPL Horizons as the object identifier (default: %(default)s).",
    )
    parser.add_argument(
        "--ast-columns",
        default=",".join(DEFAULT_AST_COLUMNS),
        help=(
            "Comma-separated columns to write to the asteroid list file "
            f"(default: {','.join(DEFAULT_AST_COLUMNS)})."
        ),
    )
    parser.add_argument(
        "--ast-output",
        default="ast_VRI.txt",
        help="Output file for the filtered asteroid list (default: %(default)s).",
    )
    parser.add_argument(
        "--particles-output",
        default="Particles.el",
        help="Output file for orbital elements (default: %(default)s).",
    )
    parser.add_argument(
        "--epoch",
        type=float,
        default=2460800.5,
        help="Julian Date epoch used in JPL Horizons queries (default: %(default)s).",
    )
    parser.add_argument(
        "--location",
        default="500@10",
        help="Observer/center location passed to JPL Horizons (default: %(default)s).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional cap on the number of filtered objects to process.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.0,
        help="Seconds to wait between Horizons queries to reduce API pressure.",
    )
    parser.add_argument(
        "--failed-ids-output",
        default="failed_horizons_ids.txt",
        help="File where failed Horizons identifiers are written (default: %(default)s).",
    )
    return parser.parse_args()


def validate_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    """Validate that required columns exist in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate.
    required : Iterable[str]
        Set of column names that must be present.

    Raises
    ------
    ValueError
        If any required columns are missing.
    """
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in input CSV: {missing}")


def normalize_ast_columns(raw_columns: str) -> List[str]:
    """Parse a comma-separated string of column names into a cleaned list.

    Parameters
    ----------
    raw_columns : str
        Comma-separated column names.

    Returns
    -------
    List[str]
        List of stripped, non-empty column names.

    Raises
    ------
    ValueError
        If no valid column names are provided.
    """
    columns = [c.strip() for c in raw_columns.split(",") if c.strip()]
    if not columns:
        raise ValueError("--ast-columns must contain at least one column name.")
    return columns


def load_and_filter_table(
    csv_path: Path,
    a_max: float,
    sort_by: str,
    limit: int | None,
) -> pd.DataFrame:
    """Load a CSV file and filter asteroids by semi-major axis.

    Parameters
    ----------
    csv_path : Path
        Path to the input CSV file.
    a_max : float
        Maximum semi-major axis cutoff.
    sort_by : str
        Column name to sort by.
    limit : int or None
        Optional limit on number of rows to return.

    Returns
    -------
    pd.DataFrame
        Filtered and sorted DataFrame.
    """
    df = pd.read_csv(csv_path, low_memory=False)
    validate_columns(df, REQUIRED_INPUT_COLUMNS | {sort_by})

    filtered = df[df["a"] < a_max].copy()
    filtered = filtered.sort_values(by=sort_by, ascending=True)

    if limit is not None:
        if limit <= 0:
            raise ValueError("--limit must be a positive integer.")
        filtered = filtered.head(limit)

    return filtered


def write_asteroid_list(df: pd.DataFrame, ast_columns: List[str], output_path: Path) -> None:
    """Write a filtered asteroid list to a text file.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing asteroid data.
    ast_columns : List[str]
        List of columns to write.
    output_path : Path
        Output file path.
    """
    validate_columns(df, ast_columns)
    df.loc[:, ast_columns].to_csv(output_path, sep=" ", header=False, index=False)


def query_horizons_elements(
    asteroid_ids: list[str],
    location: str,
    epoch: float,
    sleep_seconds: float,
) -> tuple[astropy.table.Table, list[str]]:
    """Query JPL Horizons for orbital elements of multiple asteroids.

    Parameters
    ----------
    asteroid_ids : list[str]
        List of asteroid identifiers for Horizons query.
    location : str
        Observer location code.
    epoch : float
        Julian Date for elements.
    sleep_seconds : float
        Seconds to wait between queries to avoid rate limiting.

    Returns
    -------
    tuple
        (elements_table, failed_ids) - Combined table of elements and list of failed IDs.
    """
    tables = []
    failed_ids: list[str] = []

    for idx, asteroid_id in enumerate(asteroid_ids, start=1):
        try:
            result = Horizons(id=str(asteroid_id), location=location, epochs=epoch).elements()
            tables.append(result)
        except Exception as exc:  # noqa: BLE001
            failed_ids.append(str(asteroid_id))
            print(f"[warning] Horizons query failed for {asteroid_id}: {exc}", file=sys.stderr)

        if sleep_seconds > 0 and idx < len(asteroid_ids):
            time.sleep(sleep_seconds)


    if not tables:
        raise RuntimeError("All Horizons queries failed; no Particles.el file can be created.")

    combined = astropy.table.vstack(tables)
    return combined, failed_ids


def _as_float(value) -> float:
    """Convert an astropy scalar/quantity-like value to float."""
    if hasattr(value, "value"):
        return float(value.value)
    return float(value)


def write_particles_file(elements: astropy.table.Table, output_path: Path) -> None:
    """Write orbital elements to a REBOUND-compatible Particles.el file.

    Parameters
    ----------
    elements : astropy.table.Table
        Table containing orbital elements (a, e, incl, Omega, w, M).
    output_path : Path
        Output file path.
    """
    with output_path.open("w", encoding="utf-8") as handle:
        for row in elements:
            values = [_as_float(row[col]) for col in PARTICLE_COLUMNS]
            handle.write(
                f"{values[0]:.15f} {values[1]:.16f} {values[2]:.14f} "
                f"{values[3]:.15f} {values[4]:.15f} {values[5]:.15f}\n"
            )


def write_failed_ids(failed_ids: list[str], output_path: Path) -> None:
    """Write list of failed asteroid IDs to a file.

    Parameters
    ----------
    failed_ids : list[str]
        List of asteroid IDs that failed to query.
    output_path : Path
        Output file path.
    """
    if not failed_ids:
        return
    with output_path.open("w", encoding="utf-8") as handle:
        for asteroid_id in failed_ids:
            handle.write(f"{asteroid_id}\n")


def main() -> int:
    args = parse_args()
    start = time.time()

    input_csv = Path(args.input_csv)
    ast_output = Path(args.ast_output)
    particles_output = Path(args.particles_output)
    failed_ids_output = Path(args.failed_ids_output)
    ast_columns = normalize_ast_columns(args.ast_columns)

    try:
        filtered = load_and_filter_table(
            csv_path=input_csv,
            a_max=args.a_max,
            sort_by=args.sort_by,
            limit=args.limit,
        )
        validate_columns(filtered, {args.designation_column})
        if filtered.empty:
            raise RuntimeError("No asteroid rows matched the requested filtering criteria.")

        write_asteroid_list(filtered, ast_columns, ast_output)

        asteroid_ids = filtered[args.designation_column].astype(str).tolist()
        elements, failed_ids = query_horizons_elements(
            asteroid_ids=asteroid_ids,
            location=args.location,
            epoch=args.epoch,
            sleep_seconds=args.sleep,
        )
        write_particles_file(elements, particles_output)
        write_failed_ids(failed_ids, failed_ids_output)

    except Exception as exc:  # noqa: BLE001
        print(f"[error] {exc}", file=sys.stderr)
        return 1

    elapsed = time.time() - start
    print(f"Filtered asteroid count: {len(filtered)}")
    print(f"Asteroid list written to: {ast_output}")
    print(f"Orbital elements written to: {particles_output}")
    if failed_ids:
        print(f"Failed Horizons IDs written to: {failed_ids_output}")
        print(f"Horizons failures: {len(failed_ids)}")
    print(f"Execution time: {elapsed:.2f} seconds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
