#!/usr/bin/env python3
"""Split asteroid input files into simulation-ready batch directories."""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path


DEFAULT_AST_FILE = "ast_VRI.txt"
DEFAULT_PARTICLES_FILE = "Particles.el"
DEFAULT_TEMPLATE_FILE = "Earth_inner_ensemble.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split ast_VRI.txt and Particles.el into batch directories."
    )
    parser.add_argument(
        "--ast-file",
        default=DEFAULT_AST_FILE,
        help="Path to the asteroid list file (default: %(default)s).",
    )
    parser.add_argument(
        "--particles-file",
        default=DEFAULT_PARTICLES_FILE,
        help="Path to the orbital-elements file (default: %(default)s).",
    )
    parser.add_argument(
        "--template-script",
        default=DEFAULT_TEMPLATE_FILE,
        help="Simulation script copied into each output directory (default: %(default)s).",
    )
    parser.add_argument(
        "--lines-per-batch",
        type=int,
        default=100,
        help="Maximum number of objects per split directory (default: %(default)s).",
    )
    parser.add_argument(
        "--output-prefix",
        default="split_",
        help="Prefix used for generated directories (default: %(default)s).",
    )
    parser.add_argument(
        "--copied-script-name",
        default="Earth_inner_ensemble.py",
        help="Filename to use for the copied template script inside each split directory.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if ast-file and particles-file do not contain the same number of lines.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete an existing output directory before re-creating it.",
    )
    return parser.parse_args()


def read_lines(path: Path) -> list[str]:
    """Read all lines from a file.

    Parameters
    ----------
    path : Path
        Path to the file to read.

    Returns
    -------
    list[str]
        List of lines (including newlines).
    """
    with path.open("r", encoding="utf-8") as handle:
        return handle.readlines()


def ensure_inputs(ast_lines: list[str], particle_lines: list[str], strict: bool) -> None:
    """Validate that asteroid and particle input files are properly formatted.

    Parameters
    ----------
    ast_lines : list[str]
        Lines from the asteroid list file.
    particle_lines : list[str]
        Lines from the particles file.
    strict : bool
        If True, requires both files to have the same number of lines.

    Raises
    ------
    ValueError
        If validation fails.
    """
    if not ast_lines:
        raise ValueError("The asteroid list file is empty.")
    if not particle_lines:
        raise ValueError("The particles file is empty.")
    if strict and len(ast_lines) != len(particle_lines):
        raise ValueError(
            "Input files do not contain the same number of lines. "
            f"ast-file={len(ast_lines)}, particles-file={len(particle_lines)}"
        )


def prepare_directory(path: Path, overwrite: bool) -> None:
    """Create a directory, optionally removing it first if it exists.

    Parameters
    ----------
    path : Path
        Directory path to create.
    overwrite : bool
        If True, removes existing directory before creating.
    """
    if path.exists() and overwrite:
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def split_into_dirs(
    ast_file: Path,
    particles_file: Path,
    template_script: Path,
    lines_per_batch: int,
    output_prefix: str,
    copied_script_name: str,
    strict: bool,
    overwrite: bool,
) -> int:
    """Split asteroid and particle files into batch directories for parallel simulations.

    Each batch directory receives a subset of the asteroid list, particle elements,
    and a copy of the template simulation script.

    Parameters
    ----------
    ast_file : Path
        Path to the asteroid list file.
    particles_file : Path
        Path to the particles file with orbital elements.
    template_script : Path
        Path to the template simulation script.
    lines_per_batch : int
        Maximum number of lines per batch directory.
    output_prefix : str
        Prefix for batch directory names.
    copied_script_name : str
        Filename for the copied template script.
    strict : bool
        If True, requires equal line counts in input files.
    overwrite : bool
        If True, overwrites existing directories.

    Returns
    -------
    int
        Number of batch directories created.
    """
    if lines_per_batch <= 0:
        raise ValueError("--lines-per-batch must be a positive integer.")
    if not template_script.exists():
        raise FileNotFoundError(f"Template script not found: {template_script}")

    ast_lines = read_lines(ast_file)
    particle_lines = read_lines(particles_file)
    ensure_inputs(ast_lines, particle_lines, strict=strict)

    n_batches = max(
        math.ceil(len(ast_lines) / lines_per_batch),
        math.ceil(len(particle_lines) / lines_per_batch),
    )

    for batch_index in range(n_batches):
        start = batch_index * lines_per_batch
        stop = (batch_index + 1) * lines_per_batch

        directory = Path(f"{output_prefix}{batch_index + 1}")
        prepare_directory(directory, overwrite=overwrite)

        (directory / ast_file.name).write_text(
            "".join(ast_lines[start:stop]),
            encoding="utf-8",
        )
        (directory / particles_file.name).write_text(
            "".join(particle_lines[start:stop]),
            encoding="utf-8",
        )
        shutil.copy(template_script, directory / copied_script_name)

    return n_batches


def main() -> int:
    args = parse_args()

    try:
        n_batches = split_into_dirs(
            ast_file=Path(args.ast_file),
            particles_file=Path(args.particles_file),
            template_script=Path(args.template_script),
            lines_per_batch=args.lines_per_batch,
            output_prefix=args.output_prefix,
            copied_script_name=args.copied_script_name,
            strict=args.strict,
            overwrite=args.overwrite,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"[error] {exc}")
        return 1

    print(
        f"Created {n_batches} directories containing {Path(args.ast_file).name}, "
        f"{Path(args.particles_file).name}, and {args.copied_script_name}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
