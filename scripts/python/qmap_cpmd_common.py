#!/usr/bin/env python3
"""Shared utilities for the QMAP–DMSO CPMD analysis scripts.

The routines are deliberately dependency-light and stream multi-frame XYZ files,
so trajectories of several hundred megabytes do not need to be loaded into RAM.
All atom indices stored in the repository configuration are 1-based; they are
converted to 0-based indices internally.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Iterable, Sequence
import csv
import json
import math
import numpy as np

AU_TIME_FS = 0.024188843265857
DEBYE_PER_E_ANGSTROM = 4.803204712570263
LIGHT_SPEED_CM_S = 2.99792458e10

@dataclass(frozen=True)
class XYZFrame:
    index: int
    symbols: tuple[str, ...]
    coordinates: np.ndarray
    comment: str


def read_xyz_frames(path: str | Path) -> Iterator[XYZFrame]:
    """Yield frames from a standard or CPMD-style multi-frame XYZ file.

    Extra columns after x, y and z (e.g. velocities) are ignored. Blank lines
    between frames are tolerated. A malformed or truncated frame raises a
    descriptive ValueError rather than silently returning partial data.
    """
    path = Path(path)
    with path.open('r', encoding='utf-8', errors='replace') as handle:
        frame = 0
        while True:
            line = handle.readline()
            while line and not line.strip():
                line = handle.readline()
            if not line:
                break
            try:
                n_atoms = int(line.split()[0])
            except Exception as exc:
                raise ValueError(f'{path}: expected atom count before frame {frame}, got {line!r}') from exc
            comment = handle.readline().rstrip('\n')
            symbols: list[str] = []
            coords = np.empty((n_atoms, 3), dtype=float)
            for i in range(n_atoms):
                atom_line = handle.readline()
                if not atom_line:
                    raise ValueError(f'{path}: truncated frame {frame}; expected {n_atoms} atoms')
                fields = atom_line.split()
                if len(fields) < 4:
                    raise ValueError(f'{path}: invalid atom line in frame {frame}: {atom_line!r}')
                symbols.append(fields[0].capitalize())
                try:
                    coords[i] = [float(fields[1]), float(fields[2]), float(fields[3])]
                except ValueError as exc:
                    raise ValueError(f'{path}: invalid coordinates in frame {frame}: {atom_line!r}') from exc
            yield XYZFrame(frame, tuple(symbols), coords, comment)
            frame += 1


def first_xyz_frame(path: str | Path) -> XYZFrame:
    try:
        return next(read_xyz_frames(path))
    except StopIteration as exc:
        raise ValueError(f'No XYZ frame found in {path}') from exc


def load_mapping(path: str | Path) -> dict:
    with Path(path).open('r', encoding='utf-8') as handle:
        return json.load(handle)


def index0(value: int) -> int:
    return int(value) - 1


def indices0(values: Iterable[int]) -> np.ndarray:
    return np.asarray([int(v) - 1 for v in values], dtype=int)


def cell_vector(cell: float | Sequence[float]) -> np.ndarray:
    c = np.asarray(cell if np.ndim(cell) else [cell, cell, cell], dtype=float)
    if c.shape != (3,) or np.any(c <= 0):
        raise ValueError('cell must be a positive scalar or a length-3 vector')
    return c


def minimum_image(delta: np.ndarray, cell: float | Sequence[float]) -> np.ndarray:
    c = cell_vector(cell)
    return np.asarray(delta, dtype=float) - c * np.round(np.asarray(delta, dtype=float) / c)


def mic_vector(a: np.ndarray, b: np.ndarray, cell: float | Sequence[float]) -> np.ndarray:
    """Vector from point a to point b under the minimum-image convention."""
    return minimum_image(np.asarray(b) - np.asarray(a), cell)


def mic_distance(a: np.ndarray, b: np.ndarray, cell: float | Sequence[float]) -> float:
    return float(np.linalg.norm(mic_vector(a, b, cell)))


def angle_deg(a: np.ndarray, vertex: np.ndarray, c: np.ndarray,
              cell: float | Sequence[float] | None = None) -> float:
    v1 = np.asarray(a) - np.asarray(vertex)
    v2 = np.asarray(c) - np.asarray(vertex)
    if cell is not None:
        v1 = minimum_image(v1, cell)
        v2 = minimum_image(v2, cell)
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-14 or n2 < 1e-14:
        return float('nan')
    cosine = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    return math.degrees(math.acos(cosine))


def dihedral_deg(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray,
                 cell: float | Sequence[float] | None = None) -> float:
    """Signed four-point torsion in degrees, in [-180, 180]."""
    b0 = np.asarray(p1) - np.asarray(p0)
    b1 = np.asarray(p2) - np.asarray(p1)
    b2 = np.asarray(p3) - np.asarray(p2)
    if cell is not None:
        b0 = minimum_image(b0, cell)
        b1 = minimum_image(b1, cell)
        b2 = minimum_image(b2, cell)
    n = np.linalg.norm(b1)
    if n < 1e-14:
        return float('nan')
    u = b1 / n
    v = b0 - np.dot(b0, u) * u
    w = b2 - np.dot(b2, u) * u
    nv, nw = np.linalg.norm(v), np.linalg.norm(w)
    if nv < 1e-14 or nw < 1e-14:
        return float('nan')
    x = np.dot(v, w)
    y = np.dot(np.cross(u, v), w)
    return math.degrees(math.atan2(y, x))


def frame_time_ps(frame_index: int, timestep_au: float, saved_every_steps: int) -> float:
    return frame_index * timestep_au * AU_TIME_FS * saved_every_steps / 1000.0


def validate_symbols(symbols: Sequence[str], mapping: dict) -> None:
    expected = int(mapping['system']['atom_count'])
    if len(symbols) != expected:
        raise ValueError(f'Trajectory has {len(symbols)} atoms, expected {expected}')
    for symbol, ranges in mapping['element_blocks'].items():
        for start, stop in ranges:
            found = set(symbols[start - 1:stop])
            if found != {symbol}:
                raise ValueError(f'Expected {symbol} at indices {start}–{stop}, found {sorted(found)}')


def unwrap_connected(coords: np.ndarray, atom_indices_1based: Sequence[int],
                     edges_1based: Sequence[Sequence[int]],
                     cell: float | Sequence[float]) -> np.ndarray:
    """Return unwrapped coordinates for one connected molecule.

    The output follows ``atom_indices_1based``. Connectivity is traversed from
    the first atom, adding minimum-image bond vectors. It is suitable for
    molecular dipoles and internal coordinates, but does not unwrap the entire
    cluster into a unique macroscopic image.
    """
    ids = [int(i) for i in atom_indices_1based]
    position = {atom: k for k, atom in enumerate(ids)}
    adjacency = {atom: [] for atom in ids}
    for a, b in edges_1based:
        if a in adjacency and b in adjacency:
            adjacency[a].append(b)
            adjacency[b].append(a)
    out = np.full((len(ids), 3), np.nan, dtype=float)
    root = ids[0]
    out[position[root]] = coords[root - 1]
    queue = [root]
    visited = {root}
    while queue:
        a = queue.pop(0)
        for b in adjacency[a]:
            if b in visited:
                continue
            out[position[b]] = out[position[a]] + mic_vector(coords[a - 1], coords[b - 1], cell)
            visited.add(b)
            queue.append(b)
    if len(visited) != len(ids):
        missing = sorted(set(ids) - visited)
        raise ValueError(f'Connectivity does not reach atoms: {missing}')
    return out


def autocorrelation_fft(series: np.ndarray, unbiased: bool = True) -> np.ndarray:
    """Autocorrelation along axis 0 using zero-padded FFT."""
    x = np.asarray(series, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    n = x.shape[0]
    nfft = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x, n=nfft, axis=0)
    corr = np.fft.irfft(f * np.conjugate(f), n=nfft, axis=0)[:n].real
    if unbiased:
        corr /= np.arange(n, 0, -1)[:, None]
    else:
        corr /= n
    return corr.squeeze()


def read_dipole_csv(path: str | Path, prefix: str = 'cluster') -> tuple[np.ndarray, np.ndarray]:
    """Read time and Cartesian dipole columns from a CSV file.

    Accepted components include ``cluster_mu_x_D``/``mu_x_D`` and analogous
    y/z names. Time may be ``time_ps`` or is inferred from ``frame`` only if
    the caller has already converted it.
    """
    with Path(path).open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError('Dipole CSV has no header')
        names = {n.lower(): n for n in reader.fieldnames}
        def choose(candidates: Sequence[str]) -> str:
            for c in candidates:
                if c.lower() in names:
                    return names[c.lower()]
            raise ValueError(f'Missing one of columns: {candidates}')
        tcol = choose(['time_ps', 'time (ps)', 'time'])
        xcol = choose([f'{prefix}_mu_x_D', f'{prefix}_mux_D', f'{prefix}_x_D', 'mu_x_D', 'mux', 'px'])
        ycol = choose([f'{prefix}_mu_y_D', f'{prefix}_muy_D', f'{prefix}_y_D', 'mu_y_D', 'muy', 'py'])
        zcol = choose([f'{prefix}_mu_z_D', f'{prefix}_muz_D', f'{prefix}_z_D', 'mu_z_D', 'muz', 'pz'])
        times, vectors = [], []
        for row in reader:
            times.append(float(row[tcol]))
            vectors.append([float(row[xcol]), float(row[ycol]), float(row[zcol])])
    t = np.asarray(times, dtype=float)
    mu = np.asarray(vectors, dtype=float)
    if len(t) < 2:
        raise ValueError('At least two dipole samples are required')
    if not np.all(np.diff(t) > 0):
        raise ValueError('time_ps must be strictly increasing')
    return t, mu
