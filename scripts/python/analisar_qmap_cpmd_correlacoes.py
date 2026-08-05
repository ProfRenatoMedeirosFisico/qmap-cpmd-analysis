#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analisar_qmap_cpmd_correlacoes.py

Objetivo:
  Ler trajetória .xyz multi-frame ou .pdb, identificar automaticamente QMAP e DMSO
  em um sistema QMAP + DMSO, calcular:
    1) dipolo do QMAP
    2) menor distância O(QMAP)–O(DMSO) por frame
    3) dihedral C2-C1-C3-C4, se os índices forem informados
  e gerar figuras de correlação em 600 dpi.

Uso básico:
  python analisar_qmap_cpmd_correlacoes.py --traj traj.xyz --charges cargas.txt

Com dihedral:
  python analisar_qmap_cpmd_correlacoes.py --traj traj.xyz --charges cargas.txt --dihedral 12 8 15 16

Índices:
  Por padrão, use índices iniciando em 1, como em VMD/Avogadro.
  Para usar índices iniciando em 0, adicione: --zero-based

Formatos aceitos para cargas:
  1) Uma carga por linha, na ordem dos átomos:
       -0.12
        0.08
        ...
  2) Duas colunas: indice carga
       1 -0.12
       2  0.08
       ...
  3) CSV/TSV por frame com colunas:
       frame,atom,charge
     ou matriz com nframes linhas e natoms colunas.

Observação crítica:
  A identificação automática separa DMSO por átomos de enxofre (S) e seus vizinhos próximos.
  Isso funciona para QMAP + DMSO em XYZ/PDB comum, mas deve ser validado visualmente.
  O dihedral C2-C1-C3-C4 NÃO deve ser inventado automaticamente para publicação:
  informe os quatro índices corretos com --dihedral.
"""

import argparse
import csv
import os
import sys

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Leitura de trajetórias
# -----------------------------

def read_xyz(path):
    frames = []
    symbols_ref = None

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.rstrip() for ln in f if ln.strip()]

    i = 0
    while i < len(lines):
        try:
            n = int(lines[i].split()[0])
        except Exception:
            raise ValueError(f"Erro lendo XYZ perto da linha {i+1}: esperado número de átomos.")

        block = lines[i + 2:i + 2 + n]
        if len(block) != n:
            raise ValueError("XYZ incompleto: número de linhas menor que o número de átomos.")

        symbols = []
        coords = []
        for ln in block:
            parts = ln.split()
            if len(parts) < 4:
                raise ValueError(f"Linha XYZ inválida: {ln}")
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

        if symbols_ref is None:
            symbols_ref = symbols
        elif symbols != symbols_ref:
            print("AVISO: símbolos mudaram entre frames. Verifique a trajetória.", file=sys.stderr)

        frames.append(np.array(coords, dtype=float))
        i += n + 2

    return symbols_ref, frames


def read_pdb(path):
    # Leitor simples para PDB com MODEL/ENDMDL. Se não houver MODEL, lê como um frame.
    frames = []
    symbols_ref = []
    current_coords = []
    current_symbols = []

    def infer_symbol(atom_name, element_field):
        e = element_field.strip()
        if e:
            return e.capitalize()
        name = atom_name.strip()
        letters = ''.join([c for c in name if c.isalpha()])
        if not letters:
            return "X"
        if len(letters) >= 2 and letters[:2].capitalize() in ["Cl", "Br", "Bi", "Si"]:
            return letters[:2].capitalize()
        return letters[0].upper()

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            rec = ln[:6].strip()
            if rec == "MODEL":
                current_coords = []
                current_symbols = []
            elif rec in ("ATOM", "HETATM"):
                atom_name = ln[12:16]
                element_field = ln[76:78] if len(ln) >= 78 else ""
                sym = infer_symbol(atom_name, element_field)
                x = float(ln[30:38])
                y = float(ln[38:46])
                z = float(ln[46:54])
                current_symbols.append(sym)
                current_coords.append([x, y, z])
            elif rec == "ENDMDL":
                if current_coords:
                    frames.append(np.array(current_coords, dtype=float))
                    if not symbols_ref:
                        symbols_ref = current_symbols[:]
                    current_coords = []
                    current_symbols = []

    if current_coords:
        frames.append(np.array(current_coords, dtype=float))
        if not symbols_ref:
            symbols_ref = current_symbols[:]

    if not frames:
        raise ValueError("Nenhum frame encontrado no PDB.")

    return symbols_ref, frames


def read_trajectory(path):
    ext = os.path.splitext(path.lower())[1]
    if ext == ".xyz":
        return read_xyz(path)
    if ext == ".pdb":
        return read_pdb(path)
    raise ValueError("Formato não reconhecido. Use .xyz ou .pdb.")


# -----------------------------
# Cargas
# -----------------------------

def load_charges(path, n_atoms, n_frames):
    if path is None:
        print("AVISO: arquivo de cargas não fornecido. O dipolo não será calculado.", file=sys.stderr)
        return None

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        raw = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]

    if not raw:
        raise ValueError("Arquivo de cargas vazio.")

    header = [h.strip().lower() for h in raw[0].replace(";", ",").replace("\t", ",").split(",")]
    if {"frame", "atom", "charge"}.issubset(set(header)):
        idx_frame = header.index("frame")
        idx_atom = header.index("atom")
        idx_charge = header.index("charge")
        arr = np.full((n_frames, n_atoms), np.nan, dtype=float)
        for ln in raw[1:]:
            parts = [p.strip() for p in ln.replace(";", ",").replace("\t", ",").split(",")]
            fr = int(parts[idx_frame])
            at = int(parts[idx_atom])
            q = float(parts[idx_charge])
            if at >= 1 and at <= n_atoms:
                at -= 1
            arr[fr, at] = q
        if np.isnan(arr).any():
            raise ValueError("Arquivo frame,atom,charge possui valores faltantes.")
        return arr

    matrix = []
    ok_matrix = True
    for ln in raw:
        parts = ln.replace(";", " ").replace(",", " ").split()
        try:
            vals = [float(x) for x in parts]
        except Exception:
            ok_matrix = False
            break
        matrix.append(vals)

    if ok_matrix:
        if len(matrix) == n_atoms and all(len(row) == 1 for row in matrix):
            q = np.array([row[0] for row in matrix], dtype=float)
            return np.tile(q, (n_frames, 1))

        if len(matrix) == n_atoms and all(len(row) >= 2 for row in matrix):
            col1 = [int(row[0]) for row in matrix]
            if set(col1) == set(range(1, n_atoms + 1)) or set(col1) == set(range(n_atoms)):
                q = np.zeros(n_atoms, dtype=float)
                one_based = set(col1) == set(range(1, n_atoms + 1))
                for row in matrix:
                    idx = int(row[0]) - 1 if one_based else int(row[0])
                    q[idx] = float(row[1])
                return np.tile(q, (n_frames, 1))

        if len(matrix) == n_frames and all(len(row) == n_atoms for row in matrix):
            return np.array(matrix, dtype=float)

    raise ValueError(
        "Não consegui interpretar o arquivo de cargas. Use uma carga por linha, "
        "duas colunas indice/carga, matriz nframes x natoms, ou CSV frame,atom,charge."
    )


# -----------------------------
# Geometria
# -----------------------------

def distance(a, b):
    return float(np.linalg.norm(a - b))


def compute_dihedral(p0, p1, p2, p3):
    # retorna graus no intervalo [-180, 180]
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b1_norm = np.linalg.norm(b1)
    if b1_norm < 1e-12:
        return np.nan
    b1 = b1 / b1_norm

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    if np.linalg.norm(v) < 1e-12 or np.linalg.norm(w) < 1e-12:
        return np.nan

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def dipole_debye(coords, charges):
    # coords em Å, carga em e; 1 e·Å = 4.80320427 Debye
    mu_vec = np.sum(coords * charges[:, None], axis=0)
    return float(np.linalg.norm(mu_vec) * 4.80320427), mu_vec * 4.80320427


# -----------------------------
# Identificação automática QMAP/DMSO
# -----------------------------

def identify_dmso_and_qmap(symbols, coords0):
    symbols_clean = [s.capitalize() for s in symbols]
    s_indices = [i for i, s in enumerate(symbols_clean) if s == "S"]

    dmso_atoms = set()
    dmso_oxygens = []
    dmso_groups = []

    for si in s_indices:
        o_candidates = [i for i, s in enumerate(symbols_clean) if s == "O"]
        if not o_candidates:
            continue
        o_near = min(o_candidates, key=lambda i: distance(coords0[si], coords0[i]))

        c_candidates = [i for i, s in enumerate(symbols_clean) if s == "C"]
        c_sorted = sorted(c_candidates, key=lambda i: distance(coords0[si], coords0[i]))
        c_near = c_sorted[:2]

        group = {si, o_near, *c_near}

        h_candidates = [i for i, s in enumerate(symbols_clean) if s == "H"]
        for ci in c_near:
            h_sorted = sorted(h_candidates, key=lambda i: distance(coords0[ci], coords0[i]))
            for hi in h_sorted[:3]:
                if distance(coords0[ci], coords0[hi]) < 1.35:
                    group.add(hi)

        dmso_groups.append(sorted(group))
        dmso_atoms.update(group)
        dmso_oxygens.append(o_near)

    all_atoms = set(range(len(symbols)))
    qmap_atoms = sorted(all_atoms - dmso_atoms)
    qmap_oxygens = [i for i in qmap_atoms if symbols_clean[i] == "O"]

    phenolic_o = None
    if qmap_oxygens:
        h_qmap = [i for i in qmap_atoms if symbols_clean[i] == "H"]
        if h_qmap:
            phenolic_o = min(
                qmap_oxygens,
                key=lambda oi: min(distance(coords0[oi], coords0[hi]) for hi in h_qmap)
            )
        else:
            phenolic_o = qmap_oxygens[0]

    return {
        "dmso_groups": dmso_groups,
        "dmso_atoms": sorted(dmso_atoms),
        "dmso_oxygens": sorted(set(dmso_oxygens)),
        "qmap_atoms": qmap_atoms,
        "qmap_oxygens": qmap_oxygens,
        "phenolic_o": phenolic_o,
    }


# -----------------------------
# Plot
# -----------------------------

def moving_average_xy(x, y, window=25):
    x = np.asarray(x)
    y = np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    if len(y) < window:
        return x, y
    kernel = np.ones(window) / window
    ys = np.convolve(y, kernel, mode="valid")
    xs = x[window // 2: window // 2 + len(ys)]
    return xs, ys


def save_scatter(x, y, xlabel, ylabel, title, outpng, window=25):
    plt.figure(figsize=(7.0, 5.2))
    plt.scatter(x, y, s=14, alpha=0.55, edgecolors="none")
    xs, ys = moving_average_xy(x, y, window=window)
    if len(xs) > 2:
        plt.plot(xs, ys, linewidth=2.0)
    plt.xlabel(xlabel, fontsize=13)
    plt.ylabel(ylabel, fontsize=13)
    plt.title(title, fontsize=13)
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(outpng, dpi=600)
    plt.close()


def save_timeseries(t, series, ylabel, title, outpng):
    plt.figure(figsize=(7.5, 4.5))
    plt.plot(t, series, linewidth=1.4)
    plt.xlabel("Frame", fontsize=13)
    plt.ylabel(ylabel, fontsize=13)
    plt.title(title, fontsize=13)
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(outpng, dpi=600)
    plt.close()


# -----------------------------
# Principal
# -----------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--traj", required=True, help="Trajetória .xyz ou .pdb")
    parser.add_argument("--charges", default=None, help="Arquivo de cargas fixas ou por frame")
    parser.add_argument("--dihedral", nargs=4, type=int, default=None,
                        help="Quatro índices para o dihedral C2 C1 C3 C4")
    parser.add_argument("--zero-based", action="store_true",
                        help="Usar índices começando em 0 para --dihedral")
    parser.add_argument("--outdir", default="figuras_qmap_correlacoes",
                        help="Diretório de saída")
    parser.add_argument("--smooth-window", type=int, default=25,
                        help="Janela da média móvel nos gráficos")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    symbols, frames = read_trajectory(args.traj)
    n_atoms = len(symbols)
    n_frames = len(frames)

    print(f"Trajetória lida: {args.traj}")
    print(f"Número de átomos: {n_atoms}")
    print(f"Número de frames: {n_frames}")

    ident = identify_dmso_and_qmap(symbols, frames[0])

    print("\nIdentificação automática:")
    print(f"  DMSO encontrados pelo número de S: {len(ident['dmso_groups'])}")
    print(f"  Átomos atribuídos ao DMSO: {len(ident['dmso_atoms'])}")
    print(f"  Átomos atribuídos ao QMAP: {len(ident['qmap_atoms'])}")
    print(f"  O(DMSO): {[i+1 for i in ident['dmso_oxygens']]}")
    print(f"  O(QMAP): {[i+1 for i in ident['qmap_oxygens']]}")
    print(f"  O fenólico provável do QMAP: {None if ident['phenolic_o'] is None else ident['phenolic_o']+1}")

    if ident["phenolic_o"] is None or not ident["dmso_oxygens"]:
        raise RuntimeError("Não consegui identificar O(QMAP) ou O(DMSO). Verifique símbolos/estrutura.")

    charges = load_charges(args.charges, n_atoms, n_frames)
    can_dipole = charges is not None

    qmap_atoms = np.array(ident["qmap_atoms"], dtype=int)
    phenolic_o = ident["phenolic_o"]
    dmso_oxygens = ident["dmso_oxygens"]

    dihedral_idx = None
    if args.dihedral is not None:
        dihedral_idx = np.array(args.dihedral, dtype=int)
        if not args.zero_based:
            dihedral_idx -= 1
        if np.any(dihedral_idx < 0) or np.any(dihedral_idx >= n_atoms):
            raise ValueError("Índices do dihedral fora do intervalo.")

    rows = []
    dipoles = []
    min_oo = []
    dihedrals = []

    for fr, coords in enumerate(frames):
        rmin = min(distance(coords[phenolic_o], coords[o]) for o in dmso_oxygens)
        min_oo.append(rmin)

        if can_dipole:
            q = charges[fr, qmap_atoms]
            c = coords[qmap_atoms].copy()
            c -= np.mean(c, axis=0)
            mu, muvec = dipole_debye(c, q)
        else:
            mu, muvec = np.nan, np.array([np.nan, np.nan, np.nan])
        dipoles.append(mu)

        if dihedral_idx is not None:
            d = compute_dihedral(
                coords[dihedral_idx[0]],
                coords[dihedral_idx[1]],
                coords[dihedral_idx[2]],
                coords[dihedral_idx[3]],
            )
        else:
            d = np.nan
        dihedrals.append(d)

        rows.append([fr, mu, muvec[0], muvec[1], muvec[2], rmin, d])

    dipoles = np.array(dipoles, dtype=float)
    min_oo = np.array(min_oo, dtype=float)
    dihedrals = np.array(dihedrals, dtype=float)

    csv_path = os.path.join(args.outdir, "qmap_correlacoes.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["frame", "dipole_QMAP_D", "mu_x_D", "mu_y_D", "mu_z_D",
                    "min_OQMAP_ODMSO_A", "dihedral_C2_C1_C3_C4_deg"])
        w.writerows(rows)

    print(f"\nTabela salva em: {csv_path}")

    save_timeseries(
        np.arange(n_frames), min_oo,
        "min O(QMAP)–O(DMSO) distance (Å)",
        "Minimum O(QMAP)–O(DMSO) distance along the trajectory",
        os.path.join(args.outdir, "timeseries_min_OO.png")
    )

    if can_dipole:
        save_timeseries(
            np.arange(n_frames), dipoles,
            "QMAP dipole moment (D)",
            "QMAP dipole moment along the trajectory",
            os.path.join(args.outdir, "timeseries_dipole_qmap.png")
        )

        save_scatter(
            min_oo, dipoles,
            "min O(QMAP)–O(DMSO) distance (Å)",
            "QMAP dipole moment (D)",
            "Correlation between first-shell O–O contact and QMAP dipole",
            os.path.join(args.outdir, "correlation_OO_distance_vs_dipole.png"),
            window=args.smooth_window
        )

        if dihedral_idx is not None:
            save_timeseries(
                np.arange(n_frames), dihedrals,
                "Dihedral C2–C1–C3–C4 (deg)",
                "Backbone dihedral along the trajectory",
                os.path.join(args.outdir, "timeseries_dihedral.png")
            )

            save_scatter(
                dihedrals, dipoles,
                "Dihedral C2–C1–C3–C4 (deg)",
                "QMAP dipole moment (D)",
                "Correlation between backbone dihedral and QMAP dipole",
                os.path.join(args.outdir, "correlation_dihedral_vs_dipole.png"),
                window=args.smooth_window
            )
        else:
            print("\nAVISO: --dihedral não informado. A figura dihedral vs dipole não foi gerada.")
            print("Informe os quatro índices C2 C1 C3 C4 para gerar essa figura.")
    else:
        print("\nAVISO: sem arquivo de cargas. Figuras envolvendo dipolo não foram geradas.")

    report_path = os.path.join(args.outdir, "relatorio_execucao.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("Relatório de execução - QMAP CPMD correlações\n")
        f.write("=============================================\n\n")
        f.write(f"Trajetória: {args.traj}\n")
        f.write(f"Número de átomos: {n_atoms}\n")
        f.write(f"Número de frames: {n_frames}\n")
        f.write(f"DMSO encontrados: {len(ident['dmso_groups'])}\n")
        f.write(f"Átomos QMAP: {len(ident['qmap_atoms'])}\n")
        f.write(f"O(DMSO), índices 1-based: {[i+1 for i in ident['dmso_oxygens']]}\n")
        f.write(f"O(QMAP), índices 1-based: {[i+1 for i in ident['qmap_oxygens']]}\n")
        f.write(f"O fenólico provável, 1-based: {None if phenolic_o is None else phenolic_o+1}\n")
        f.write("\nArquivos principais gerados:\n")
        f.write("  qmap_correlacoes.csv\n")
        f.write("  timeseries_min_OO.png\n")
        if can_dipole:
            f.write("  timeseries_dipole_qmap.png\n")
            f.write("  correlation_OO_distance_vs_dipole.png\n")
            if dihedral_idx is not None:
                f.write("  timeseries_dihedral.png\n")
                f.write("  correlation_dihedral_vs_dipole.png\n")

    print(f"Relatório salvo em: {report_path}")
    print("\nConcluído.")


if __name__ == "__main__":
    main()
