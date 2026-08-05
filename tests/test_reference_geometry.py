import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/python'))
from qmap_cpmd_common import first_xyz_frame,load_mapping,validate_symbols,mic_distance

def test_reference_geometry():
    m=load_mapping(ROOT/'config/atom_indices.json')
    f=first_xyz_frame(ROOT/'structures/GEOMETRY_reference.xyz')
    validate_symbols(f.symbols,m)
    assert len(m['qmap']['atom_indices'])==31
    assert len(m['dmso_molecules'])==21
    assert 0.8 < mic_distance(f.coordinates[198],f.coordinates[11],16.0) < 1.2
    assert m['qmap']['labels']['N1']==198
