"""Fit Cartesian coordinates to boron planes and the original bond lengths."""
import numpy as np
from rdkit import Chem


def _norm(v):
    return np.linalg.norm(v, axis=1)


def _unit(v):
    return v / np.maximum(_norm(v)[:, None], 1e-12)


def _phase(xyz, torsions):
    a, b, c, d = xyz[torsions].transpose(1, 0, 2)
    axis = _unit(c - b)
    u, v = a - b, d - c
    u = _unit(u - np.einsum('ij,ij->i', u, axis)[:, None] * axis)
    v = _unit(v - np.einsum('ij,ij->i', v, axis)[:, None] * axis)
    return np.column_stack((np.einsum('ij,ij->i', u, v),
                            np.einsum('ij,ij->i', np.cross(axis, u), v)))


def _heights(xyz, centers):
    b, a, c, d = xyz[centers].transpose(1, 0, 2)
    normal = np.cross(c - a, d - a)
    if np.any(_norm(normal) < 1e-8):
        raise ValueError('Boron neighbors do not define a plane.')
    return np.einsum('ij,ij->i', b - a, _unit(normal))


def project_boron_geometry(mol, start, torsions, bond_tolerance):
    from scipy.optimize import least_squares
    from scipy.sparse import csr_matrix

    n = len(start)
    neighbors = [[a.GetIdx() for a in atom.GetNeighbors()] for atom in mol.GetAtoms()]
    bonds = np.array([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])
    centers = np.array([(i, *adj) for i, adj in enumerate(neighbors)
                        if mol.GetAtomWithIdx(i).GetAtomicNum() == 5 and len(adj) == 3])
    torsion_atoms = np.array([t[:4] for t in torsions], dtype=int).reshape(-1, 4)
    target = np.radians([t[4] for t in torsions])
    phase = np.column_stack((np.cos(target), np.sin(target)))
    angles = np.array([(a, i, b) for i, adj in enumerate(neighbors)
                       for j, a in enumerate(adj) for b in adj[j+1:]])

    def distances(xyz, pairs):
        return _norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]])

    def angle_cosines(xyz):
        a, b, c = xyz[angles].transpose(1, 0, 2)
        return np.einsum('ij,ij->i', _unit(a - b), _unit(c - b))

    lengths = distances(start, bonds)
    initial_angles = angle_cosines(start)
    # Exclude bonded and 1–3 pairs; do not require curing existing short contacts.
    excluded = np.eye(n, dtype=bool)
    for i, adj in enumerate(neighbors):
        excluded[i, adj] = True
        for j in adj:
            excluded[i, neighbors[j]] = True
    pairs = np.column_stack(np.where(np.triu(~excluded, 1)))
    table = Chem.GetPeriodicTable()
    radii = np.array([table.GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    initial_contacts = distances(start, pairs)
    floors = np.minimum(0.72 * radii[pairs].sum(axis=1), initial_contacts)
    keep = initial_contacts < 4.0
    coords = start.copy()

    for _ in range(3):
        contacts = pairs[keep]
        blocks = (bonds, np.repeat(torsion_atoms, 2, axis=0), centers, angles, contacts)
        row_indices, col_indices = [np.arange(n * 3)], [np.arange(n * 3)]
        offset = n * 3
        for atoms in blocks:
            cols = (3 * atoms[..., None] + np.arange(3)).reshape(len(atoms), atoms.shape[1] * 3)
            row_indices.append(np.repeat(np.arange(offset, offset + len(atoms)), cols.shape[1]))
            col_indices.append(cols.ravel())
            offset += len(atoms)
        rows, cols = np.concatenate(row_indices), np.concatenate(col_indices)
        sparsity = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(offset, n * 3))

        def residual(flat):
            xyz = flat.reshape(n, 3)
            return np.concatenate(((0.15 * (xyz - start)).ravel(),
                                   150 * (distances(xyz, bonds) - lengths),
                                   120 * (_phase(xyz, torsion_atoms) - phase).ravel(),
                                   150 * _heights(xyz, centers),
                                   1.5 * (angle_cosines(xyz) - initial_angles),
                                   30 * np.maximum(floors[keep] - distances(xyz, contacts), 0)))

        fit = least_squares(residual, coords.ravel(), jac_sparsity=sparsity, max_nfev=150,
                            ftol=1e-6, gtol=1e-6, tr_options={'atol': 1e-8, 'btol': 1e-8})
        coords = fit.x.reshape(n, 3)
        new = (distances(coords, pairs) < floors - 0.03) & ~keep
        if not np.any(new):
            break
        keep |= new

    drift = np.abs(distances(coords, bonds) - lengths).max()
    measured = _phase(coords, torsion_atoms)
    delta = np.degrees(np.arctan2(measured[:, 1], measured[:, 0]) - target)
    error = np.abs((delta + 180) % 360 - 180).max(initial=0)
    height = np.abs(_heights(coords, centers)).max()
    clash = (floors - distances(coords, pairs)).max(initial=0)
    if (not np.isfinite(coords).all() or not np.allclose(_norm(measured), 1, atol=1e-6) or
            drift > bond_tolerance or error > 0.5 or height > 0.001 or clash > 0.03):
        raise ValueError(f'Bond drift {drift:.4f} A; torsion error {error:.3f} deg; '
                         f'B plane deviation {height:.4f} A; contact deficit {clash:.4f} A.')
    return coords, drift, error, height
