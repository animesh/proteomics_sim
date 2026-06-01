
import numpy as np
from peptides import generate_peptides
from chrom_assign import assign_chromatography
from fragmentation import fragments
from dda import acquire as dda_acquire
from dia import acquire as dia_acquire
from scoring import score_against_library


def _fragment_vector(masses):
    vector = {}
    for m in masses:
        key = int(round(m))
        vector[key] = vector.get(key, 0) + 1
    return vector


def _ms2_abundance(abundance, masses):
    return float(abundance) / len(masses) if masses else 0.0


class Simulator:
    def run(self, n_peptides=1000, window=20, gradient_min=22, topn=10, scan_interval=0.6):
        peps = generate_peptides(n_peptides)
        peps = assign_chromatography(peps, gradient_min * 60.0)

        time_axis = np.arange(0, gradient_min * 60.0 + scan_interval, scan_interval)

        library = {}
        vectors = {}
        for p in peps:
            masses = fragments(p[1])
            library[p[0]] = {
                'sequence': p[1],
                'masses': masses,
                'vector': _fragment_vector(masses)
            }
            vectors[p[0]] = library[p[0]]['vector']

        dda_selections = dda_acquire(peps, time_axis, topn=topn)
        selected_by_id = {}
        for t, peptide, intensity in dda_selections:
            if peptide[0] not in selected_by_id:
                selected_by_id[peptide[0]] = peptide

        dda_entries = []
        for peptide in selected_by_id.values():
            scores = score_against_library(vectors[peptide[0]], vectors)
            masses = library[peptide[0]]['masses']
            dda_entries.append({
                'peptide': peptide,
                'scores': scores,
                'n_fragments': len(masses),
                'ms2_abundance': _ms2_abundance(peptide[5], masses)
            })

        dia_windows = dia_acquire(peps, time_axis, window=window)
        dia_entries = {}
        for low, observations in dia_windows.items():
            window_vector = {}
            for obs in observations:
                peptide = obs['peptide']
                for key, value in vectors[peptide[0]].items():
                    window_vector[key] = window_vector.get(key, 0) + value
            window_scores = score_against_library(window_vector, vectors)
            dia_entries[low] = []
            for obs in observations:
                peptide = obs['peptide']
                dia_entries[low].append({
                    'peptide': peptide,
                    'scores': window_scores,
                    'n_fragments': len(library[peptide[0]]['masses']),
                    'trace_points': len(obs['trace']),
                    'chrom_area': obs['area'],
                    'ms2_abundance': _ms2_abundance(
                        peptide[5],
                        library[peptide[0]]['masses']
                    )
                })

        chromatogram = {
            'peptides': peps,
            'time_axis': time_axis
        }

        return {
            'library': library,
            'peptides': peps,
            'chromatogram': chromatogram,
            'ms1': chromatogram,
            'ms2': {
                'dda': dda_entries,
                'dia': dia_entries
            },
            'dda': dda_selections,
            'dia': dia_windows,
            'metrics': {
                'dda_unique_precursors': len(selected_by_id),
                'dia_unique_precursors': len({obs['peptide'][0] for w in dia_windows.values() for obs in w}),
                'dia_windows': len(dia_windows)
            }
        }
