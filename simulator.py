import numpy as np

from peptides import generate_peptides
from fragmentation import fragments
from dda import acquire as dda_acquire
from dia import acquire as dia_acquire
from scoring import score_against_library
from chromatography import assign_chromatography


def _fragment_vector(masses):
    vec = {}
    for m in masses:
        key = int(round(m))
        vec[key] = vec.get(key, 0) + 1
    return vec


def _ms2_abundance(abundance, masses):
    if not masses:
        return 0.0
    return float(abundance) / len(masses)


class Simulator:
    def run(self, n_peptides=1000, window=20, gradient_min=22, topn=10, scan_interval=0.6):
        peps = generate_peptides(n_peptides)
        peps_chrom = assign_chromatography(
            peps,
            rt_min=0.0,
            rt_max=gradient_min * 60.0,
            intensity=1.0
        )

        library = {}
        for p in peps:
            masses = fragments(p[1])
            library[p[0]] = {
                'sequence': p[1],
                'masses': masses,
                'vector': _fragment_vector(masses)
            }
        vectors = {pid: info['vector'] for pid, info in library.items()}

        time_axis = np.arange(0.0, gradient_min * 60.0 + scan_interval, scan_interval)

        dda_candidates = dda_acquire(peps_chrom, time_axis, topn=topn)
        selected_unique = []
        seen = set()
        for peptide in dda_candidates:
            if peptide[0] not in seen:
                seen.add(peptide[0])
                selected_unique.append(peptide)

        dda = []
        for peptide in selected_unique:
            scores = score_against_library(vectors[peptide[0]], vectors)
            masses = library[peptide[0]]['masses']
            n_fragments = len(masses)
            dda.append({
                'peptide': peptide,
                'scores': scores,
                'n_fragments': n_fragments,
                'ms2_abundance': _ms2_abundance(peptide[5], masses)
            })

        dia_bins = dia_acquire(peps_chrom, window)
        dia = {}
        for low, peptides in dia_bins.items():
            window_spectrum = {}
            for peptide in peptides:
                for key, value in vectors[peptide[0]].items():
                    window_spectrum[key] = window_spectrum.get(key, 0) + value
            window_scores = score_against_library(window_spectrum, vectors)
            window_items = []
            for peptide in peptides:
                masses = library[peptide[0]]['masses']
                n_fragments = len(masses)
                window_items.append({
                    'peptide': peptide,
                    'scores': window_scores,
                    'n_fragments': n_fragments,
                    'ms2_abundance': round(_ms2_abundance(peptide[5], masses), 4)
                })
            dia[low] = window_items

        chromatogram = {
            'peptides': peps_chrom,
            'time_axis': time_axis.tolist()
        }

        return {
            'peptides': peps,
            'chromatogram': chromatogram,
            'ms1': chromatogram,
            'dda': dda,
            'dia': dia,
            'ms2': {
                'dda': dda,
                'dia': dia
            },
            'library': library
        }
