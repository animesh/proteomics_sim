import argparse
from pathlib import Path

from simulator import Simulator


def plot_chromatogram(r, output_path, label_count=20):
    import matplotlib.pyplot as plt
    import numpy as np
    from chromatography import peak

    time_axis = np.array(r['chromatogram']['time_axis'])
    traces = []
    for p in r['ms1']['peptides']:
        traces.append((p, peak(p[4], p[5], time_axis)))

    total_intensity = np.zeros_like(time_axis)
    for _, trace in traces:
        total_intensity += trace

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(time_axis, total_intensity, color='black', lw=1.5, label='Total MS1 intensity')

    peptide_points = sorted(traces, key=lambda item: item[0][5], reverse=True)[:label_count]
    for p, trace in peptide_points:
        rt = p[4]
        intensity = float(peak(rt, p[5], np.array([rt]))[0])
        ax.scatter([rt], [intensity], color='red', s=30, zorder=3)
        ax.annotate(
            p[1],
            xy=(rt, intensity),
            xytext=(0, 8),
            textcoords='offset points',
            ha='center',
            va='bottom',
            fontsize=8,
            rotation=45,
            color='blue',
        )

    ax.set_title(f'Chromatogram with peptide labels (n={len(r["ms1"]["peptides"])})')
    ax.set_xlabel('Retention time (s)')
    ax.set_ylabel('MS1 intensity proxy')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    print(f'Saved chromatogram plot to {output_path}')


parser = argparse.ArgumentParser(description='Run a toy DDA/DIA peptide simulation.')
parser.add_argument('--n_peptides', type=int, default=3, help='Number of peptides to generate')
parser.add_argument('--gradient_min', type=float, default=1.0, help='Gradient length in minutes')
parser.add_argument('--window', type=float, default=2.0, help='DIA window width in Th')
parser.add_argument('--topn', type=int, default=1, help='Top N precursors to select for DDA')
parser.add_argument('--plot', action='store_true', help='Generate MS1 chromatogram plot with peptide labels')
parser.add_argument('--plot_file', type=str, default=None, help='Output file path for the chromatogram plot')
args = parser.parse_args()

if args.plot and args.plot_file is None:
    sanitized_gradient = str(args.gradient_min).replace('.', 'p')
    sanitized_window = str(args.window).replace('.', 'p')
    args.plot_file = (
        f'plots/example_n{args.n_peptides}_g{sanitized_gradient}_w{sanitized_window}'
        f'_topn{args.topn}_chromatogram.png'
    )

r = Simulator().run(
    n_peptides=args.n_peptides,
    window=args.window,
    gradient_min=args.gradient_min,
    topn=args.topn
)

library = r['library']

print('Generated peptides:')
for p in r['peptides']:
    print(p)

ms1_peptides = r['ms1']['peptides']
ms1_total_abundance = sum(p[5] for p in ms1_peptides)
ms1_mean_abundance = ms1_total_abundance / len(ms1_peptides) if ms1_peptides else 0.0
print('\nMS1 summary:')
print(f'  precursor count: {len(ms1_peptides)}')
print(f'  total abundance: {ms1_total_abundance:.4f}')
print(f'  mean abundance: {ms1_mean_abundance:.4f}')
print('\nMS1 peptides:')
for p in ms1_peptides[: min(10, len(ms1_peptides))]:
    print(f'  id={p[0]} seq={p[1]} mz={p[3]:.4f} rt={p[4]:.4f} abundance={p[5]:.4f}')

print('\nDDA selected peptides with MS2 proxy and scoring:')
for item in r['ms2']['dda']:
    peptide = item['peptide']
    best_match_id, best_match_score = max(item['scores'].items(), key=lambda kv: kv[1])
    best_match_seq = library[best_match_id]['sequence']
    print(
        f'  id={peptide[0]} seq={peptide[1]} mz={peptide[3]:.4f} rt={peptide[4]:.4f} '
        f'abundance={peptide[5]:.4f} n_fragments={item["n_fragments"]} '
        f'ms2_abundance={item["ms2_abundance"]:.4f} top_match={best_match_seq} ({best_match_score:.4f})'
    )
    rounded_scores = {k: round(v, 4) for k, v in item['scores'].items()}
    print(f'    scores={rounded_scores}')
msa_total = sum(item['peptide'][5] for item in r['ms2']['dda'])
ms2_total = sum(item['ms2_abundance'] * item['n_fragments'] for item in r['ms2']['dda'])
print(f'\nDDA MS2 summary: {len(r["ms2"]["dda"])} DDA spectra selected')
print(f'  total MS1 abundance selected: {msa_total:.4f}')
print(f'  total MS2 abundance proxy: {ms2_total:.4f}')
print(f'  mean MS2 abundance per fragment: {sum(item["ms2_abundance"] for item in r["ms2"]["dda"])/len(r["ms2"]["dda"]):.4f}')

sorted_ms1 = sorted(ms1_peptides, key=lambda p: p[4])
coeluting_pair = None
for p1, p2 in zip(sorted_ms1, sorted_ms1[1:]):
    dt = p2[4] - p1[4]
    if dt <= 8.0:
        coeluting_pair = (dt, p1, p2)
        break

if coeluting_pair:
    dt, p1, p2 = coeluting_pair
    selected_ids = {item['peptide'][0] for item in r['ms2']['dda']}
    print('\nDDA coelution check:')
    print(f'  closest RT separation: {dt:.4f} seconds')
    print(f'  peptide1 selected: {p1[0] in selected_ids} id={p1[0]} seq={p1[1]} rt={p1[4]:.4f} abundance={p1[5]:.4f}')
    print(f'  peptide2 selected: {p2[0] in selected_ids} id={p2[0]} seq={p2[1]} rt={p2[4]:.4f} abundance={p2[5]:.4f}')

print('\nDIA windows with MS2 proxy per window and scoring:')
ms2_dia_total = 0.0
for low, peptides in sorted(r['ms2']['dia'].items()):
    print(f'  window {low} - {low + args.window}: {len(peptides)} peptides')
    if peptides:
        window_scores = peptides[0]['scores']
        top_n = min(len(peptides), 2)
        top_matches = sorted(window_scores.items(), key=lambda kv: kv[1], reverse=True)[:top_n]
        top_matches_str = ', '.join(
            f"{library[pid]['sequence']} ({score:.4f})" for pid, score in top_matches
        )
        print(f'    window top_matches: {top_matches_str}')
    for item in peptides:
        peptide = item['peptide']
        ms2_dia_total += item['ms2_abundance'] * item['n_fragments']
        print(
            f'    id={peptide[0]} seq={peptide[1]} mz={peptide[3]:.4f} rt={peptide[4]:.4f} '
            f'abundance={peptide[5]:.4f} n_fragments={item["n_fragments"]} '
            f'ms2_abundance={item["ms2_abundance"]:.4f}'
        )
        rounded_scores = {k: round(v, 4) for k, v in item['scores'].items()}
        print(f'      scores={rounded_scores}')
print(f'\nDIA MS2 summary: {len(r["ms2"]["dia"])} DIA windows generated')
print(f'  total MS2 abundance proxy across DIA peptides: {ms2_dia_total:.4f}')

if args.plot:
    plot_chromatogram(r, args.plot_file)
