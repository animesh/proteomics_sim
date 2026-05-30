import os

import matplotlib.pyplot as plt
import numpy as np
from simulator import ProteomicsStudy
from digest import generate_peptides
from chromatography import peak
from fragmentation import fragment_efficiency

n_peptides = 1000
n_replicates = 10
gradient_min = 22
dia_window = 2.0
seed = 42

output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

summary_filename = os.path.join(
    output_dir,
    f"summary_n{n_peptides}_rep{n_replicates}_g{gradient_min}_w{dia_window:.1f}.png"
)
stacks_filename = os.path.join(
    output_dir,
    f"chromatogram_stacks_n{n_peptides}_rep{n_replicates}_g{gradient_min}_w{dia_window:.1f}.png"
)

study = ProteomicsStudy(seed=seed)
results = study.run(
    n_peptides=n_peptides,
    n_replicates=n_replicates,
    gradient_min=gradient_min,
    dia_window=dia_window
)

ms1 = results["ms1"]
dda = results["dda"]
dia = results["dia"]
cv = results["cv"]
cv_dia = results.get("cv_dia")

print("Simulation complete")
print(f"MS1 peptides : {cv.shape[0]}")
print(f"DDA selections : {dda.shape[0]}")
print(f"DIA windows : {dia[["low", "high"]].drop_duplicates().shape[0]}")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes[0, 0].hist(cv["cv"].dropna(), bins=40, color="tab:blue", alpha=0.8)
axes[0, 0].set_title("MS1 CV distribution")
axes[0, 0].set_xlabel("CV (%)")
axes[0, 0].set_ylabel("Peptides")

dda_counts = dda.groupby("peptide_id").size()
axes[0, 1].hist(dda_counts, bins=30, color="tab:orange", alpha=0.8)
axes[0, 1].set_title("DDA selections per peptide")
axes[0, 1].set_xlabel("Selections")
axes[0, 1].set_ylabel("Peptides")

dia_mean = dia.groupby(["low", "high"])["signal"].mean().reset_index()
axes[1, 0].plot(dia_mean["low"], dia_mean["signal"], marker="o", linestyle="-", markersize=3)
axes[1, 0].set_title("Mean DIA window signal")
axes[1, 0].set_xlabel("Window low m/z")
axes[1, 0].set_ylabel("Mean signal")

if cv_dia is not None and not cv_dia["cv"].dropna().empty:
    axes[1, 1].hist(cv_dia["cv"].dropna(), bins=40, color="tab:green", alpha=0.8)
    axes[1, 1].set_title("DIA window CV distribution")
    axes[1, 1].set_xlabel("CV (%)")
    axes[1, 1].set_ylabel("Windows")
else:
    axes[1, 1].text(0.5, 0.5, "DIA CV not available", ha="center", va="center")
    axes[1, 1].set_xticks([])
    axes[1, 1].set_yticks([])

fig.tight_layout()
fig.savefig(summary_filename, dpi=150)
print(f"Saved summary plot to {summary_filename}")
plt.close(fig)

# Build stacked replicate chromatograms for DDA and DIA
gradient_length_sec = gradient_min * 60
pep = generate_peptides(
    n=n_peptides,
    gradient_length_sec=gradient_length_sec,
    seed=seed
)

time_axis = np.arange(0, gradient_length_sec, 0.6)
rng = np.random.default_rng(seed)

all_dda_traces = []
all_dia_traces = []
dda_picks = []
for rep in range(n_replicates):
    pep_rep = pep.copy()
    pep_rep["rt"] += rng.normal(0, 0.4, len(pep_rep))
    pep_rep["abundance"] *= rng.lognormal(0, 0.15, len(pep_rep))

    dda_trace = []
    pick_times = []
    exclusion = {}
    for t in time_axis:
        active = []
        for r in pep_rep.itertuples():
            inten = peak(r.rt, r.abundance, t)
            if inten > 100:
                active.append((r, inten))
        active.sort(key=lambda x: x[1], reverse=True)

        selected_intensity = 0.0
        selected_any = False
        selected_count = 0
        for r, inten in active:
            if r.peptide_id in exclusion and t < exclusion[r.peptide_id]:
                continue
            exclusion[r.peptide_id] = t + 5
            selected_intensity += inten
            selected_any = True
            selected_count += 1
            if selected_count >= 10:
                break
        dda_trace.append(selected_intensity)
        if selected_any:
            pick_times.append(t)

    dia_trace = []
    for t in time_axis:
        total_signal = 0.0
        for r in pep_rep.itertuples():
            total_signal += peak(r.rt, r.abundance, t) * fragment_efficiency(30, r.optimal_ce, 5)
        dia_trace.append(total_signal)

    all_dda_traces.append(np.array(dda_trace))
    all_dia_traces.append(np.array(dia_trace))
    dda_picks.append(np.array(pick_times))

stack_fig, (ax_dda, ax_dia) = plt.subplots(1, 2, figsize=(18, 8), sharex=False, sharey=False)

dda_offset = max(trace.max() for trace in all_dda_traces) * 1.2
for idx, trace in enumerate(all_dda_traces):
    ax_dda.plot(time_axis / 60, trace + idx * dda_offset, color="tab:blue", alpha=0.7)
    ax_dda.scatter(
        dda_picks[idx] / 60,
        np.full_like(dda_picks[idx], idx * dda_offset),
        s=12,
        color="black",
        marker="|",
        label="DDA pick" if idx == 0 else ""
    )
ax_dda.set_title("DDA replicate chromatograms with pick markers")
ax_dda.set_xlabel("Retention time (min)")
ax_dda.set_ylabel("Replicate number")
ax_dda.set_yticks([idx * dda_offset for idx in range(n_replicates)])
ax_dda.set_yticklabels([str(idx + 1) for idx in range(n_replicates)])
ax_dda.legend(loc="upper right")
ax_dda.text(0.02, 0.95, "DDA", transform=ax_dda.transAxes, fontsize=12, fontweight="bold", va="top")

bins = np.arange(400, 1300 + dia_window, dia_window)
window_counts = [np.sum((pep["mz"] >= low) & (pep["mz"] < low + dia_window)) for low in bins[:-1]]
dia_matrix = np.tile(window_counts, (n_replicates, 1))

im = ax_dia.imshow(
    dia_matrix,
    aspect="auto",
    cmap="YlGn",
    origin="lower",
    extent=[400, 1300, 1, n_replicates],
    interpolation="nearest"
)
ax_dia.set_title("DIA window counts per replicate")
ax_dia.set_xlabel("Window low m/z")
ax_dia.set_ylabel("Replicate number")
ax_dia.set_yticks(np.arange(1, n_replicates + 1))
ax_dia.set_xticks(np.linspace(400, 1300, 10))
ax_dia.text(0.02, 0.95, "DIA", transform=ax_dia.transAxes, fontsize=12, fontweight="bold", va="top")
stack_fig.colorbar(im, ax=ax_dia, label="Number of precursors")

stack_fig.tight_layout()
stack_fig.savefig(stacks_filename, dpi=150)
print(f"Saved stacked chromatogram plot to {stacks_filename}")
plt.close(stack_fig)
