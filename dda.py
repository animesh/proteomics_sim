from chromatography import peak


def acquire(peptides, time_axis, topn=10):
    selections = []

    for t in time_axis:
        active = []

        for p in peptides:
            intensity = peak(p[4], p[5], t)
            if intensity > 0:
                active.append((p, intensity))

        active.sort(key=lambda x: x[1], reverse=True)
        selections.extend([p for p, _ in active[:topn]])

    return selections