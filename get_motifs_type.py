def motifs_GC(n=4):
    base_list = ['A', 'T', 'G', 'C']
    motif_set_n = dict((k, []) for k in range(1, n + 1))
    motif_set_n[1] = base_list
    for k in range(2, n + 1):
        motif_set_n[k] = [i + j for i in base_list for j in motif_set_n[k - 1]]
    return motif_set_n
