class Orthogroup:
    def __init__(self, group_name, genes):
        self.group_name = group_name
        genes_spl = [g[1:].split(' ')[0] if g.startswith(' ') else g.split(' ')[0] for g in genes ]
        

        self.genes = genes_spl

    def __repr__(self):
        return f"<obj: {self.group_name}>"

    def __str__(self):
        return f"{self.group_name}: {(','.join(self.genes))}"


class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def get_nuc(self, pos):
        return self.seq[pos]


class NoGenesInOrthogroupError(Exception):
    pass