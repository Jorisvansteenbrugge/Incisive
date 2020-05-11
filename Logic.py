from Custom_Objects import Orthogroup, Sequence, NoGenesInOrthogroupError
from statistics import mean, stdev
from Bio import SeqIO
from subprocess import call
from collections import Counter


def get_genes(col, strip_names=True):
    genes = col.split(',')
    if '=' in genes[0] and strip_names:
        genes = [gene.split('=')[-1] for gene in genes]
    return genes


class OrthoParser:

    def __init__(self, orthogroups_file, genome):
        self.orthofile = orthogroups_file
        self.genome = genome
        self.orthogroups = []

        self.parse_orthogroups()

    def parse_orthogroups(self, sep='\t'):
        with open(self.orthofile) as orthofile:
            header = orthofile.readline().strip().split(sep)
            genome_col = header.index(self.genome)

            for line in orthofile:
                line = line.strip().split(sep)
                group_name = line[0]

                try:
                    genes = get_genes(line[genome_col], strip_names=False)
                    self.orthogroups.append(Orthogroup(group_name, genes))
                except IndexError:  # genome has no genes in the orthogroup
                    pass

    def get_orthogroups(self):
        return self.orthogroups


class MSA_runner:

    def __init__(self, genome_protein_fasta, memory_in="/dev/shm/alignment_input.fa",
                 alignment_out="/dev/shm/alignment_MSA.fa"):

        self.memory_in = memory_in
        self.alignment_out = alignment_out
        self.protein_dict = self._make_fasta_dict(genome_protein_fasta)
        self.protein_sequences = []

    @staticmethod
    def _make_fasta_dict(protein_fasta, strip_ids = True):
        if strip_ids:
            return {record.id.replace(" ", ""): str(record.seq) for record in SeqIO.parse(protein_fasta, 'fasta')}
        else:
            return {record.id: str(record.seq) for record in SeqIO.parse(protein_fasta, 'fasta')}

    def import_orthogroup_sequences(self, orthogroup):
        for gene in orthogroup.genes:
            try:
                seq = self.protein_dict[gene]
            except KeyError:
                raise NoGenesInOrthogroupError

            fasta_entry = f">{gene}\n{seq}"
            self.protein_sequences.append(fasta_entry)

    def _sequences_to_file(self):
        with open(self.memory_in, 'w') as outfile:
            outfile.write('\n'.join(self.protein_sequences) + "\n")

    def align_orthogroup_sequences(self, tool='muscle'):
        self._sequences_to_file()

        cmd = f"{tool} < {self.memory_in} > {self.alignment_out} 2> /dev/null"
        # print(f"Running: {cmd}")

        call(cmd, shell=True)  # passing to speedup testing

        return self.alignment_out


class AlignmentParser:

    def __init__(self, alignment_file):
        self.alignment_file = alignment_file
        self.sequences = []

    def import_alignment(self):
        for record in SeqIO.parse(self.alignment_file, 'fasta'):
            self.sequences.append(Sequence(record.id.replace(" ", ""), str(record.seq)))

    def _get_longest_seq_len(self):
        return max([len(seq) for seq in self.sequences])

    def get_nuc_at_pos(self, pos):
        return [seq.get_nuc(pos) for seq in self.sequences]

    def to_skip(self, pos, num_sequences, seq_cutoff_fraction):
        gap_c = self.get_nuc_at_pos(pos).count('-')
        if gap_c > (num_sequences * seq_cutoff_fraction):
            return True  # Skip
        else:
            return False

    def get_majority_allel(self, pos):
        alleles = self.get_nuc_at_pos(pos)
        counts = dict(Counter(alleles))

        majority_allele = max(counts)
        return (majority_allele, counts[majority_allele])

    @staticmethod
    def _to_z_score(counts, mu, sigma):
        z_scores = [(c - mu) / sigma for c in counts]
        return z_scores

    def get_positional_scores(self, cleanup=True, seq_cutoff_fraction=0.3):
        num_sequences = len(self.sequences)
        max_length = self._get_longest_seq_len()

        majority_alleles = []
        majority_counts = []
        majority_scores = []

        for i in range(max_length):
            if cleanup and self.to_skip(i, num_sequences, seq_cutoff_fraction):
                continue

            majority_allele, count = self.get_majority_allel(i)
            positional_score = (count / num_sequences)

            majority_alleles.append(majority_allele)
            majority_counts.append(count)
            majority_scores.append(positional_score)

        mean_score = mean(majority_counts)
        stdev_score = stdev(majority_counts)

        z_scores = self._to_z_score(majority_counts, mean_score, stdev_score)
        composite_z_score = sum(z_scores)

        return composite_z_score
