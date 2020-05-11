#!/usr/bin/env python

from Logic import OrthoParser, MSA_runner, AlignmentParser
from Custom_Objects import NoGenesInOrthogroupError
from os import path


######################

orthofile = "/home/joris/nemaNAS/steen176/Ava/data/orthofinder_output_ORFs_75_110/Results_Mar04/Orthogroups/Orthogroups.tsv"
alignment_outfolder = "/home/joris/nemaNAS/steen176/Ava/data/Incisive_out_75_110"

genome = 'meloidogyne_incognita.PRJEB8714.WBPS14.genomic_softmasked.fa.orfs_75_110'
genome_protein_fasta = "/home/joris/nemaNAS/steen176/Ava/data/genome_sequences/orfs/orfs_75_110/meloidogyne_incognita.PRJEB8714.WBPS14.genomic_softmasked.fa.orfs_75_110.fa"
genome_CDS_fasta = "/home/joris/nemaNAS/steen176/Ava/data/genome_sequences/cds/incognita_orf_75_110_cds.fa"




def perform_MSA(orthogroup, fasta) -> str:
    alignment_outfile = f"{alignment_outfolder}/{orthogroup.group_name}_{path.basename(fasta)}_alignment.fasta"
    msa = MSA_runner(fasta, memory_in="/dev/shm/input.fa", alignment_out=alignment_outfile)
    try:
        msa.import_orthogroup_sequences(orthogroup)
        alignment_file = msa.align_orthogroup_sequences(tool="kalign")

        return alignment_file
    except NoGenesInOrthogroupError:
        return False


def calc_score(alignment_file):
    alignment_parser = AlignmentParser(alignment_file)
    alignment_parser.import_alignment()
    comp_z_score = alignment_parser.get_positional_scores()

    return comp_z_score


def query_single_group(query):
    ortho_parser = OrthoParser(orthofile, genome)
    groups = [group for group in ortho_parser.get_orthogroups() if group.group_name == query]
    for group in groups:
        alignment_file = perform_MSA(group)
        if not alignment_file:
            continue
        try:
            comp_z_score = calc_score(alignment_file)

            print(f"{group.group_name}\t{comp_z_score}")
        except ValueError:
            pass
        except ZeroDivisionError:
            pass


if __name__ == "__main__":
    header = "Orthogroup\tProteinZscore\tNuclZscore"
    ortho_parser = OrthoParser(orthofile, genome)
    for group in ortho_parser.get_orthogroups():
        prot_alignment_file = perform_MSA(group, genome_protein_fasta)
        nucl_alignment_file = perform_MSA(group, genome_CDS_fasta)

        if not prot_alignment_file or not nucl_alignment_file:
            continue

        try:
            prot_comp_z_score = calc_score(prot_alignment_file)
            nucl_comp_z_score = calc_score(nucl_alignment_file)

            print(f"{group.group_name}\t{prot_comp_z_score}\t{nucl_comp_z_score}")
        except ValueError:
            pass
        except ZeroDivisionError:
            pass


