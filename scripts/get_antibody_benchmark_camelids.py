import sys
import os
from Bio import SeqIO, BlastIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import os

def is_sequence_in_matches(query_file: str, matches: List[str], seqres_path: str, perc_id_cutoff: float, cover_cutoff: float) -> bool:
    """
    Check if a query sequence is present in a list of matches.
    
    Args:
        query_file (str): Path to the query fasta file.
        matches (List[str]): List of match IDs.
        seqres_path (str): Path to the seqres fasta file.
        perc_id_cutoff (float): Minimum percentage identity required for a match.
        cover_cutoff (float): Minimum coverage percentage required for a match.
    
    Returns:
        bool: True if the query sequence is present in the matches, False otherwise.
    """
    with open(query_file) as query:
        query_seq = SeqIO.read(query, 'fasta')
        for match in matches:
            match_file = match + '.fasta'
            with open(match_file) as fasta:
                match_seqs = SeqIO.parse(fasta, 'fasta')
                for match_seq in match_seqs:
                    alignments = pairwise2.align.globalxx(query_seq.seq, match_seq.seq)
                    identity = float(pairwise2.format_alignment(*alignments[0]).split()[3].split('/')[0]) / float(pairwise2.format_alignment(*alignments[0]).split()[3].split('/')[1])
                    cov_perc = float(pairwise2.format_alignment(*alignments[0]).split()[3].split('/')[1]) / len(query_seq.seq)
                    if identity < perc_id_cutoff and cov_perc < cover_cutoff:
                        return False
    return True



def match_sequences(data_dir: str, blast_path: str, test_pdb: str='', perc_id_cutoff: float=0.98, cover_cutoff: float=0.8, res_cutoff: float=3.25) -> None:
    """
    Find matching sequences in the given data directory and BLAST database.
    
    Args:
        data_dir (str): Path to the data directory.
        blast_path (str): Path to the BLAST executable.
        test_pdb (str, optional): Test PDB ID to limit the search to. Defaults to ''.
        perc_id_cutoff (float, optional): Minimum percentage identity required for a match. Defaults to 0.98.
        cover_cutoff (float, optional): Minimum coverage percentage required for a match. Defaults to 0.8.
        res_cutoff (float, optional): Maximum resolution cutoff for PDB structures. Defaults to 3.25.
    """

    sabdab_path = os.path.join(data_dir, "all_antibody_structures_8_2020.txt")
    seqres_path = os.path.join(data_dir, "pdb_seqres_8_2020.fa")
    pdb_res_path = os.path.join(data_dir, "pdb_resolutions_8_2020.txt")

    with open(pdb_res_path) as file:
        res_lines = file.read().splitlines()

    pdbres = {}
    no_res_yet = True
    for line in res_lines:
        if line.startswith("IDCODE"):
            no_res_yet = False
        elif no_res_yet:
            continue
        else:
            pdb, *_, res = line.split()
            pdbres[pdb.lower()] = float(res)

    with open(sabdab_path) as file:
        sabdab_lines = file.read().splitlines()

    pdb_processed = {}
    for line in sabdab_lines:
        fields = line.split("\t")
        pdb = fields[0]
        if pdb_processed.get(pdb, None):
            continue
        
        pdb_processed[pdb] = 1
        hchain, lchain, model, ag_chain, ag_type, *_, res, method = fields

        if model == '1' or pdbres.get(pdb, 0) <= 0 or pdbres[pdb] > res_cutoff:
            continue
        if hchain == 'NA' or lchain != 'NA' or ag_chain == 'NA':
            continue
        if not (ag_type == 'protein' or ag_type.startswith('protein')):
            continue
        if test_pdb and pdb.lower() != test_pdb.lower():
            continue

        ag_chains = ag_chain.split('|')
        ag_chains = [item.strip() for item in ag_chains]

        query_seq_records = list(SeqIO.parse(seqres_path, 'fasta'))
        chain_sequences = get_sequences(pdb, [hchain] + ag_chains, query_seq_records)

        chain_matches = match_chain_sequences(chain_sequences, query_seq_records, blast_path, perc_id_cutoff, cover_cutoff)

        for chain_match in chain_matches:
            if chain_match == '':
                continue

        unbound_ab_matches = get_unbound_matches(chain_matches[0], chain_matches, seqres_path, pdbres, perc_id_cutoff, cover_cutoff)
        unbound_ag_matches = get_unbound_matches(chain_matches[1:], chain_matches, seqres_path, pdbres, perc_id_cutoff, cover_cutoff)
        print(f'{pdb}\t{unbound_ab_matches}\t{unbound_ag_matches}')


def get_sequences(pdb_id: str, chains: List[str], seq_records: List[SeqRecord]) -> List[str]:
    """
    Get the sequences for the specified PDB chains.
    
    Args:
        pdb_id (str): PDB ID.
        chains (List[str]): List of chain IDs.
        seq_records (List[SeqRecord]): List of SeqRecords.
    
    Returns:
        List[str]: List of sequences corresponding to the specified chains.
    """
    sequences = []
    for chain in chains:
        for record in seq_records:
            if f'{pdb_id}_{chain}' in record.id:
                sequences.append(str(record.seq))
                break
    return sequences


def match_chain_sequences(sequences: List[str], query_seq_records: List[SeqRecord], blast_path: str, perc_id_cutoff: float, cover_cutoff: float) -> List[str]:
    """
    Match the chain sequences against the query sequences using BLAST.
    
    Args:
        sequences (List[str]): List of chain sequences.
        query_seq_records (List[SeqRecord]): List of query SeqRecords.
        blast_path (str): Path to the BLAST executable.
        perc_id_cutoff (float): Minimum percentage identity required for a match.
        cover_cutoff (float): Minimum coverage percentage required for a match.
    
    Returns:
        List[str]: List of chain matches.
    """

    chain_matches = []

    for seq in sequences:

        SeqIO.write(SeqRecord(Seq(seq)), "temp2.fa", "fasta")

        blast_command = f'{blast_path} -query temp2.fa -db {query_seq_records}'
        blast_output = os.popen(blast_command).read()

        blast_records = BlastIO.parse(StringIO(blast_output), 'blast-text')

        for blast_record in blast_records:

            hit_def = blast_record.alignments[0].hit_def
            identity = blast_record.alignments[0].hsps[0].identities
            align_length = blast_record.alignments[0].hsps[0].align_length

            curr_pdb = hit_def.split('_')[0]
            curr_chn = hit_def.split('_')[1]

            perc_id = identity / align_length
            cov_perc = align_length / len(seq)

            if perc_id >= perc_id_cutoff and cov_perc >= cover_cutoff and curr_pdb != pdb:
                chain_matches.append(f'{curr_pdb}{curr_chn}')

    return chain_matches


def get_unbound_matches(matches: List[str], all_matches: List[str], seqres_path: str, pdbres: dict, perc_id_cutoff: float, cover_cutoff: float) -> str:
    """
    Get the unbound matches for the given matches.
    
    Args:
        matches (List[str]): List of chain matches.
        all_matches (List[str]): List of all chain matches.
        seqres_path (str): Path to the seqres fasta file.
        pdbres (dict): Dictionary containing PDB resolutions.
        perc_id_cutoff (float): Minimum percentage identity required for a match.
        cover_cutoff (float): Minimum coverage percentage required for a match.
    
    Returns:
        str: Unbound matches.
    """

    unbound_matches = ''

    for match in matches:

        pdb = match[:4]
        if pdbres.get(pdb, 0) <= 0 or pdbres[pdb] > res_cutoff:
            continue

        all_seqs = get_sequences(pdb, [match[4:]], seqres_objs)

        num_matches = 0
        mismatch_found = 0
        
        for seq in all_seqs:
            for i, _ in enumerate(all_matches):
                SeqIO.write(SeqRecord(Seq(seq)), 'temp_query.fa', 'fasta')
                if not is_sequence_in_matches('temp_query.fa', all_matches[i], seqres_path, perc_id_cutoff, cover_cutoff):
                    mismatch_found = 1
                    break
                num_matches += 1

        if mismatch_found == 0 and num_matches == len(all_matches[i]):
            unbound_matches += f' {pdb}'

    return unbound_matches.strip()

# Test the functions using your data_dir and blast_exe paths
data_dir = "/home/pierceb/databases/"
blast_exe = "/ibbr/ncbi-blast/bin/blastp"
match_sequences(data_dir, blast_exe)