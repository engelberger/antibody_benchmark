import subprocess

PERCID_CUTOFF = 0.98
COV_CUTOFF = 0.8
RES_CUTOFF = 3.25

data_dir = "/home/pierceb/databases/"
blast_exe = "/ibbr/ncbi-blast/bin/blastp"

sabdab_list = data_dir + "all_antibody_structures_8_2020.txt"
seqres_file = data_dir + "pdb_seqres_8_2020.fa"
pdb_res_file = data_dir + "pdb_resolutions_8_2020.txt"

def run_blast(query_file, db_file):
    blast_cmd = [blast_exe, "-query", query_file, "-db", db_file]
    result = subprocess.run(blast_cmd, capture_output=True, text=True)
    return result.stdout

def get_chain_sequence(pdb, chain):
    with open(seqres_file, "r") as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if lines[i].startswith(">") and lines[i][1:5] == pdb and lines[i][6] == chain:
                seq = ""
                i += 1
                while i < len(lines) and not lines[i].startswith(">"):
                    seq += lines[i].strip()
                    i += 1
                return seq
    return ""

def main(test_pdb=""):
    with open(pdb_res_file, "r") as file:
        res_lines = file.readlines()

    pdbres = {}
    no_res_yet = True
    for line in res_lines:
        if line[0:6] == "IDCODE":
            no_res_yet = False
        if no_res_yet:
            continue
        pdb, _, res = line.split()
        pdbres[pdb.lower()] = float(res)

    with open(sabdab_list, "r") as file:
        list_lines = file.readlines()

    for line in list_lines:
        fields = line.split("\t")
        pdb = fields[0]
        if pdb in pdb_processed:
            continue
        pdb_processed.add(pdb)
        hchain = fields[1]
        lchain = fields[2]
        model = int(fields[3])
        ag_chain = fields[4]
        ag_type = fields[5]
        if model == 1:
            continue
        res = float(fields[16])
        method = fields[17]
        if res > RES_CUTOFF or not method.startswith("X-RAY"):
            continue
        if pdbres[pdb.lower()] > RES_CUTOFF or pdbres[pdb.lower()] <= 0:
            continue
        if hchain == "NA" or lchain == "NA" or ag_chain == "NA":
            continue
        if not ag_type.lower().startswith("protein"):
            continue
        if test_pdb != "" and pdb != test_pdb.lower():
            continue

        # Get the chain sequences
        ag_chains = ag_chain.split("|")
        ag_chains = [chain.strip() for chain in ag_chains]

        # Get the sequences
        seqs = []
        chains = [hchain, lchain] + ag_chains

        chain_matches = []
        for chain in chains:
            chain_matches.append("")
            seq = get_chain_sequence(pdb, chain)
            if not seq:
                print(f"Error finding sequence for {pdb} {chain}")
                continue

            with open("temp.fa", "w") as file:
                file.write(f"> {pdb} {chain}\n")
                file.write(f"{seq}\n")

            # Run BLAST
            blout = run_blast("temp.fa", seqres_file)
            blastout_lines = blout.split("\n")
            
            tot_num = 0
            match_num = 0
            query_offset = -1
            subj_offset = -1
            curr_chn = ""
            curr_pdb = ""
            for line in blastout_lines:
                if line.startswith(">"):
                    if "mol:protein" in line:
                        _, pdb_chn, _, _ = line.split("_")
                    else:
                        pdb_chn = line.split(" ")[0][1:] + line.split(" ")[1][1]
                    curr_pdb = pdb_chn[:4]
                    curr_chn = pdb_chn[4]
                elif "Identities =" in line:
                    match_num, tot_num = map(int, line.split()[2].split("/"))
                    perc_id = match_num / tot_num
                    if perc_id < PERCID_CUTOFF:
                        continue
                    cov_perc = tot_num / len(seq)
                    if cov_perc < COV_CUTOFF:
                        continue
                    if curr_pdb == pdb:
                        continue

                    # Looks like we have a winner
                    if chain_matches[-1]:
                        chain_matches[-1] += "\t"
                    chain_matches[-1] += curr_pdb + curr_chn
        
        # Check if at least one chain does not match
        if all(not chain_match for chain_match in chain_matches):
            continue

        # Check for non-antibody or non-antigen chains
        hchain_matches = chain_matches[0].split("\t") if chain_matches[0] else []
        lchain_matches = chain_matches[1].split("\t") if chain_matches[1] else []
        ub_mab_matches = ""
        mab_checked = set()

        for hchns in hchain_matches:
            hpdb = hchns[:4]
            hchn = hchns[4]
            if hpdb in mab_checked:
                continue
            mab_checked.add(hpdb)
            if pdbres[hpdb] > RES_CUTOFF or pdbres[hpdb] <= 0:
                continue
            hchains = ""
            lchains = ""
            for lchns in lchain_matches:
                lpdb = lchns[:4]
                lchn = lchns[4]
                if lpdb == hpdb:
                    lchains += lchn
            if not lchains:
                continue
            for hchns2 in hchain_matches:
                hpdb2 = hchns2[:4]
                hchn2 = hchns2[4]
                if hpdb2 == hpdb:
                    hchains += hchn2
            mabchains = hchains + lchains
            num_matches = 0
            mismatch_found = False

            for line in seqres_lines:
                if line.startswith(">"):
                    pdb2 = line[1:5]
                    chain = line[6]
                    if pdb2 != hpdb:
                        continue
                    num_matches += 1
                    if chain in mabchains:
                        continue
                    else:
                        mismatch_found = True
                        break

            if mismatch_found:
                continue
            if num_matches != len(mabchains):
                print("Error: weird chain matches")
                continue
            if ub_mab_matches:
                ub_mab_matches += " "
            ub_mab_matches += hpdb

        ub_ag_matches = ""
        ag_checked = set()
        tmp_matches = chain_matches[2].split("\t") if chain_matches[2] else []

        for agchns in tmp_matches:
            agpdb = agchns[:4]
            agchn = agchns[4]
            if agpdb in ag_checked:
                continue
            ag_checked.add(agpdb)
            if pdbres[agpdb] > RES_CUTOFF or pdbres[agpdb] <= 0:
                continue
            allagchns = ""
            for j in range(2, len(chain_matches)):
                curr_agchns = ""
                tmp_matches = chain_matches[j].split("\t") if chain_matches[j] else []
                for agchns2 in tmp_matches:
                    agpdb2 = agchns2[:4]
                    agchn2 = agchns2[4]
                    if agpdb2 == agpdb:
                        curr_agchns += agchn2
                if not curr_agchns:
                    allagchns = ""
                    break
                allagchns += curr_agchns
            if not allagchns:
                continue
            num_matches = 0
            mismatch_found = False
            for line in seqres_lines:
                if line.startswith(">"):
                    pdb2 = line[1:5]
                    chain = line[6]
                    if pdb2 != agpdb:
                        continue
                    num_matches += 1
                    if chain in allagchns:
                        continue
                    else:
                        mismatch_found = True
                        break
            if mismatch_found:
                continue
            if num_matches != len(allagchns):
                print("This is weird")
            if ub_ag_matches:
                ub_ag_matches += " "
            ub_ag_matches += agpdb

        print(f"{pdb}\t{ub_mab_matches}\t{ub_ag_matches}")


if __name__ == "__main__":
    pdb_processed = set()
    main()