#!/usr/bin/python3

import sys
import os

def main():
    file = sys.argv[1]
    if file == "":
        sys.exit("usage: get_cases.py antibody_antigen_triplet_list")

    pdb_res_file = "/Users/bpierce/databases/pdb_resolutions_8_2020.txt"
    current_cases = "bm5.5.txt"

    curr_bm = []

    file_lines = read_file_to_list(file)
    file_lines = [line.rstrip() for line in file_lines]

    possible_cases = {}
    for line in file_lines:
        fields = line.split("\t")
        if fields[1] != "" and fields[2] != "":
            possible_cases[fields[0]] = fields[1] + "\t" + fields[2]

    res_lines = read_file_to_list(pdb_res_file)
    res_lines = [line.rstrip() for line in res_lines]

    pdbres = {}
    no_res_yet = True
    for i, line in enumerate(res_lines):
        if line.startswith("IDCODE"):
            i += 2
            no_res_yet = False
        if no_res_yet:
            continue
        pdb, stuff, res = line.split(" ")
        pdbres[pdb.lower()] = res

    bm5_lines = read_file_to_list(current_cases)
    bm5_lines = [line.rstrip() for line in bm5_lines]

    bm5_matches = 0
    for line in bm5_lines:
        if line in possible_cases:
            bm5_matches += 1
            curr_bm.append(line)
    print(str(bm5_matches) + " out of " + str(len(bm5_lines)) + " BM5 matches")

    all_comps = list(possible_cases.keys())
    for comp in all_comps:
        mabs, ags = possible_cases[comp].split("\t")
        mab_comps = mabs.split(" ")
        top_mab = mab_comps[0]
        top_mab_res = pdbres[mab_comps[0]]
        for mab_comp in mab_comps[1:]:
            if pdbres[mab_comp] < top_mab_res:
                top_mab = mab_comp
                top_mab_res = pdbres[mab_comp]
        ag_comps = ags.split(" ")
        top_ag = ag_comps[0]
        top_ag_res = pdbres[ag_comps[0]]
        for ag_comp in ag_comps[1:]:
            if pdbres[ag_comp] < top_ag_res:
                top_ag = ag_comp
                top_ag_res = pdbres[ag_comp]
        possible_cases[comp] = top_mab + "\t" + top_ag

    new_cases = 0
    new_case_comps = []
    for comp in all_comps:
        match_found = False
        for bm in curr_bm:
            if comp == bm or possible_cases[comp] == possible_cases[bm]:
                match_found = True
                break
        if not match_found:
            curr_bm.append(comp)
            new_case_comps.append(comp)
            new_cases += 1
    print(str(new_cases) + " new cases found!")
    for case in new_case_comps:
        print(case + "\t" + possible_cases[case])


def read_file_to_list(filepath):
    if not os.path.exists(filepath):
        sys.exit("unable to open file: " + filepath)
    with open(filepath, 'r') as file:
        return file.readlines()


if __name__ == "__main__":
    main()