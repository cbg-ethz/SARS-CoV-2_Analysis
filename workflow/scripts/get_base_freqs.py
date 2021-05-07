#!/usr/bin/env python3

#################################################################################################################################################################
###     This program extracts position specific nucleotide frequencies from all vcf files in the given file hierarchy.
###     *  Can exclude deletion read counts by setting "noDel=True".
###     *  Columns in vcf file: CHROM    POS	ID	REF	ALT	QUAL	FILTER	INFO
###     *  Example line (ShoRaH):
###        NC_045512.2	241	.	T	C	4.77121	PASS	Freq1=1;Freq2=1;Freq3=0;Post1=1.324;Post2=1.0286;Post3=0;Fvar=2793;Rvar=2818;Ftot=2794;Rtot=2819;Pval=1;Qval=1
###     *  Example line (LoFreq):
#################################################################################################################################################################

import sys
import os
import glob

import pandas as pd
from tqdm import tqdm

INDEL_COLS = 10


def get_base_freqs(
    fname_list,
    fname_samples,
    out_file,
    delThreshold=0,
    inThreshold=0,
    indelCoverageThreshold=0,
):

    df_samples = pd.read_csv(fname_samples)
    usedSamples = df_samples["accession"].tolist()
    print(f"Parsing {len(usedSamples)}/{len(fname_list)} samples")

    cols = [
        "SAMPLE",
        "POS",
        "REF_BASE",
        "ALT_BASE",
        "(ADJUSTED_)READ_COUNT",
        "A_freq",
        "C_freq",
        "G_freq",
        "T_freq",
        "DEL_freq",
        "IN_freq",
    ]
    cols.extend([f"IN{i}_freq" for i in range(INDEL_COLS)])
    print(out_file)
    out_stream = open(out_file, "w")
    out_stream.write("\t".join(cols) + "\n")

    # all files matching this pattern are processed
    # pattern = "samples/*/*/variants/SNVs/snvs.vcf"
    pattern = os.path.join(os.path.dirname(fname_list[0]), "*.vcf")
    fileList = glob.glob(pattern)

    for file in tqdm(fileList):
        sampleName = os.path.basename(file)
        # Remove file extensions to get sample name
        if sampleName.endswith(".vcf"):
            sampleName = sampleName[:-4]
        if sampleName.startswith("snvs_"):
            sampleName = sampleName[5:]

        if sampleName not in usedSamples:
            continue

        baseCounts = {}  # pos -> base counts, refBase, altBase
        pos = 0
        with open(file) as f:
            lines = f.read().splitlines()
            for line in lines:
                # Skip vcf header lines
                if line.startswith("#"):
                    continue

                # Get the base read counts from the line
                elems = line.split("\t")

                assert int(elems[1]) >= pos, "VCF file needs to be sorted"

                pos = int(elems[1])
                refBase = elems[3]
                altBase = elems[4]

                # LoFreq input
                if "DP4" in elems[7]:
                    info = elems[7].split(";")
                    # Counts for ref-forward bases, ref-reverse, alt-forward,
                    #   and alt-reverse bases
                    counts = [
                        int(j)
                        for i in info
                        for j in i[4:].split(",")
                        if i.startswith("DP4=")
                    ]
                    Fvar = counts[2]
                    Rvar = counts[3]
                    Ftot = Fvar + counts[0]
                    Rtot = Rvar + counts[1]
                # ShoRaH input
                else:
                    Fvar = elems[7].split("Fvar=")[1].split(";")[0]
                    Rvar = elems[7].split("Rvar=")[1].split(";")[0]
                    Ftot = elems[7].split("Ftot=")[1].split(";")[0]
                    Rtot = elems[7].split("Rtot=")[1].split(";")[0]

                total = int(Ftot) + int(Rtot)
                altTotal = int(Fvar) + int(Rvar)

                # Deletion in LoFreq
                if len(refBase) > len(altBase):
                    for i_pos, i_base in enumerate(refBase[1:], 1):
                        try:
                            baseCounts[pos + i_pos]["-"] += altTotal
                            baseCounts[pos + i_pos]["total"] += total
                        except KeyError:
                            baseCounts[pos + i_pos] = {
                                "A": 0,
                                "C": 0,
                                "G": 0,
                                "T": 0,
                                "-": altTotal,
                                "+0": 0,
                                "total": total,
                                "ref": i_base,
                                "alt": ["-"],
                            }
                    continue
                # Mutation, Insertion, or Deletion in SHoRaH
                else:
                    try:
                        baseCounts[pos]["ref"] = refBase
                        baseCounts[pos]["alt"].append(altBase)
                    # New position
                    except KeyError:
                        baseCounts[pos] = {
                            "A": 0,
                            "C": 0,
                            "G": 0,
                            "T": 0,
                            "-": 0,
                            "+0": 0,
                            "total": total,
                            "ref": refBase,
                            "alt": [altBase],
                        }
                        # Insertion (LoFreq and ShoRaH)
                        if len(altBase) > len(refBase):
                            altBase = "+0"
                        baseCounts[pos][altBase] = altTotal
                    # Adjust read counts for already present position
                    else:
                        # Insertion (LoFreq and ShoRaH)
                        if len(altBase) > len(refBase):
                            in_no = sum(
                                [len(i) > 1 for i in baseCounts[pos]["alt"][:-1]]
                            )
                            # If there is already an insetion, add another col
                            altBase = f"+{in_no}"
                            baseCounts[pos][altBase] = 0

                        baseCounts[pos][altBase] += altTotal
                        baseCounts[pos]["total"] += total

                # Checks that read counts in refering to same position make sense
                if baseCounts[pos][altBase] != 0:
                    # Only deletion events: can overlap
                    if all(["-" == i for i in baseCounts[pos]["alt"]]):
                        pass
                    # Different events: shouldnt overlap
                    elif baseCounts[pos][altBase] != altTotal:
                        print(
                            f"Warning: pos {pos}: different alt   base count: "
                            f"{baseCounts[pos][altBase]} = {altTotal}"
                        )

                # DEPRECATED for LoFreq input
                # if total != baseCounts[pos]['total'] != 0:
                #     print("Warning: pos {}: different total base count: {} = {}" \
                #         .format(pos, baseCounts[pos]['total'], total))

        # Assign difference between total read count and sum of alt read counts as ref count
        for pos, pos_data in sorted(baseCounts.items()):
            # Get number of insertions
            in_list = [i for i in pos_data["alt"] if len(i) > 1]
            in_no = len(in_list)
            # Get sum of alternative read counts
            altCounts = (
                pos_data["A"]
                + pos_data["C"]
                + pos_data["G"]
                + pos_data["T"]
                + pos_data["-"]
            )
            for in_i in range(in_no):
                altCounts += pos_data[f"+{in_i}"]
            pos_data[pos_data["ref"]] = pos_data["total"] - altCounts

            # Check that read counts make sense
            if pos_data[pos_data["ref"]] < 0:
                print(f"ERROR: NEGATIVE REF COUNT AT POSITION {pos} IN {file}:")
                print(baseCounts[pos])

            useDel = True
            useIn = [True for i in range(in_no)]
            # Check if indels pass read coverage threshold
            if pos_data["total"] <= indelCoverageThreshold:
                useDel = False
                useIn = [False for i in range(in_no)]
            else:
                # Check if deletion pass the relative presence threshold
                if pos_data["-"] / pos_data["total"] <= delThreshold:
                    useDel = False
                # Check if insertion pass the relative presence threshold
                for in_i in range(in_no):
                    if pos_data[f"+{in_i}"] / pos_data["total"] <= inThreshold:
                        useIn[in_i] = False

            # If deletion is not considered, update total considered read counts
            if not useDel:
                pos_data["total"] -= pos_data["-"]
                pos_data["-"] = 0
                pos_data["alt"] = [i for i in pos_data["alt"] if i != "-"]
            # If insertion is not considered, update total considered read counts
            # If insertion is not considered, update total considered read counts
            for in_i, in_flag in enumerate(useIn):
                if not in_flag:
                    pos_data["total"] -= pos_data[f"+{in_i}"]
                    pos_data[f"+{in_i}"] = 0
                    pos_data["alt"] = [i for i in pos_data["alt"] if i != in_list[in_i]]

            # Skip position if no reads are left
            if pos_data["total"] == 0:
                continue
            # Skip position if no SNP reads are left
            if pos_data[pos_data["ref"]] == pos_data["total"]:
                continue

            # compute variant frequencies from base counts
            new_line = [
                sampleName,
                pos,
                pos_data["ref"],
                ",".join(pos_data["alt"]),
                pos_data["total"],
                pos_data["A"] / pos_data["total"],
                pos_data["C"] / pos_data["total"],
                pos_data["G"] / pos_data["total"],
                pos_data["T"] / pos_data["total"],
                pos_data["-"] / pos_data["total"],
            ]
            new_line.extend((INDEL_COLS + 1) * [0])
            # Check if new insertion line was added
            in_freq_sum = 0
            valid_in_no = sum(useIn)
            if 0 < valid_in_no <= INDEL_COLS:
                # Iterate over second, third/ fourth... insertion
                in_col = 0
                for in_i, in_i_flag in enumerate(useIn):
                    if in_i_flag:
                        in_freq = pos_data[f"+{in_i}"] / pos_data["total"]
                        new_line[11 + in_col] = in_freq
                        in_freq_sum += in_freq
                        in_col += 1
            elif valid_in_no > INDEL_COLS:
                print(f"ERROR: >{INDEL_COLS} INSERTIONS AT POS {pos} IN {file}:")
                print(baseCounts[pos])
            new_line[10] = in_freq_sum

            out_stream.write("\t".join([str(i) for i in new_line]) + "\n")

    out_stream.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname_list", type=str, nargs="+")
    parser.add_argument("-s", "--fname_samples", type=str, default="")
    parser.add_argument("-o", "--out_file", type=str, default="test_out.csv")
    parser.add_argument("-d", "--delThreshold", type=float, default=0)
    parser.add_argument("-i", "--inThreshold", type=float, default=0)
    parser.add_argument("-c", "--indelCoverageThreshold", type=int, default=0)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    if "snakemake" in globals():
        get_base_freqs(
            snakemake.input["fname_list"],
            snakemake.input["fname_samples"],
            snakemake.output["fname"],
            delThreshold=snakemake.config["params"]["deletion_threshold"],
            inThreshold=snakemake.config["params"]["deletion_threshold"],
            indelCoverageThreshold=snakemake.config["params"][
                "deletion_coverage_threshold"
            ],
        )
    else:
        import argparse

        args = parse_args()
        if not args.fname_samples:
            import tempfile

            def clean_fname(fname):
                if fname.endswith(".vcf"):
                    fname = fname[:-4]
                if fname.startswith("snvs_"):
                    fname = fname[5:]
                return fname

            samples = tempfile.NamedTemporaryFile(delete=False)
            sample_list = [clean_fname(os.path.basename(i)) for i in args.fname_list]
            samples.write(str.encode("accession\n" + "\n".join(sample_list)))
            samples.close()
            args.fname_samples = samples.name
        get_base_freqs(
            args.fname_list,
            args.fname_samples,
            args.out_file,
            delThreshold=args.delThreshold,
            inThreshold=args.inThreshold,
            indelCoverageThreshold=args.indelCoverageThreshold,
        )
