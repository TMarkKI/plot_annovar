import pandas as pd
import argparse


def calculate_depth(depth_file, coverage=0):
    results = dict()
    df = pd.read_csv(depth_file, sep="\t", header=None, names=["contig", "pos", "depth"])
    for contig, tmp in df.groupby("contig"):
        above = tmp[tmp["depth"] >= coverage].shape[0]
        results[contig] = {f"{coverage}x": above}
    depth_df = pd.DataFrame.from_dict(results, orient="index")
    return depth_df

def get_contig_length(depth_file):
    results = dict()
    df = pd.read_csv(depth_file, sep="\t", header=None, names=["contig", "pos", "depth"])
    for contig, tmp in df.groupby("contig"):
        results[contig] = {"length": tmp.shape[0]}
    len_df = pd.DataFrame.from_dict(results, orient="index")
    return len_df

def output(o_df):
    pd.o_df.to_csv()

def main(args):
    depth_file = args.depth
    master = get_contig_length(depth_file)
    for cov in [1,5,10,15,20,30]:
        master = master.join(calculate_depth(depth_file, coverage=cov))
        master[f"coverage % ({cov}x"] = master[f"{cov}x"] / master["length"] * 100
    #master["coverage % (10x)"] = master["10x"] / master["length"] * 100
    #master["coverage % (15x)"] = master["15x"] / master["length"] * 100
    #master["coverage % (20x)"] = master["20x"] / master["length"] * 100
    #master["coverage % (30x)"] = master["30x"] / master["length"] * 100
    print(master)
    master.to_csv(args.output, sep="\t")
    #output(master) = args.output

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--depth")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    main(args)
