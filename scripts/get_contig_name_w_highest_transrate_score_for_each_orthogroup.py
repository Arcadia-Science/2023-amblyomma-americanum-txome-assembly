import argparse
from collections import defaultdict

def read_orthogroups(filepath):
    orthogroups = defaultdict(list)

    with open(filepath, 'r') as file:
        for line in file:
            orthogroup, contig_names = line.strip().split(": ")
            for contig in contig_names.split():
                orthogroups[orthogroup].append(contig)

    return orthogroups

def process_contigs_csv(filepaths, orthogroups):
    contig_scores = {}
    contig_coverage = {}

    for filepath in filepaths:
        with open(filepath, 'r') as file:
            header = file.readline().strip().split(',')
            score_idx = header.index("score")
            coverage_idx = header.index("coverage")

            for line in file:
                parts = line.strip().split(',')
                contig_scores[parts[0]] = float(parts[score_idx])
                contig_coverage[parts[0]] = float(parts[coverage_idx])

    selected_contigs = []

    for contigs in orthogroups.values():
        # right now, all isoseq contig names start with "SRR". 
        # in the future we can have a renaming script that that adds isoseq to the contig name, and solve for that
        isoseq_contigs = [contig for contig in contigs if contig.startswith("SRR")]
        if isoseq_contigs:
            selected_contigs.extend(isoseq_contigs)
        else:
            contigs.sort(key=lambda x: (contig_scores.get(x, 0), contig_coverage.get(x, 0)), reverse=True)
            selected_contigs.append(contigs[0])

    return selected_contigs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process orthogroups and contigs. For each orthogroup, select the contig with the desired criteria.")
    parser.add_argument("output_path", help="Path for the output file")
    parser.add_argument("orthogroups_path", help="Path to the orthogroups file")
    parser.add_argument("contigs_csv_paths", nargs='+', help="Paths to the contigs.csv files")
    args = parser.parse_args()

    orthogroups = read_orthogroups(args.orthogroups_path)
    selected_contigs = process_contigs_csv(args.contigs_csv_paths, orthogroups)

    with open(args.output_path, 'w') as out_file:
        for contig in selected_contigs:
            out_file.write(f"{contig}\n")
