import sys

def parse_diamond_gene_column(line):
    """
    Parse the gene_name column of a diamond dataframe line.
    
    Parameters:
    - line (str): A single line from the diamond output file
    
    Returns:
    - tuple: Parsed values (contig_name, gene)
    """
    # Splitting the line
    values = line.strip().split('\t')
    
    contig_name = values[0]
    gene_name_parts = values[1].split('|')
    gene = gene_name_parts[-1].split('_')[0]
    
    return contig_name, gene

def main():
    output_file = sys.argv[1]
    orthomerged_file = sys.argv[2]
    raw_files = sys.argv[3:]

    orthomerged_genes = set()
    with open(orthomerged_file, 'r') as f:
        for line in f:
            _, gene = parse_diamond_gene_column(line)
            orthomerged_genes.add(gene)

    raw_contig_genes = {}
    for raw_file in raw_files:
        with open(raw_file, 'r') as f:
            for line in f:
                contig_name, gene = parse_diamond_gene_column(line)
                if gene not in orthomerged_genes:
                    # Use bitscore to decide which gene to associate with a contig.
                    bitscore = float(line.strip().split('\t')[-1])
                    if contig_name not in raw_contig_genes or bitscore > raw_contig_genes[contig_name][1]:
                        raw_contig_genes[contig_name] = (gene, bitscore)

    # Filter out contig_names already in orthomerged
    contig_output = [contig for contig in raw_contig_genes.keys() if contig not in orthomerged_genes]

    # Save unique contig names to output file
    with open(output_file, 'w') as f:
        for contig in contig_output:
            f.write(contig + '\n')

if __name__ == "__main__":
    main()
