import pandas as pd
import glob
import sys


def parse_diamond_gene_column(diamond_df):
    """
    Parse the gene_name column of a diamond dataframe.
    
    Parameters:
    - diamond_df (pd.DataFrame): DataFrame containing diamond output data.
    
    Returns:
    - pd.DataFrame: Updated DataFrame with the gene_name parsed.
    """
    # parse the gene_name column
    diamond_df[['tmp', 'hit', 'gene']] = diamond_df['gene_name'].str.split('|', expand=True)
    diamond_df[['gene', 'tmp2']] = diamond_df['gene'].str.split('_', expand=True)
    
    # drop unnecessary columns
    diamond_df = diamond_df.drop(columns=['tmp', 'tmp2', 'hit'])
    
    return diamond_df

def main():
    # get arguments from command line
    output_file = sys.argv[1]
    orthomerged_file = sys.argv[2]
    raw_files = sys.argv[3:]

    # diamond blastx results column names 
    diamond_colnames = ["contig_name", "gene_name", "pident", "length",
                        "mismatches", "gap", "qstart", "qend", "sstart", "send",
                        "evalue", "bitscore"]

    # read in orthomerged diamond results
    orthomerged = pd.read_csv(orthomerged_file, sep='\t', names=diamond_colnames)
    orthomerged = parse_diamond_gene_column(orthomerged)

    # read in raw assembly diamond results
    dfs = [pd.read_csv(f, sep='\t', names=diamond_colnames) for f in raw_files]
    raw = pd.concat(dfs, ignore_index=True)
    raw = parse_diamond_gene_column(raw)

    # determine which genes are annotated in raw assemblies but not in orthomerged
    newbies = raw[~raw['gene'].isin(orthomerged['gene'])]
    newbies = newbies.sort_values('bitscore', ascending=False).groupby('contig_name').head(1)
    newbies = newbies[~newbies['contig_name'].isin(orthomerged['contig_name'])]

    # drop duplicate contig names 
    unique_contig_names = newbies['contig_name'].drop_duplicates()
    # write a file containing the new contig_names 
    unique_contig_names.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
