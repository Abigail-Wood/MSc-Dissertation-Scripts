import argparse
import pandas as pd

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-q','--queries', required = True,
        help='Input list of query names with column header "query".')
    parser.add_argument(
        '-f','--fasta', required = True,
        help='Input Fasta file containing protein sequences.')
    parser.add_argument(
        '-o','--output', required = True,
        help='Output name for Fasta file containing protein sequences.')
    return parser.parse_args()


def extract_fasta_sequence(queries, fasta, output):
    """ Extract T. trichiura protein sequence from a Fasta file.
    Arguments:
    queries - file output with list of query names. 
    fasta - fasta file source.
    """
    df = pd.read_csv(queries)
    query_list = list(df['query'])
    collect_sequence = False
    with open(fasta, 'r') as f, open(output,'w+') as o:
        for line in f.readlines():
            if line.startswith(">") and line[1:-1] in query_list:
                o.write("\n" + f"{line}") 
                collect_sequence = True
            elif collect_sequence == True:
                if line.startswith(">"):
                    collect_sequence = False
                else:
                    o.write((line.strip()))
    return 0

def main():
    args = parse_args()
    extract_fasta_sequence(args.queries, args.fasta, args.output)

if __name__ == "__main__":
    main()
