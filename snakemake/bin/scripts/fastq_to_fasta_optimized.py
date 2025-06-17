import gzip
import sys

def fastq_to_fasta(input_fastq, output_fasta):
    with gzip.open(input_fastq, 'rt') if input_fastq.endswith('.gz') else open(input_fastq, 'r') as fq, open(output_fasta, 'w') as fa:
        while True:
            header = fq.readline().strip()
            if not header:
                break
            seq = fq.readline().strip()
            fq.readline()  # skip plus line
            fq.readline()  # skip quality line
            fa.write(f">{header[1:]}\n{seq}\n")

if __name__ == "__main__":
    input_fastq = sys.argv[1]
    output_fasta = sys.argv[2]
    fastq_to_fasta(input_fastq, output_fasta)
