import pandas as pd
import gzip

def create_tss_dictionary(gtf_file, output_csv):
    """
    Create a BioMart-style TSS dictionary from GENCODE GTF
    
    GTF format: chr, source, feature, start, end, score, strand, frame, attributes
    """
    
    genes = []
    
    # Read GTF (handle gzip)
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    mode = 'rt' if gtf_file.endswith('.gz') else 'r'
    
    with open_func(gtf_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type != 'gene':
                continue
            
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            
            # Parse attributes
            attr_dict = {}
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                key, value = attr.split(' ', 1)
                attr_dict[key] = value.strip('"')
            
            gene_id = attr_dict.get('gene_id', '')
            gene_name = attr_dict.get('gene_name', '')
            gene_type = attr_dict.get('gene_type', '')
            
            # TSS is start position (or end if minus strand)
            # Format to match H5 keys: chr{chrom}_{TSS}_{TSS+2}
            # Handle chromosome naming (GTF may have 'chr' prefix already)
            chr_name = chrom if chrom.startswith('chr') else f"chr{chrom}"
            
            if strand == '+':
                tss = start
                tss_enformer = f"{chr_name}_{tss}_{tss+2}"
            else:
                tss = end
                tss_enformer = f"{chr_name}_{tss-2}_{tss}"
            
            genes.append({
                'ensembl_gene_id': gene_id,
                'gene_name': gene_name,
                'chromosome': chrom,
                'TSS': tss,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_type': gene_type,
                'TSS_enformer_input': tss_enformer
            })
    
    # Create DataFrame and save
    df = pd.DataFrame(genes)
    df = df.drop_duplicates(subset=['ensembl_gene_id'])
    df.to_csv(output_csv, index=False)
    
    print(f"Created {output_csv} with {len(df)} genes")
    print(f"Columns: {list(df.columns)}")

if __name__ == "__main__":
    # Download GENCODE GTF first:
    # wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
    
    create_tss_dictionary('gencode.v49.annotation.gtf.gz', 'coding_gene_biomart_TSS_dictionary.csv')
