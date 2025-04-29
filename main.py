import polars as pl
from pathlib import Path


def is_XorY(name:str)->bool:
    """Checks whether a chromosome is X or Y chromosome by its name."""
    return 'x' in name.lower() or 'y' in name.lower()

def extract_line_metadata(line:str) -> dict:

    # get chromosome name
    chr_name = line.split('\t')[0]
    # get metadata
    metadata = line.split('\t')[-1].split('; ')
    headers = [m.split()[0] for m in metadata]
    values =  [m.split()[1].replace('"','') for m in metadata]
    info = dict(zip(headers,values))
    # add chr name to dict
    info['chr'] = chr_name

    return info
    

def parse_gtf(gtf_file:str):
    
    tsv_file = gtf_file.replace('gtf','tsv')
    with open(tsv_file,'w') as tsv:
        tsv.write("chr\tgene_id\tgene_name\n") # header line
        with open(gtf_file,'r') as file:
            for line in file:
                if line.startswith('##'): # ignore headers
                    continue
                elif line.startswith('\n'): # last line
                    pass
                else:
                    chr_name = line.split('\t')[0]
                    if is_XorY(chr_name):
                        info = extract_line_metadata(line=line)
                        to_write = f"{info['chr']}\t{info['gene_id']}\t{info['gene_name']}\n"
                        tsv.write(to_write)


    return


def main():
    parse_gtf('gencode.v47.annotation.gtf')



if __name__ == "__main__":
    main()
