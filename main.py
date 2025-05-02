import polars as pl
from argparse import ArgumentParser


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
    """Parse gtf file into  tsv file."""
    tsv_file = gtf_file.replace('gtf','tsv')
    with open(tsv_file,'w') as tsv:
        tsv.write("chr\tgene_id\tgene_name\tline\n") # header line
        with open(gtf_file,'r') as file:
            for i,line in enumerate(file):
                if line.startswith('##'): # ignore headers
                    continue
                elif line.startswith('\n'): # last line
                    pass
                else:
                    chr_name = line.split('\t')[0]
                    if is_XorY(chr_name):
                        info = extract_line_metadata(line=line)
                        to_write = f"{info['chr']}\t{info['gene_id']}\t{info['gene_name']}\t{i}\n"
                        tsv.write(to_write)
    
    return


def find_PAR(tsv_file:str,*,selector:str='gene_id', sort:bool=True)->pl.DataFrame:
    """find regions to add PAR suffix"""
    df = pl.read_csv(tsv_file,separator="\t")
    unq = df.unique(pl.all().exclude('line'))
    gene_ids = unq.select(selector).unique(selector).to_series().to_list()

    # find repeated gene_id for X and Y chromosomes
    repeated = []
    for gene_id in gene_ids:
        filtered = unq.filter(pl.col(selector)==gene_id)
        if filtered.height > 1: # means there is a repeat (exist both on x and y chr)
            repeated.append(gene_id)
    
    # find line to change gene_id (add suffix) 
    to_change = df.filter(pl.col(selector).is_in(repeated),pl.col('chr').is_in(['Y','chrY']))
    # add the suffix
    modified = to_change.with_columns(pl.col(selector).add('_PAR_Y').alias('gene_id_par'))

    if sort:
        modified = modified.sort('line')

    return modified

def replace_PAR(gtf_file:str)->None:
    """Replace gene id with gene id PAR."""

    tsv_file = gtf_file.replace('gtf','tsv') # change suffix to tsv
    parse_gtf(gtf_file=gtf_file) # generates a tsv file 
    
    # modify the gene ids
    modified = find_PAR(tsv_file)
    line_numbers = modified.select('line').to_series().to_list() # line to be altered

    par_gtf_file = gtf_file.replace('.gtf','.par.gtf')
    # IO operations
    with open(par_gtf_file,'w') as par: # open PAR gtf to write
        with open(gtf_file,'r') as gtf: # open gtf to read
            for i,line in enumerate(gtf):
                if i in line_numbers:
                    # if line is to be changed
                    row = modified.filter(pl.col('line')==i)
                    print(row)
                    gene_id = row.select('gene_id').item()  # retrieve original gene id
                    gene_id_par = row.select('gene_id_par').item() # retrieve PAR gene id
                    new_line = line.replace(gene_id,gene_id_par) # replace the gene id with PAR gene id
                else:
                    new_line = line
                # write the line
                par.write(new_line)

    return


def main():
    replace_PAR('madeup.gtf')

if __name__ == "__main__":
    main()
