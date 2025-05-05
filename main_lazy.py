import polars as pl
import argparse
from typing import Literal

class FileTypeCouldNotBeInferredError(Exception):
    def __init__(self, message="The file type of the supposed gtf or gff3 file could not be inferred."):
        self.message = message
        super().__init__(self.message)


def is_XorY(name:str)->bool:
    """Checks whether a chromosome is X or Y chromosome by its name."""
    return 'x' in name.lower() or 'y' in name.lower()

def infer_file_type(fname:str)-> Literal['gtf','gff3']:
    '''Infer the file type of the gtf or gff file.'''
    file_type = None
    with open(fname,'r') as file:
        for line in file:
            if line.startswith('##'): # ignore headers
                continue
            else:
                attributes = line.split('\t')[-1]
                if "=" in attributes:
                    file_type = 'gff3'
                elif '"' in attributes:
                    file_type = 'gtf'
                else:
                    raise FileTypeCouldNotBeInferredError("The file type of the supposed gtf or gff3 file could not be inferred.")
                
    return file_type


def _extract_line_metadata_gtf(line:str) -> dict:

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


def _extract_line_metadata_gff3(line:str) -> dict:

    # get chromosome name
    chr_name = line.split('\t')[0]
    # get metadata
    metadata = line.split('\t')[-1].split(';') # no spaces for gff3
    headers = [m.split('=')[0] for m in metadata]
    values =  [m.split('=')[1] for m in metadata]
    info = dict(zip(headers,values))
    # add chr name to dict
    info['chr'] = chr_name

    # rename some keys
    info['gene_id'] = info.pop('ID')
    info['gene_name'] = info.pop('Name')

    return info


def extract_line_metadata(line:str, file_type:Literal['gff3','gtf']) -> dict:

    if file_type == 'gtf':
        info = _extract_line_metadata_gtf(line=line)
    elif file_type == 'gff3':
        info = _extract_line_metadata_gff3(line=line)
    else:
        raise ValueError(f'File type should be either `gff3` or `gtf` not `{file_type}`')

    return info
    

def parse_gxf(gxf_file:str):
    """Parse gxf file into  tsv file."""
    extension = gxf_file.split(".")[-1]
    tsv_file = gxf_file.replace(extension,'tsv')
    file_type = infer_file_type(fname=gxf_file)
    with open(tsv_file,'w') as tsv:
        tsv.write("chr\tgene_id\tgene_name\tline\n") # header line
        with open(gxf_file,'r') as file:
            for i,line in enumerate(file):
                if line.startswith('##'): # ignore headers
                    continue
                elif line.startswith('\n'): # last line
                    pass
                else:
                    chr_name = line.split('\t')[0]
                    if is_XorY(chr_name):
                        info = extract_line_metadata(line=line,file_type=file_type)
                        to_write = f"{info['chr']}\t{info['gene_id']}\t{info['gene_name']}\t{i}\n"
                        tsv.write(to_write)
    
    return


def find_par(tsv_file:str,*,selector:str='gene_id', sort:bool=True)->pl.LazyFrame:
    """find regions to add PAR suffix"""
    lf = pl.scan_csv(tsv_file,separator="\t")
    unq = lf.unique(pl.all().exclude('line'))
    gene_ids = unq.select(selector).unique(selector).collect().to_series().to_list()

    # find repeated gene_id for X and Y chromosomes
    repeated = []
    for gene_id in gene_ids:
        filtered = unq.filter(pl.col(selector)==gene_id).collect()
        if filtered.height > 1: # means there is a repeat (exist both on x and y chr)
            repeated.append(gene_id)
    
    # find line to change gene_id (add suffix) 
    to_change = lf.filter(pl.col(selector).is_in(repeated),pl.col('chr').is_in(['Y','chrY']))
    # add the suffix
    modified = to_change.with_columns(pl.col(selector).add('_PAR_Y').alias('gene_id_par'))

    if sort:
        modified = modified.sort('line')

    return modified

def add_par(gxf_file:str)->None:
    """Add PAR suffix to gene_id."""
    extension = gxf_file.split(".")[-1]
    tsv_file = gxf_file.replace(extension,'tsv')# change suffix to tsv
    parse_gxf(gxf_file=gxf_file) # generates a tsv file 
    
    # modify the gene ids
    modified = find_par(tsv_file)
    print(f"Found {modified.collect().height}",sep="\t")
    line_numbers = modified.select('line').collect().to_series().to_list() # line to be altered

    par_gxf_file = gxf_file.replace(extension,f'par.{extension}')
    # IO operations
    with open(par_gxf_file,'w') as par: # open PAR gxf to write
        with open(gxf_file,'r') as gxf: # open gxf to read
            for i,line in enumerate(gxf):
                if i in line_numbers:
                    # if line is to be changed
                    row = modified.filter(pl.col('line')==i)
                    gene_id = row.select('gene_id').collect().item()  # retrieve original gene id
                    gene_id_par = row.select('gene_id_par').collect().item() # retrieve PAR gene id
                    new_line = line.replace(gene_id,gene_id_par) # replace the gene id with PAR gene id
                else:
                    new_line = line
                # write the line
                par.write(new_line)

    return


def main():
    parser = argparse.ArgumentParser(description="gxfpar")
    parser.add_argument("gxf_files",nargs="+",help="gtf and gff file path(s)")
    args = parser.parse_args()

    # files 
    n_files = len(args.gxf_files)
    for i,file in enumerate(args.gxf_files,start=1):
        print(f"{i}/{n_files}\tFile:{file}", end="\t")
        add_par(file)
        print("DONE...")
        


if __name__ == "__main__":
    main()
