from pathlib import Path

Path('rawdata/hg19').mkdir(parents=True, exist_ok=True)


CHR_NAMES = [f'chr{i}' for i in list(range(1, 22)) + ['X', 'Y']]
GP_URL = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/{}.fa.gz'



rule all:
    input:
        expand('output/{chrom_name}.hybrid_init.bedgraph', chrom_name=CHR_NAMES)

rule download_hg19_chromosomes:
    input:
        expand('rawdata/hg19/{chr_name}.fa.gz',
        chr_name=CHR_NAMES)


rule download_chromosome_file:
    output:
        'rawdata/hg19/{chr_name}.fa.gz'
    params:
        url = lambda wildcards: GP_URL.format(wildcards['chr_name'])
    shell:'''
    wget {params.url} -O {output}
    '''


rule unzip_chrom_file:
    input:
        'rawdata/hg19/{chr_name}.fa.gz'
    output:
        'rawdata/hg19/{chr_name}.fa'
    shell:'''
    gzip -d {input}
    '''


rule profile_hybrid_init_energies:
    conda:
        'envs/python.yml'
    input:
        'rawdata/hg19/{chrom_name}.fa'
    output:
        'output/{chrom_name}.hybrid_init.bedgraph'
    params:
        chr_name = lambda wildcards: wildcards['chrom_name'],
        window_size=300,
        dist_between_windows=50,  # effectively the initiation region
    script:"scripts/profile_hybrid_energies.py"



