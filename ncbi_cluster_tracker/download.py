import os
import re
import requests
import shutil

import pandas as pd
import tqdm

from ncbi_cluster_tracker.logger import logger

# Can't auto-translate these to folders on NCBI FTP
TAXGROUP_TO_ORGANISM = {
    'Acinetobacter baumannii': 'Acinetobacter',
    'Aeromonas hydrophila': 'Aeromonas',
    'Campylobacter jejuni': 'Campylobacter',
    'E.coli and Shigella': 'Escherichia_coli_Shigella',
    'Elizabethkingia anophelis': 'Elizabethkingia',
    'Klebsiella pneumoniae': 'Klebsiella',
    'Listeria monocytogenes': 'Listeria',
    'Morganella morganii': 'Morganella',
    'Providencia alcalifaciens': 'Providencia',
    'Salmonella enterica': 'Salmonella',
    'Serratia marcescens': 'Serratia',
}


def download_cluster_files(clusters_df: pd.DataFrame):
    clusters_df['snp_url'] = clusters_df.apply(build_snp_url, axis=1)
    urls = clusters_df['snp_url'].to_list()
    download_snps(urls)


def build_ftp_base_url(taxgroup_name: str) -> str:
    url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Results/'
    url += get_organism_from_taxgroup(taxgroup_name)
    url += '/latest_snps/'
    return url


def get_organism_from_taxgroup(taxgroup_name: str) -> str:
    if taxgroup_name in TAXGROUP_TO_ORGANISM:
        organism = f'{TAXGROUP_TO_ORGANISM[taxgroup_name]}'
    else:
        organism = f'{taxgroup_name.replace(' ', '_')}'
    return organism


def build_snp_url(row: pd.Series) -> str:
    url = build_ftp_base_url(row['taxgroup_name'])
    url += f'SNP_trees/{row['cluster']}.tar.gz'
    return url


def get_latest_kmer_group_acc(taxgroup_name: str) -> str:
    # BigQuery doesn't give us the specific kmer_group_acc, so we'll have to 
    # pull it from the XML descriptor file on the FTP
    url = build_ftp_base_url(taxgroup_name)
    response = requests.get(url)
    if not response.ok:
        raise Exception('Could not find k-mer group accession: {url}')
    acc_match = re.search(r'PDG(\d+)\.(\d+)', response.text)
    if acc_match is not None:
        acc = acc_match.group(0) 
    else:
        raise Exception(f'Could not find k-mer group accession: {url}')
    return acc


def build_tree_viewer_url(
        row: pd.Series,
        taxgroup_to_kmer_group_acc: dict[str, str] = {}
) -> str:
    taxgroup_name = row['taxgroup_name']
    if row['cluster'] not in taxgroup_to_kmer_group_acc:
        taxgroup_to_kmer_group_acc[taxgroup_name] = (
            get_latest_kmer_group_acc(taxgroup_name)
        )
    url = 'https://www.ncbi.nlm.nih.gov/pathogens/tree#'
    url += get_organism_from_taxgroup(taxgroup_name)
    url += '/'
    url += taxgroup_to_kmer_group_acc[taxgroup_name]
    url += '/'
    url += row['cluster']
    return url 


def download_snps(urls: list[str]) -> None:
    out_subdir = os.path.join(os.environ['NCT_OUT_SUBDIR'], 'snps')
    logger.info(f'Downloading SNP cluster data to {out_subdir}...')
    os.makedirs(out_subdir, exist_ok=True)
    for url in tqdm.tqdm(urls):
        with requests.Session() as session:
            for _ in range(3):
                response = session.get(url)
                if response.ok:
                    break
                # try incrementing cluster version
                regex = r'(?<=\.)(\d+)(?=\.tar\.gz)'
                url = re.sub(regex, lambda x: str(int(x.group())+1), url)
            else:
                raise Exception(f'Could not find cluster: {url}')

            destination = os.path.join(out_subdir, os.path.basename(url))
            with open(destination, 'wb') as f:
                f.write(response.content)
            shutil.unpack_archive(destination, out_subdir)
