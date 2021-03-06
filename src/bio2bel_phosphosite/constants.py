# -*- coding: utf-8 -*-

import os

from bio2bel.utils import get_connection, get_data_dir

MODULE_NAME = 'phosphosite'
DATA_DIR = get_data_dir(MODULE_NAME)
DEFAULT_CACHE_CONNECTION = get_connection(MODULE_NAME)

PROTEIN_NAMESPACE = 'UNIPROT'

PHOSPHORYLATION_URL = 'https://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz'
PHOSPHORYLATION_PATH = os.path.join(DATA_DIR, 'Phosphorylation_site_dataset.gz')

ACETYLATION_URL = 'https://www.phosphosite.org/downloads/Acetylation_site_dataset.gz'
ACETYLATION_PATH = os.path.join(DATA_DIR, 'Acetylation_site_dataset.gz')

O_GALNAC_URL = 'https://www.phosphosite.org/downloads/O-GalNAc_site_dataset.gz'
O_GALNAC_PATH = os.path.join(DATA_DIR, 'O-GalNAc_site_dataset.gz')

O_GLCNAC_URL = 'https://www.phosphosite.org/downloads/O-GlcNAc_site_dataset.gz'
O_GLCNAC_PATH = os.path.join(DATA_DIR, 'O-GlcNAc_site_dataset.gz')

SUMOYLATION_URL = 'https://www.phosphosite.org/downloads/Sumoylation_site_dataset.gz'
SUMOYLATION_PATH = os.path.join(DATA_DIR, 'Sumoylation_site_dataset.gz')

UBIQUITINATION_URL = 'https://www.phosphosite.org/downloads/Ubiquitination_site_dataset.gz'
UBIQUITINATION_PATH = os.path.join(DATA_DIR, 'Ubiquitination_site_dataset.gz')

REGULATORY_SITES_URL = 'https://www.phosphosite.org/downloads/Regulatory_sites.gz'
REGULATORY_SITES_PATH = os.path.join(DATA_DIR, 'Regulatory_sites.gz')

DISEASE_ASSOCIATED_SITES_URL = 'https://www.phosphosite.org/downloads/Disease-associated_sites.gz'
DISEASE_ASSOCIATED_SITES_PATH = os.path.join(DATA_DIR, 'Disease-associated_sites.gz')

PTMVAR_URL = 'https://www.phosphosite.org/downloads/PTMVar.xlsx.zip'
PTMVAR_PATH = os.path.join(DATA_DIR, 'PTMVar.xlsx.zip')
