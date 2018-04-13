# -*- coding: utf-8 -*-

import os

from bio2bel.utils import get_connection, get_data_dir

MODULE_NAME = 'phosphosite'
DATA_DIR = get_data_dir(MODULE_NAME)
DEFAULT_CACHE_CONNECTION = get_connection(MODULE_NAME)

PHOSPHORYLATION_URL = 'https://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz'
PHOSPHORYLATION_PATH = os.path.join(DATA_DIR, 'Phosphorylation_site_dataset.gz')

O_GALNAC_URL = 'https://www.phosphosite.org/downloads/O-GalNAc_site_dataset.gz'
O_GALNAC_PATH = os.path.join(DATA_DIR, 'O-GalNAc_site_dataset.gz')

SUMOYLATION_URL = 'https://www.phosphosite.org/downloads/Sumoylation_site_dataset.gz'
SUMOYLATION_PATH = os.path.join(DATA_DIR, 'Sumoylation_site_dataset.gz')

UBIQUITINATION_URL = 'https://www.phosphosite.org/downloads/Ubiquitination_site_dataset.gz'
UBIQUITINATION_PATH = os.path.join(DATA_DIR, 'Ubiquitination_site_dataset.gz')
