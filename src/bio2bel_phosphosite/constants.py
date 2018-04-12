# -*- coding: utf-8 -*-

import os

from bio2bel.utils import get_connection, get_data_dir

MODULE_NAME = 'phosphosite'
DATA_DIR = get_data_dir(MODULE_NAME)
DEFAULT_CACHE_CONNECTION = get_connection(MODULE_NAME)

PHOSPHORYLATION_URL = 'https://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz'
PHOSPHORYLATION_PATH = os.path.join(DATA_DIR, 'Phosphorylation_site_dataset.gz')
