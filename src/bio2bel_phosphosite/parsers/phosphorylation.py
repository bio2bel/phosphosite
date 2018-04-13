# -*- coding: utf-8 -*-

from bio2bel_phosphosite.parser import make_modification_df_getter
from ..constants import PHOSPHORYLATION_PATH, PHOSPHORYLATION_URL

__all__ = [
    'get_phosphorylation_df',
]

get_phosphorylation_df = make_modification_df_getter(PHOSPHORYLATION_URL, PHOSPHORYLATION_PATH)
