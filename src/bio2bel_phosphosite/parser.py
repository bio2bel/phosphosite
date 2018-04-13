# -*- coding: utf-8 -*-

import pandas as pd

from bio2bel import make_downloader
from bio2bel_phosphosite.constants import (
    O_GALNAC_PATH, O_GALNAC_URL, PHOSPHORYLATION_PATH, PHOSPHORYLATION_URL, SUMOYLATION_PATH,
    SUMOYLATION_URL, UBIQUITINATION_PATH, UBIQUITINATION_URL,
)

__all__ = [
    'get_phosphorylation_df',
    'get_o_galnac_df',
    'get_sumoylation_df',
    'get_ubiquinitation_df',
]


def make_modification_df_getter(data_url, data_path):
    download_function = make_downloader(data_url, data_path)

    def get_modifications_df(url=None, cache=True, force_download=False):
        """Gets the modifications site flat file

        :param Optional[str] url: The URL (or file path) to download. Defaults to the ChEBI data.
        :param bool cache: If true, the data is downloaded to the file system, else it is loaded from the internet
        :param bool force_download: If true, overwrites a previously cached file
        :rtype: pandas.DataFrame
        """
        if url is None and cache:
            url = download_function(force_download=force_download)

        return pd.read_csv(
            url or data_url,
            skiprows=2,
            sep='\t'
        )

    return get_modifications_df


get_phosphorylation_df = make_modification_df_getter(PHOSPHORYLATION_URL, PHOSPHORYLATION_PATH)
get_ubiquinitation_df = make_modification_df_getter(UBIQUITINATION_URL, UBIQUITINATION_PATH)
get_o_galnac_df = make_modification_df_getter(O_GALNAC_URL, O_GALNAC_PATH)
get_sumoylation_df = make_modification_df_getter(SUMOYLATION_URL, SUMOYLATION_PATH)
