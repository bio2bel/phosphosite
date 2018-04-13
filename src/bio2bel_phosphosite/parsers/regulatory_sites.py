# -*- coding: utf-8 -*-

import pandas as pd

from bio2bel import make_downloader
from ..constants import REGULATORY_SITES_PATH, REGULATORY_SITES_URL

__all__ = [
    'download_regulatory_sites'
]

download_regulatory_sites = make_downloader(REGULATORY_SITES_URL, REGULATORY_SITES_PATH)


def get_regulatory_sites_df(url=None, cache=True, force_download=False):
    """Gets the modifications site flat file

    :param Optional[str] url: The URL (or file path) to download.
    :param bool cache: If true, the data is downloaded to the file system, else it is loaded from the internet
    :param bool force_download: If true, overwrites a previously cached file
    :rtype: pandas.DataFrame
    """
    if url is None and cache:
        url = download_regulatory_sites(force_download=force_download)

    return pd.read_csv(
        url or REGULATORY_SITES_URL,
        skiprows=2,
        sep='\t'
    )
