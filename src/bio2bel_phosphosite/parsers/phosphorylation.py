# -*- coding: utf-8 -*-

import pandas as pd

from bio2bel import make_downloader
from ..constants import PHOSPHORYLATION_PATH, PHOSPHORYLATION_URL

__all__ = [
    'download_phosphorylation',
    'get_phosphorylation_df',
]

download_phosphorylation = make_downloader(PHOSPHORYLATION_URL, PHOSPHORYLATION_PATH)


def get_phosphorylation_df(url=None, cache=True, force_download=False):
    """Gets the phosphorylations site flat file

    :param Optional[str] url: The URL (or file path) to download. Defaults to the ChEBI data.
    :param bool cache: If true, the data is downloaded to the file system, else it is loaded from the internet
    :param bool force_download: If true, overwrites a previously cached file
    :rtype: pandas.DataFrame
    """
    if url is None and cache:
        url = download_phosphorylation(force_download=force_download)

    return pd.read_csv(
        url or PHOSPHORYLATION_URL,
        skiprows=2,
        sep='\t'
    )
