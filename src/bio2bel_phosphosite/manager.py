# -*- coding: utf-8 -*-

import logging

from bio2bel import AbstractManager
from .constants import MODULE_NAME
from .models import Base, Modification, Protein, Species
from .parsers.phosphorylation import get_phosphorylation_df

__all__ = ['Manager']

log = logging.getLogger(__name__)


class Manager(AbstractManager):
    module_name = MODULE_NAME
    flask_admin_models = [Protein, Species, Modification]

    @property
    def base(self):
        return Base

    def populate(self, url=None):
        df = get_phosphorylation_df(url=url)
