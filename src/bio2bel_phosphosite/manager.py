# -*- coding: utf-8 -*-

import logging
import time
from typing import Optional

from sqlalchemy import func
from tqdm import tqdm

from bio2bel import AbstractManager
from .constants import MODULE_NAME
from .models import Base, Modification, Protein, Species
from .parser import (
    get_acetylation_df, get_o_galnac_df, get_o_glcnac_df, get_phosphorylation_df, get_sumoylation_df,
    get_ubiquinitation_df,
)

__all__ = ['Manager']

log = logging.getLogger(__name__)


def _parse_mod(s):
    """Parses the modification string. Follows the format Letter + Integer + Dash + Code

    :param s:
    :rtype: str
    """
    residue = s[0]
    position, modification_type = s[1:].split('-')

    return residue, int(position), modification_type


class Manager(AbstractManager):
    module_name = MODULE_NAME
    flask_admin_models = [Protein, Species, Modification]

    def __init__(self, connection=None):
        super().__init__(connection=connection)

        self.name_to_species = {}
        self.uniprot_id_to_protein = {}

    @property
    def base(self):
        return Base

    def get_species_by_name(self, name) -> Optional[Species]:
        return self.session.query(Species).filter(Species.name == name).one_or_none()

    def get_or_create_species(self, name) -> Species:
        species = self.name_to_species.get(name)
        if species is not None:
            return species

        species = self.get_species_by_name(name)
        if species is not None:
            self.name_to_species[name] = species
            return species

        species = self.name_to_species[name] = Species(name=name)
        self.session.add(species)
        return species

    def get_protein_by_uniprot_id(self, uniprot_id) -> Optional[Protein]:
        return self.session.query(Protein).filter(Protein.uniprot_id == uniprot_id).one_or_none()

    def get_or_create_protein(self, uniprot_id, **kwargs) -> Protein:
        protein = self.uniprot_id_to_protein.get(uniprot_id)
        if protein is not None:
            return protein

        protein = self.get_protein_by_uniprot_id(uniprot_id)
        if protein is not None:
            self.uniprot_id_to_protein[uniprot_id] = protein
            return protein

        protein = self.uniprot_id_to_protein[uniprot_id] = Protein(
            uniprot_id=uniprot_id,
            **kwargs
        )
        self.session.add(protein)
        return protein

    def _populate_modification_df(self, df):
        log.info('building models')
        for organism_name, organism_df in df.groupby('ORGANISM'):

            species = self.get_or_create_species(organism_name)

            organism_group = organism_df.groupby(['GENE', 'PROTEIN', 'ACC_ID'])
            organism_group_it = tqdm(organism_group, total=len(organism_group), desc=organism_name)

            for (gene_name, protein_name, uniprot_id), protein_df in organism_group_it:

                protein = self.get_or_create_protein(
                    uniprot_id,
                    gene_name=gene_name,
                    protein_name=protein_name,
                    species=species,
                )

                for mod in protein_df.MOD_RSD:
                    residue, position, modification_type = _parse_mod(mod)

                    modification = Modification(
                        protein=protein,
                        residue=residue,
                        position=position,
                        type=modification_type,
                    )

                    self.session.add(modification)

        t = time.time()
        log.info('committing models')
        self.session.commit()
        log.info('done committing models in %.2f seconds', time.time() - t)

    def count_residues(self):
        """

        :rtype: dict[str,int]
        """
        return dict(self.session.query(Modification.residue, func.count(Modification.residue)).group_by(
            Modification.residue).all())

    def count_modification_types(self):
        """

        :rtype: dict[str,int]
        """
        return dict(
            self.session.query(Modification.type, func.count(Modification.type)).group_by(Modification.type).all())

    def count_proteins(self):
        return self._count_model(Protein)

    def count_species(self):
        return self._count_model(Species)

    def count_modifications(self):
        return self._count_model(Modification)

    def summarize(self):
        return dict(
            proteins=self.count_proteins(),
            species=self.count_species(),
            modifications=self.count_modifications(),
            residues=self.count_residues(),
            modification_types=self.count_modification_types(),
        )

    def populate(self, phosphorylation_url=None, sumoylation_url=None, ubiquitination_url=None, o_galnac_url=None,
                 o_glcnac_url=None, acetylation_url=None):
        """Downloads and populates data

        :param phosphorylation_url:
        :param sumoylation_url:
        :param ubiquitination_url:
        :param o_galnac_url:
        :param o_glcnac_url:
        :param acetylation_url:
        """
        log.info('phosphorylation')
        phosphorylation_df = get_phosphorylation_df(url=phosphorylation_url)
        self._populate_modification_df(phosphorylation_df)

        log.info('acetylation')
        acetylation_df = get_acetylation_df(url=acetylation_url)
        self._populate_modification_df(acetylation_df)

        log.info('sumoylation')
        sumoylation_df = get_sumoylation_df(url=sumoylation_url)
        self._populate_modification_df(sumoylation_df)

        log.info('ubiquitination')
        ubiquination_df = get_ubiquinitation_df(url=ubiquitination_url)
        self._populate_modification_df(ubiquination_df)

        log.info('o-galnac-ation')
        o_galnac_df = get_o_galnac_df(url=o_galnac_url)
        self._populate_modification_df(o_galnac_df)

        log.info('o-glcnac-ation')
        o_glcnac_df = get_o_glcnac_df(url=o_glcnac_url)
        self._populate_modification_df(o_glcnac_df)
