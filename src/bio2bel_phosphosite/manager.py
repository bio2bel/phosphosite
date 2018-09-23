# -*- coding: utf-8 -*-

import logging
from typing import List, Mapping, Optional

import time
from sqlalchemy import and_, func
from tqdm import tqdm

from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from pybel import BELGraph
from .constants import MODULE_NAME, PROTEIN_NAMESPACE
from .models import Base, Modification, ModificationType, Mutation, MutationEffect, Protein, Species
from .parsers import (
    get_acetylation_df, get_o_galnac_df, get_o_glcnac_df, get_phosphorylation_df, get_ptmvar_df, get_sumoylation_df,
    get_ubiquinitation_df,
)

__all__ = ['Manager']

log = logging.getLogger(__name__)

_pmod_map = {
    'ac': 'Ac',
    'p': 'Ph',
    'ga': 'OGlyco',  # TODO make more specific!
    'gl': 'OGlyco',  # TODO make more specific!
    'sm': 'Sumo',
    'ub': 'Ub',
    'Acetylation': 'Ac',
    'Ubiquitylation': 'Ub',
    'Phosphorylation': 'Ph',
    'Sumoylation': 'Sumo',
    'Succinylation': 'Succ',
    'Neddylation': 'Nedd',
    'Methylation': 'Me',
}

_ptmvar_rows = ['UPID', 'ACC_ID', 'dbSNP', 'WT_AA', 'MUT_RSD#', 'VAR_AA', 'VAR_TYPE', 'MOD_RSD', 'MOD_AA', 'MOD_TYPE',
                'VAR_POSITION']


def _parse_mod(s):
    """Parses the modification string. Follows the format Letter + Integer + Dash + Code

    :param s:
    """
    residue = s[0]
    position, modification_type = s[1:].split('-')

    return residue, int(position), _pmod_map[modification_type]


class Manager(AbstractManager, BELManagerMixin, FlaskMixin):
    """Manager for PhosphoSitePlus."""

    module_name = MODULE_NAME
    flask_admin_models = [Protein, Modification, Mutation, MutationEffect, ModificationType, Species]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.name_to_modification_type = {}
        self.name_to_species = {}
        self.uniprot_id_to_protein = {}
        self.modifications = {}
        self.mutations = {}

    @property
    def _base(self):
        return Base

    def is_populated(self) -> bool:
        """Check if the database is already populated."""
        return 0 < self.count_modifications()

    def get_modification_type_by_name(self, name) -> Optional[ModificationType]:
        return self.session.query(ModificationType).filter(ModificationType.name == name).one_or_none()

    def get_or_create_modification_type(self, name) -> ModificationType:
        modification_type = self.name_to_modification_type.get(name)
        if modification_type is not None:
            return modification_type

        modification_type = self.get_modification_type_by_name(name)
        if modification_type is not None:
            self.name_to_modification_type[name] = modification_type
            return modification_type

        modification_type = self.name_to_modification_type[name] = ModificationType(name=name)
        self.session.add(modification_type)
        return modification_type

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

    def get_mutation(self, uniprot_id: str, from_aa: str, position: int, to_aa: str) -> Optional[Mutation]:
        return self.session.query(Mutation).join(Protein).filter(and_(
            Protein.uniprot_id == uniprot_id,
            Mutation.from_aa == from_aa,
            Mutation.to_aa == to_aa,
            Mutation.position == position
        )).one_or_none()

    def get_or_create_mutation(self, uniprot_id: str, from_aa: str, position: int, to_aa: str) -> Mutation:
        _tuple = (uniprot_id, from_aa, position, to_aa)
        mutation = self.mutations.get(_tuple)
        if mutation is not None:
            return mutation

        mutation = self.get_mutation(uniprot_id, from_aa, position, to_aa)
        if mutation is not None:
            self.mutations[_tuple] = mutation
            return mutation

        mutation = self.modifications[_tuple] = Mutation(
            protein=self.get_or_create_protein(uniprot_id),
            from_aa=from_aa,
            position=position,
            to_aa=to_aa,
        )
        self.session.add(mutation)
        return mutation

    def get_modification(self, uniprot_id: str, residue: str, position: int, modification_type: str) -> Optional[
                         Modification]:
        return self.session.query(Modification).join(Protein).join(ModificationType).filter(and_(
            Protein.uniprot_id == uniprot_id,
            Modification.residue == residue,
            Modification.position == position,
            ModificationType.name == modification_type
        )).one_or_none()

    def get_or_create_modification(self, uniprot_id: str, residue: str, position: int,
                                   modification_type: str) -> Modification:
        """
        :param uniprot_id: The UniProt identifier
        :param residue: Amino acid
        :param position: Position
        :param modification_type: Name of modification
        :return:
        """
        _tuple = (uniprot_id, residue, position, modification_type)
        modification = self.modifications.get(_tuple)
        if modification is not None:
            return modification

        modification = self.get_modification(uniprot_id, residue, position, modification_type)
        if modification is not None:
            self.modifications[_tuple] = modification
            return modification

        modification = self.modifications[_tuple] = Modification(
            protein=self.get_or_create_protein(uniprot_id),
            residue=residue,
            position=position,
            modification_type=self.get_or_create_modification_type(modification_type),
        )
        self.session.add(modification)
        return modification

    def _populate_modification_df(self, df):
        log.info('building models')
        for organism_name, organism_df in tqdm(df.groupby('ORGANISM'), desc='Species'):

            species = self.get_or_create_species(organism_name)

            organism_group = organism_df.groupby(['GENE', 'PROTEIN', 'ACC_ID'])
            organism_group_it = tqdm(organism_group, total=len(organism_group), desc=organism_name, leave=False)
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
                        modification_type=self.get_or_create_modification_type(modification_type),
                    )
                    self.session.add(modification)

        t = time.time()
        log.info('committing models')
        self.session.commit()
        log.info('done committing models in %.2f seconds', time.time() - t)

    def count_residues(self) -> Mapping[str, int]:
        """Count the frequency of modification on each residue type."""
        return dict(
            self.session
                .query(Modification.residue, func.count(Modification.residue))
                .group_by(Modification.residue)
                .all()
        )

    def count_modification_types(self) -> Mapping[str, int]:
        """Count the number of each modification type."""
        return dict(
            self.session
                .query(ModificationType.name, func.count(ModificationType.name))
                .join(Modification)
                .group_by(ModificationType.name)
                .all()
        )

    def count_proteins(self) -> int:
        return self._count_model(Protein)

    def count_species(self) -> int:
        return self._count_model(Species)

    def list_species(self) -> List[Species]:
        return self._list_model(Species)

    def count_modifications(self) -> int:
        return self._count_model(Modification)

    def count_mutations(self) -> int:
        return self._count_model(Mutation)

    def count_mutation_effects(self) -> int:
        return self._count_model(MutationEffect)

    def summarize(self):
        return dict(
            proteins=self.count_proteins(),
            species=self.count_species(),
            modifications=self.count_modifications(),
            residues=self.count_residues(),
            modification_types=self.count_modification_types(),
            mutations=self.count_mutations(),
            mutation_effects=self.count_mutation_effects(),

        )

    def list_modifications(self) -> List[Modification]:
        """List all modifications."""
        return self.session.query(Modification).all()

    def list_mutation_effects(self) -> List[MutationEffect]:
        """List all mutation effects."""
        return self.session.query(MutationEffect).all()

    def _populate_modifications(self, phosphorylation_url=None, sumoylation_url=None, ubiquitination_url=None,
                                o_galnac_url=None, o_glcnac_url=None, acetylation_url=None):
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

    def _populate_ptmvar(self, url: Optional[str] = None) -> None:
        """Download and populate the PTMVar data set."""
        df = get_ptmvar_df(url=url)

        it = tqdm(df[_ptmvar_rows].itertuples(), total=len(df.index), desc='PTMVar')

        for idx, upid, upid2, dbsnp, from_aa, mut_rsd, to_aa, var_type, mod_rsd, mod_aa, mod_type, var_position in it:

            if upid != upid2:
                log.warning('problem with line - non-matching identifiers %s and %s', upid, upid2)
                continue

            mutation = self.get_or_create_mutation(upid, from_aa, mut_rsd, to_aa, var_type=var_type, dbsnp=dbsnp)
            modification = self.get_or_create_modification(upid, residue=mod_aa, position=mod_rsd,
                                                           modification_type=_pmod_map[mod_type])

            e = MutationEffect(
                mutation=mutation,
                modification=modification,
                var_position=var_position
            )
            self.session.add(e)

        t = time.time()
        log.info('committing models')
        self.session.commit()
        log.info('done committing models in %.2f seconds', time.time() - t)

    def populate(self, phosphorylation_url=None, sumoylation_url=None, ubiquitination_url=None, o_galnac_url=None,
                 o_glcnac_url=None, acetylation_url=None, ptmvar_url=None):
        """Downloads and populates data

        :param phosphorylation_url:
        :param sumoylation_url:
        :param ubiquitination_url:
        :param o_galnac_url:
        :param o_glcnac_url:
        :param acetylation_url:
        :param ptmvar_url:
        """
        self._populate_modifications(
            phosphorylation_url=phosphorylation_url,
            sumoylation_url=sumoylation_url,
            ubiquitination_url=ubiquitination_url,
            o_galnac_url=o_galnac_url,
            o_glcnac_url=o_glcnac_url,
            acetylation_url=acetylation_url,
        )

        self._populate_ptmvar(url=ptmvar_url)

    def to_bel(self) -> BELGraph:
        """Converts PhosphoSite knowledge to BEL"""
        graph = BELGraph(
            name='PhosphositePlus Modifications',
            version='1.0.0'  # need to get from data source itself
        )

        graph.namespace_url[
            PROTEIN_NAMESPACE] = 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/uniprot/uniprot-20170813.belns'

        for m in tqdm(self.list_modifications(), total=self.count_modifications(), desc='modifications'):
            m.add_as_relation(graph)

        for me in tqdm(self.list_mutation_effects(), total=self.count_mutation_effects(), desc='mutation effects'):
            me.add_as_relation(graph)

        return graph
