# -*- coding: utf-8 -*-

"""Database model for bio2bel_phosphosite"""

from sqlalchemy import Column, ForeignKey, Index, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from pybel.constants import HAS_VARIANT
from pybel.dsl import pmod, protein, protein_substitution
from pybel.language import amino_acid_dict

__all__ = [
    'Base',
    'Species',
    'Protein',
    'Modification',
    'MutationEffect',
    'Mutation',
    'ModificationType',
]

Base = declarative_base()

TABLE_PREFIX = 'phosphosite'
SPECIES_TABLE_NAME = f'{TABLE_PREFIX}_species'
PROTEIN_TABLE_NAME = f'{TABLE_PREFIX}_protein'
MODIFICATION_TYPE_TABLE_NAME = f'{TABLE_PREFIX}_modificationType'
MODIFICATION_TABLE_NAME = f'{TABLE_PREFIX}_modification'
MUTATION_TABLE_NAME = f'{TABLE_PREFIX}_mutation'
MUTATION_MODIFICATION_TABLE_NAME = f'{TABLE_PREFIX}_mutation_modification'


class Species(Base):
    __tablename__ = SPECIES_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255), unique=True, index=True, nullable=False)

    def __repr__(self):
        return self.name


class Protein(Base):
    __tablename__ = PROTEIN_TABLE_NAME

    id = Column(Integer, primary_key=True)

    gene_name = Column(String(255))
    protein_name = Column(String(255))
    uniprot_id = Column(String(255), unique=True, index=True, nullable=False)

    species_id = Column(Integer, ForeignKey(f'{SPECIES_TABLE_NAME}.id'), nullable=True)
    species = relationship(Species)

    def __repr__(self):
        if self.gene_name:
            return f'{self.uniprot_id} ({self.gene_name})'
        return self.uniprot_id

    def as_bel(self):
        """Returns this model as a BEL entity"""
        return protein(
            namespace='UNIPROT',
            name=str(self.uniprot_id)
        )


class ModificationType(Base):
    """Represents the type of modifications"""
    __tablename__ = MODIFICATION_TYPE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255), unique=True, nullable=False, index=True)

    def __repr__(self):
        return self.name


class Modification(Base):
    """Represents a protein with a post-translational modification"""
    __tablename__ = MODIFICATION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    protein_id = Column(Integer, ForeignKey(f'{PROTEIN_TABLE_NAME}.id'), nullable=False)
    protein = relationship(Protein)

    residue = Column(String(3), doc='Amino acid residue name')
    position = Column(Integer, doc='Position in protein')

    modification_type_id = Column(Integer, ForeignKey(f'{MODIFICATION_TYPE_TABLE_NAME}.id'), nullable=False)
    modification_type = relationship(ModificationType)

    @property
    def type(self):
        return self.modification_type.name

    def as_bel(self) -> protein:
        parent = self.protein.as_bel()
        variant = pmod(
            name=self.modification_type.name,
            position=self.position,
            code=amino_acid_dict[self.residue.upper()],
        )
        return parent.with_variants(variant)

    def add_as_relation(self, graph):
        """Adds this modification to the graph

        :param pybel.BELGraph graph:
        :return: The edge hash as a string
        :rtype: str
        """
        return graph.add_unqualified_edge(self.protein.as_bel(), self.as_bel(), HAS_VARIANT)

    Index('idx_mod', 'type', 'protein_id', 'residue', 'position')

    def __repr__(self):
        return f'{self.protein} {self.residue}{self.position} {self.modification_type}'


class Mutation(Base):
    """Keeps track of proteins that have an amino acid substitution"""
    __tablename__ = MUTATION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    protein_id = Column(Integer, ForeignKey(f'{PROTEIN_TABLE_NAME}.id'), nullable=False)
    protein = relationship(Protein)

    dbsnp = Column(String(255), nullable=True, doc='The dbSNP RS### identifier of the mutation')
    var_type = Column(String(32))  # either 'Unclassified', 'Polymorphism', 'Disease', 'missense'

    from_aa = Column(String(1))
    position = Column(Integer)
    to_aa = Column(String(1))

    def as_bel(self) -> protein:
        """Return this mutated protein"""
        parent = self.protein.as_bel()
        modification = protein_substitution(self.from_aa, self.position, self.to_aa)
        return parent.with_variants(modification)

    def add_as_relation(self, graph):
        """Adds this modification to the graph

        :param pybel.BELGraph graph:
        :return: The edge hash as a string
        :rtype: str
        """
        return graph.add_unqualified_edge(self.protein.as_bel(), self.as_bel(), HAS_VARIANT)

    def __repr__(self):
        return f'{self.protein} {self.from_aa}{self.position}{self.to_aa}'


class MutationEffect(Base):
    __tablename__ = MUTATION_MODIFICATION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    mutation_id = Column(Integer, ForeignKey(f'{MUTATION_TABLE_NAME}.id'), nullable=False)
    mutation = relationship(Mutation)

    modification_id = Column(Integer, ForeignKey(f'{MODIFICATION_TABLE_NAME}.id'), nullable=False)
    modification = relationship(Modification)

    var_position = Column(Integer, doc='Distance of mutation to modification position')
