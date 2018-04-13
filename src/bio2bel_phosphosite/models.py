# -*- coding: utf-8 -*-

"""Database model for bio2bel_phosphosite"""

from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

__all__ = [
    'Base',
    'Species',
    'Protein',
    'Modification',
]

Base = declarative_base()

TABLE_PREFIX = 'phosphosite'
SPECIES_TABLE_NAME = f'{TABLE_PREFIX}_species'
PROTEIN_TABLE_NAME = f'{TABLE_PREFIX}_protein'
MODIFICATION_TABLE_NAME = f'{TABLE_PREFIX}_modification'


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

    species_id = Column(Integer, ForeignKey(f'{SPECIES_TABLE_NAME}.id'), nullable=False)
    species = relationship(Species)

    def __repr__(self):
        return self.gene_name


class Modification(Base):
    """Represents a chemical"""
    __tablename__ = MODIFICATION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    protein_id = Column(Integer, ForeignKey(f'{PROTEIN_TABLE_NAME}.id'), nullable=False)
    protein = relationship(Protein)

    residue = Column(String(3))
    position = Column(Integer)
    type = Column(String(255))
