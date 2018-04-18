# -*- coding: utf-8 -*-

"""

Format
------
The PTMVar excel document has two sheets. On the first is the legend, which is summarized below. The columns are split
into two groups, corresponding to the source of the data being either UniProt or PhosphoSitePlus.

Source: UniProt KB

COLUMN	HEADER	DESCRIPTION
A	GENE	Unique ID (UID) assigned to a gene by the HUGO Nomenclature Committee and adapted by HUPO for the cognate protein.
B	UPID	UniProtKB ID
C	CHR_LOC	The chromosomal location of the gene that encodes the mutated protein.
D	VAR	The UniProtKB/Swiss-Prot UID of the variant.
E	dbSNP	The UID of the variant from dbSNP, the NCBI database of Single Nucleotide Polymorphisms (www.ncbi.nlm.nih.gov/SNP/).
F	AA_CHANGE	The amino acid substitution caused by the nsSNP. For example, Ser16Pro (or S16P) indicates a residue change from serine (WT) to proline at residue number 16.
G	WT	The wild type (WT) amino acid that is changed by the mutation
H	RSD#	The sequence number of the residue that is mutated
I	VAR	The variant amino acid encoded by the missense mutation
J	VAR_TYPE	Variant type: Polymorphism, Disease or Unclassified
K	DISEASE(S)	Specific disease(s) associated with the mutation. Includes the abbreviation and UID from the Online Mendelian Inheritance in Man® (OMIM®). For example, the abbreviation and UID for phenylketonuria are PKU and MIM:261600.
L	MUT_SOURCE	"The sources of mutant data include HUMSAVAR (www.uniprot.org/docs/humsavar), The Cancer Genome Atlas (TCGA; cancergenome.nih.gov/), Catalogue Of Somatic Mutations In Cancer (COSMIC; cancer.sanger.ac.uk/cancergenome/), and the cBioPortal for Cancer Genomics (cBIO; www.cbioportal.org/).

Source: PhosphoSitePlus® (PSP)	M	PROTEIN	Primary name in PhosphoSitePlus® (PSP)

COLUMN	HEADER	DESCRIPTION
N	ACC_ID	Primary accession ID in PSP
O	MOD_RSD	The sequence number of the residue with the posttranslational modification (PTM)
P	MOD_AA	The type of amino acid with the PTM
Q	CONSERVATION	The conservation of the modified amino acid between species: human (H), mouse (M), rat )R)
R	MOD_TYPE	The type of PTM occuring on the modified residue, e.g.phosphorylation, ubiquitylation, etc.
S	SITE_GRP_ID	Unique identifier of a modification site and its homologous sites in all proteoforms and species
T	MOD-SITE_SEQ	WT sequence: includes the modsite (lower case letter) plus surrounding sequence (+/- 5 rsds.). The location of the mutated residue is marked with an asterisk.
U	VAR_POSITION	The mutant position relative to posttranslational modification. N-terminal residues are marked with a minus sign. For example, '-3' indicates that the mutated residue lies three residues N-terminal to the modsite. A zero indicates that the modified residue itself has been mutated, i.e. is a Class I Variant..
V	VAR_SITE_SEQ	Variant sequence: includes the modsite (lower case letter) plus surrounding sequence (+/- 5 rsds.). The mutated residue is marked with an asterisk.
W	VAR_CLASS*	CLASS I (site loss), CLASS Ia (modsite switch), or CLASS II (flanking change) *
X	LTP_LIT	The number of associated literature records derived using low-throughput experimental techniques. An LTP result may be more reliable than an MS2 result.
Y	MS2_LIT	The number of published articles using proteomic MS experiments to locate residues that are modified.
Z	CST_CS	The number of associated shotgun proteomic experiments performed at Cell Signaling Technology (CST) in which the indicated PTM was observed.


Extra information:

* VAR_CLASS
CLASS I (site loss)	[STY]→{STY} OR [KR]→{KR}
CLASS Ia (modsite switch)	[Y]→[ST] OR [TS]→[Y]
CLASS II (flanking change)	variant +/- 5 AAs from the modification site

"""

import zipfile

import pandas as pd

from bio2bel import make_downloader
from ..constants import PTMVAR_PATH, PTMVAR_URL

__all__ = [
    'download_ptmvar',
    'get_ptmvar_df',
]

_map_mod_types = {
    'Ubiquitylation',
    'Phosphorylation',
    'Acetylation',
    'Sumoylation',
    'Succinylation',
    'Neddylation',
    'Methylation'
}

download_ptmvar = make_downloader(PTMVAR_URL, PTMVAR_PATH)


def get_ptmvar_df(url=None, cache=True, force_download=False):
    """Gets the PTMVar excel sheet

    :param Optional[str] url: The URL (or file path) to download.
    :param bool cache: If true, the data is downloaded to the file system, else it is loaded from the internet
    :param bool force_download: If true, overwrites a previously cached file
    :rtype: pandas.DataFrame
    """
    if url is None and cache:
        url = download_ptmvar(force_download=force_download)

    with zipfile.ZipFile(url) as zf:
        with zf.open('PTMVar.xlsx') as f:
            rv = pd.read_excel(
                f,
                sheet_name=1,
                skiprows=6,
            )

    # remove weird forward quote
    rv['VAR_POSITION'] = rv['VAR_POSITION'].apply(lambda x: x.strip("'"))

    return rv
