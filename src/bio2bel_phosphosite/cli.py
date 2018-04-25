# -*- coding: utf-8 -*-

"""Run this script with :code:`python3 -m bio2bel_phosphosite`"""

import click

from .manager import Manager
from .models import Modification

main = Manager.get_cli()


@main.group()
def manage():
    pass


@manage.group()
def protein():
    pass


@protein.command()
@click.argument('uniprot_id')
@click.pass_obj
def get(manager, uniprot_id):
    """Get a summary of a protein by UniProt identifier"""
    p = manager.get_protein_by_uniprot_id(uniprot_id)

    if p is None:
        click.echo(f'could not find {uniprot_id}')

    unique_positions = {m.position for m in p.modifications.all()}
    click.echo(f'Unique positions modified: {len(unique_positions)}')

    for m in p.modifications.order_by(Modification.position):
        click.echo(f'{m.position} {m.residue} {m.modification_type}')


if __name__ == '__main__':
    main()
