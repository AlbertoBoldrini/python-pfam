# Pfam interface for Python
This package provides a programming interface to retrive information from [Pfam database](https://pfam.xfam.org/).

The [Pfam database](https://pfam.xfam.org/) is a large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs).

## Installation
This package is not in the `pip` repository yet. To install it you can use the following commands:
```sh
$ git clone https://github.com/AlbertoBoldrini/python-pfam
$ pip install ./python-pfam
```


## Example
In these example we use the package `beeprint` to print to the termial data structured in a readable format. 

```python
import pfam
from beeprint import pp

# Retrive informations about the Piwi family from Pfam
piwi = pfam.family('Piwi')

# This data structure contains many fields
print('Description of Piwi family: ', piwi.comment)
print('Piwi family has %d different sequences.' % piwi.curation_details.num_seqs_full)

# We use pp() to print all data contained in the piwi structure
pp(piwi)

# Retrive the list of proteins in this family (big file...)
protein_list = piwi.proteins()

# Print a very very long list
# print(protein_list)

# Print the id of the first protein in this family
print('The ID of the first protein of the Piwi family is: ', protein_list[0].id)

# Retrive information about the first protein in this family
first_protein = protein_list[0].fetch()

# Print the sequence on this protein
print('The sequence of the first protein of the Piwi family is: ', first_protein.sequence)

# Retrive the clan (superfamily) of the Piwi family
clan = piwi.clan()

# What is the name of this clan?
print ("Piwi family belongs to %s clan." % clan.entry.id)

# Print the list of families in this clan
pp(clan.families)
```