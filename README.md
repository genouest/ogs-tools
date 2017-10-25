# OGS Tools

[![Build](https://travis-ci.org/abretaud/ogs-tools.svg?branch=master)](https://travis-ci.org/abretaud/ogs-tools)

Some tools to manipulate OGS (Official Gene Sets).
These tools were written for (insect) genomes hosted on [BIPAA platform](https://bipaa.genouest.org).

We tried to make them as generic as possible, hoping that maybe they could be useful for other people. Contributions are of course welcome if it doesn't fit perfectly your use case.

Some of these tools would deserve some big refactoring, but 'It works'©.

Here's how we manage annotations at BIPAA:

- Automatic annotation with Maker
- Apollo server for manual curation, following [our guidelines](https://bipaa.genouest.org/is/how-to-annotate-a-genome/)
- Each night, our [Apollo report application](https://bipaa.genouest.org/is/how-to-annotate-a-genome/) checks that each annotated gene follows our guidelines. A report is available for each user with a list of errors and warnings that needs to be fixed, and a list of valid gene.
- Regularly we generate a new OGS by merging the automatic annotation with Apollo annotated genes (the valid ones only). We use the `ogs_merge` script to do that.
- To submit the OGS to ENA @ EBI in EMBL format, we (would like to) use `gff2embl`

## Merging OGS

This script allows to merge a base annotation with an apollo annotation.
Genes from the base annotation are replaced by overlapping apollo genes, preserving original IDs with a version suffix.

It depends on 4 software that can be installed from conda like this:

```
conda create --name ogsmerger bedops bedtools cufflinks bcbiogff
```

Then you can source the conda environment:

```
source activate ogsmerger
```

And run the script:

```
$ python ogs_merge/ogs_merge.py -h
usage: ogs_merge.py [-h] [-p PREVIOUS_GFF] [-d DELETED] [-o OUT_PREFIX]
                    genome ogs_name id_regex id_syntax base_gff apollo_gff

positional arguments:
  genome                Genome file (fasta)
  ogs_name              Name of the new OGS
  id_regex              Regex with a capturing group around the incremental
                        part of gene ids, and a second one around the version
                        suffix (e.g.
                        'GSSPF[GPT]([0-9]{8})[0-9]{3}(\.[0-9]+)?')
  id_syntax             String representing a gene id, with {id} where the
                        incremental part of the id should be placed (e.g.
                        'GSSPFG{id}001')
  base_gff              The gff from the base annotation (usually automatic
                        annotation)
  apollo_gff            The gff from the new Apollo valid annotation

optional arguments:
  -h, --help            show this help message and exit
  -p PREVIOUS_GFF, --previous_gff PREVIOUS_GFF
                        The gff from the previous annotation version (if
                        different than <base_gff>)
  -d DELETED, --deleted DELETED
                        File containing a list of mRNAs to remove
  -o OUT_PREFIX, --out_prefix OUT_PREFIX
                        Prefix for output files (default=<ogs_name>_<today's
                        date>)
```

This script generates a GFF file, the transcript/cds/protein fasta files, a file with some statistics, and a tabular file with gene ID correspondance between the automatic annotation and the new one.

## Submitting an OGS to ENA

In gff2embl, there is a script to submit annotations to EBI ENA.

```
$ python gff2embl/gff2embl.py --h
usage: gff2embl.py [-h] -g GENOME -p PROTEINS -s SPECIES -d DESCRIPTION -e
                   EMAIL -j PROJECT [--ref_title REF_TITLE]
                   [--ref_journal REF_JOURNAL] [--ref_authors REF_AUTHORS]
                   [--ref_pubmed_id REF_PUBMED_ID]
                   [--ref_consortium REF_CONSORTIUM]
                   [--no_stop_codon NO_STOP_CODON]
                   [--division {PHG,ENV,FUN,HUM,INV,MAM,VRT,MUS,PLN,PRO,ROD,SYN,TGN,UNC,VRL}]
                   [--out-format {embl-standard,embl-ebi-submit}]
                   [gff] [out]

positional arguments:
  gff                   The gff to read from
  out                   The output embl file, ready for submission to EBI ENA

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        A fasta file containing genome sequence
  -p PROTEINS, --proteins PROTEINS
                        A fasta file containing protein sequences
  -s SPECIES, --species SPECIES
                        The name of the species
  -d DESCRIPTION, --description DESCRIPTION
                        Description of the project
  -e EMAIL, --email EMAIL
                        A valid email address
  -j PROJECT, --project PROJECT
                        A valid EBI study ID (PRJXXXXXXX)
  --ref_title REF_TITLE
                        Title of the reference
  --ref_journal REF_JOURNAL
                        Journal of the reference
  --ref_authors REF_AUTHORS
                        Authors of the reference
  --ref_pubmed_id REF_PUBMED_ID
                        PubMed ID of the reference
  --ref_consortium REF_CONSORTIUM
                        Consortium name of the reference
  --no_stop_codon NO_STOP_CODON
                        Add this option if the protein sequences don't contain
                        trailing stop codons even for complete sequences
  --division {PHG,ENV,FUN,HUM,INV,MAM,VRT,MUS,PLN,PRO,ROD,SYN,TGN,UNC,VRL}
                        The taxonomic division (INV=invertebrate)
  --out-format {embl-standard,embl-ebi-submit}
                        Flavor of EMBL output format: embl-standard=standard
                        EMBL format; embl-ebi-submit=EMBL ready to submit to
                        EBI (some special formating for automatic EBI post-
                        processing)
```

EBI are using a custom embl format with some exotic validation rules.
They provide a tool that performs validation of the produced embl file + automatically reformat some parts of embl files according to their rules.
After running gff2embl, you need to get the latest validator jar from https://mvnrepository.com/artifact/uk.ac.ebi.ena.sequence/embl-api-validator

Then run:

```
java -jar embl-api-validator-1.1.xxx.jar -fix -r ogs.embl
```

This will modify the `ogs.embl` file, be sure to keep a backup.

If you only want to validate the emb file without modifying it, run it wihout the -fix option:

```
java -jar embl-api-validator-1.1.xxx.jar -r ogs.embl
```

## License

Available under the GPLv3 License

Author: Anthony Bretaudeau <anthony.bretaudeau@inra.fr>
