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

These tools can be install with pip:

```bash
pip install gogstools
```

## Merging OGS

This script allows to merge a base annotation with an apollo annotation.
Genes from the base annotation are replaced by overlapping apollo genes, preserving original IDs with a version suffix.

It depends on 4 software that can be installed from conda like this:

```bash
conda create --name ogsmerger bedops==2.4.39 bedtools cufflinks bcbiogff
```

(pinning to bedops 2.4.39 as [version 2.4.40 broke the output of gff2bed](https://github.com/bedops/bedops/issues/244))

Then you can source the conda environment:

```bash
source activate ogsmerger
```

Install gogstools with pip inside the conda env:

```bash
pip install gogstools
```

And run the script:

```bash
usage: ogs_merge [-h] [-p PREVIOUS_GFF] [-d DELETED] [-o OUT_PREFIX] [--exon_parent_regex EXON_PARENT_REGEX] [--exon_parent_to_gene EXON_PARENT_TO_GENE] [--isoform_prefix ISOFORM_PREFIX] [--use_numbers_for_isoform]
                 [--first_isoform_not_numbered] [--no_id_padding] [--id_syntax_rna ID_SYNTAX_RNA]
                 genome ogs_name id_regex id_syntax base_gff apollo_gff

positional arguments:
  genome                Genome file (fasta)
  ogs_name              Name of the new OGS
  id_regex              Regex with a capturing group around the incremental part of gene ids, and a second one around the version suffix (e.g. 'GSSPF[GPT]([0-9]{8})[0-9]{3}(\.[0-9]+)?')
  id_syntax             String representing a gene id, with {id} where the incremental part of the id should be placed (e.g. 'GSSPFG{id}001')
  base_gff              The gff from the base annotation (usually automatic annotation)
  apollo_gff            The gff from the new Apollo valid annotation

optional arguments:
  -h, --help            show this help message and exit
  -p PREVIOUS_GFF, --previous_gff PREVIOUS_GFF
                        The gff from the previous annotation version (if different than <base_gff>)
  -d DELETED, --deleted DELETED
                        File containing a list of mRNAs to remove
  -o OUT_PREFIX, --out_prefix OUT_PREFIX
                        Prefix for output files (default=<ogs_name>_<today's date>)
  --exon_parent_regex EXON_PARENT_REGEX
                        Regex matching exons' Parent ids, with a capturing group around the gene id radical (default='([a-zA-Z0-9]+)([\.0-9]+)?([-_]R[A-Z]+)?(,[a-zA-Z0-9\.\-_]*)?' )
  --exon_parent_to_gene EXON_PARENT_TO_GENE
                        Replacement string to create a gene id from exon_parent_regex first captured group (aka gene id radical), where \1 is the captured group (default='\1' )
  --isoform_prefix ISOFORM_PREFIX
                        Prefix for the isoform part of mRNA ids (default='-R')
  --use_numbers_for_isoform
                        By default, the script will name the isoforms of a gene with letters. If you use this flag, it will be numbers instead.
  --first_isoform_not_numbered
                        If you use this flag, the first isoform will not be numbered (e.g. no -RA).
  --no_id_padding       By default, the numeric part of gene ids will be padded with zeros to a fixed lenght (automatically guessed). If you use this flag, padding is disabled.
  --id_syntax_rna ID_SYNTAX_RNA
                        String representing a gene id, with {id} where the incremental part of the id should be placed (e.g. 'GSSPFR{id}001'). Default is same as id_syntax
```

This script generates a GFF file, the transcript/cds/protein fasta files, a file with some statistics, and a tabular file with gene ID correspondance between the automatic annotation and the new one.

### Usage

This script has evolved to support different kind of gene id formats. As the options are not super easy to understand from scratch, here are 2 examples.

As the creativity in GFF file formatting is a constant source of surprise and discontent, this script may of course fail due to unexpected GFF input file...

#### Default gene id formatting system

By default, this script was initially written to support GFF files with ids looking like this:

- gene features: FOOBAR000123
- mRNA feature, first isoform: FOOBAR000123-RA
- mRNA feature, second isoform: FOOBAR000123-RB
- exon feature: FOOBAR000123-RA-exon1

In this case, the default values should be fine, and you would use regex like this:

- id_regex: `FOOBAR([0-9]{6})(\.[0-9]+)?`
- id_syntax: `FOOBAR{id}`

Genes that get merged with apollo equivalent will have a version suffix appendend to their id, like this: `FOOBAR000123.1`. mRNA, exon, etc features will be derived from this id, including the version suffix.

In case of multiple consecutive merges, this version suffix will be incremented if the script detects any new change from apollo (e.g. `FOOBAR000123.2`)

#### Another gene id formatting system

Another kind of numbering scheme (from EMBL maybe? not sure):

- gene features: gene-FOOBAR123
- mRNA feature, first isoform: rna-FOOBAR123
- mRNA feature, second isoform: rna-FOOBAR123-2
- exon feature: exon-FOOBAR000123-1 or exon-FOOBAR000123-2-1

In this case, the following options should work:

```bash
ogs_merge \
--isoform_prefix "-" \
--use_numbers_for_isoform \
--first_isoform_not_numbered \
--no_id_padding \
--exon_parent_regex '[a-z]+-([a-zA-Z0-9_]+)(-[0-9]+)?([\.0-9]+)?(,[a-zA-Z0-9\.\-_]*)?' \
--exon_parent_to_gene 'gene-\1' \
--id_syntax_rna "rna-FOOBAR{id}" \
genome.fa \
MYOGS2.1 \
"gene-FOOBAR([0-9]+)(\.[0-9]+)?" \
"gene-FOOBAR{id}" \
./base_annotation.gff \
./apollo_export.gff3
```

Genes that get merged with apollo equivalent will have a version suffix appendend to their id, like this: `FOOBAR000123.1`. mRNA, exon, etc features will be derived from this id, including the version suffix.

In case of multiple consecutive merges, this version suffix will be incremented if the script detects any new change from apollo (e.g. `FOOBAR000123.2`)

## Submitting an OGS to ENA

In gff2embl, there is a script to submit annotations to EBI ENA.

```bash
$ gff2embl --h
usage: gff2embl   [-h] -g GENOME -p PROTEINS -s SPECIES -d DESCRIPTION -e
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
After running gff2embl, you need to get the [latest validator jar](https://mvnrepository.com/artifact/uk.ac.ebi.ena.sequence/embl-api-validator)

Then run:

```bash
java -jar embl-api-validator-1.1.xxx.jar -fix -r ogs.embl
```

This will modify the `ogs.embl` file, be sure to keep a backup.

If you only want to validate the emb file without modifying it, run it wihout the -fix option:

```bash
java -jar embl-api-validator-1.1.xxx.jar -r ogs.embl
```

## License

Available under the GPLv3 License

Author: Anthony Bretaudeau <anthony.bretaudeau@inrae.fr>
