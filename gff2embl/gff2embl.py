#!/usr/bin/python

from __future__ import print_function

import argparse
import datetime
import locale
import re
import sys

from BCBio import GFF

from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, Reference, SeqFeature

from writers.BipaaEmblSubmitWriter import BipaaEmblSubmitWriter
from writers.BipaaEmblWriter import BipaaEmblWriter


def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid.  we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(" ", "+").strip()
    search = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
    record = Entrez.read(search)
    return record['IdList'][0]


def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id=taxid, db="taxonomy", retmode="xml")
    return Entrez.read(search)


def get_lineage(species, email):

    Entrez.email = email

    taxid = get_tax_id(species)

    data = get_tax_data(taxid)

    lineage = data[0]['Lineage'].split('; ')

    return lineage


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", help="A fasta file containing genome sequence", type=argparse.FileType('r'), required=True)
parser.add_argument("-p", "--proteins", help="A fasta file containing protein sequences", type=argparse.FileType('r'), required=True)
parser.add_argument("-s", "--species", help="The name of the species", required=True)
parser.add_argument("-d", "--description", help="Description of the project", required=True)
parser.add_argument("-e", "--email", help="A valid email address", required=True)
parser.add_argument("-j", "--project", help="A valid EBI study ID (PRJXXXXXXX) ", required=True)
parser.add_argument("--ref_title", help="Title of the reference")
parser.add_argument("--ref_journal", help="Journal of the reference")
parser.add_argument("--ref_authors", help="Authors of the reference")
parser.add_argument("--ref_pubmed_id", help="PubMed ID of the reference")
parser.add_argument("--ref_consortium", help="Consortium name of the reference")
parser.add_argument("--no_stop_codon", help="Add this option if the protein sequences don't contain trailing stop codons even for complete sequences ")
parser.add_argument("--division", default='INV', choices=['PHG', 'ENV', 'FUN', 'HUM', 'INV', 'MAM', 'VRT', 'MUS', 'PLN', 'PRO', 'ROD', 'SYN', 'TGN', 'UNC', 'VRL'], help="The taxonomic division (INV=invertebrate)")
parser.add_argument("--out-format", choices=['embl-standard', 'embl-ebi-submit'], default='embl-ebi-submit', help="Flavor of EMBL output format: embl-standard=standard EMBL format; embl-ebi-submit=EMBL ready to submit to EBI (some special formating for automatic EBI post-processing)")
parser.add_argument('gff', help="The gff to read from", nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('out', help="The output embl file, ready for submission to EBI ENA", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
args = parser.parse_args()

# First get the lineage (and fail now if not found)
lineage = get_lineage(args.species, args.email)
if not lineage:
    raise RuntimeError("Could not find lineage information on NCBI for species '%s'" % args.species)

print('Found lineage: %s' % lineage)

# Prepare the bibliographic reference
ref = Reference()
if args.ref_pubmed_id:
    ref.pubmed_id = args.ref_pubmed_id
if args.ref_consortium:
    ref.consrtm = args.ref_consortium
if args.ref_authors:
    ref.authors = args.ref_authors
if args.ref_title:
    ref.title = args.ref_title
if args.ref_journal:
    ref.journal = args.ref_journal
else:
    now = datetime.datetime.now()
    ref_date = now.strftime("%m-%b-%Y").upper()
    # Temp switch to C to get english month abbr
    saved = locale.setlocale(locale.LC_TIME)
    try:
        locale.setlocale(locale.LC_TIME, "C")
        ref_date = now.strftime("%m-%b-%Y").upper()
    finally:
        locale.setlocale(locale.LC_TIME, saved)
    ref.journal = "Submitted (" + ref_date + ") to the INSDC."

print('Loading input GFF and fasta files...')

seq_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta", alphabet=generic_dna))

prot_seq_dict = SeqIO.to_dict(SeqIO.parse(args.proteins, "fasta"))

gff_iter = GFF.parse(args.gff, base_dict=seq_dict)

# A custom writer as the one bundled in biopython has some limitations
if args.out_format == 'embl-standard':
    SeqIO._FormatToWriter['embl'] = BipaaEmblWriter
else:
    SeqIO._FormatToWriter['embl'] = BipaaEmblSubmitWriter

print('Parsing GFF...')

# To ease debugging
convert_only = None  # A list of gene ids to convert, set to None to convert everything

for rec in gff_iter:

    # Add a source feature corresponding to current scaffold
    q = {}
    q['mol_type'] = 'genomic DNA'
    q['organism'] = args.species
    q['note'] = rec.name
    source_f = SeqFeature(FeatureLocation(0, len(rec)), type="source", qualifiers=q)
    new_feats = [source_f]

    keep_rec = False

    for f in rec.features:  # gene

        gene_quals = {}
        locus_tag = f.qualifiers['ID'][0]

        # Debugging code
        if convert_only:
            if locus_tag in convert_only:
                keep_rec = True
            else:
                keep_rec = keep_rec and False
                continue
        else:
            keep_rec = True

        locus_tag = re.sub(r"^([a-zA-Z]+)([0-9]+)$", r"\1_\2", locus_tag)  # EBI asks locus_tag to be of the form: XXXX_00000
        gene_quals['locus_tag'] = locus_tag
        gene_quals['gene'] = locus_tag
        f.qualifiers = gene_quals

        # See if there is a human readable name or symbol for the product
        product_name = locus_tag
        if 'Name' in f.qualifiers and f.qualifiers['Name'][0] and not (len(f.qualifiers['Name'][0]) == 32 and re.match("^[A-F0-9]+$", f.qualifiers['Name'][0])):
            product_name = f.qualifiers['Name'][0].strip()
        elif 'symbol' in f.qualifiers and f.qualifiers['symbol'][0]:
            product_name = f.qualifiers['symbol'][0].strip()

        all_dbxref = []
        if 'Dbxref' in f.qualifiers and f.qualifiers['Dbxref'][0]:
            for dbxref in f.qualifiers['Dbxref']:
                splitted_dbxref = dbxref.split(":")
                db = splitted_dbxref[0].strip()
                for x in splitted_dbxref[1:]:
                    all_dbxref.append(db + ':' + x.strip())

        new_feats.append(f)
        seen_cds_locs = []
        for sf in f.sub_features:
            if sf.type == "mRNA":

                cds_locs = []
                utr_feats = []
                for ssf in sf.sub_features:
                    if ssf.type in ['CDS']:
                        cds_locs.append(ssf.location)
                    elif ssf.type in ['five_prime_UTR']:
                        utr_quals = {}
                        utr_quals['gene'] = locus_tag
                        utr_quals['note'] = ("utr_id=%s" % ssf.qualifiers['ID'][0])
                        new_utr = SeqFeature(ssf.location, type="5'UTR", qualifiers=utr_quals)
                        utr_feats.append(new_utr)
                    elif ssf.type in ['three_prime_UTR']:
                        utr_quals = {}
                        utr_quals['gene'] = locus_tag
                        utr_quals['note'] = ("utr_id=%s" % ssf.qualifiers['ID'][0])
                        new_utr = SeqFeature(ssf.location, type="3'UTR", qualifiers=utr_quals)
                        utr_feats.append(new_utr)
                    elif ssf.type in ['UTR']:
                        utr_side = "5'UTR"
                        if len(cds_locs) > 0:
                            utr_side = "3'UTR"
                        utr_quals = {}
                        utr_quals['gene'] = locus_tag
                        utr_quals['note'] = ("utr_id=%s" % ssf.qualifiers['ID'][0])
                        new_utr = SeqFeature(ssf.location, type=utr_side, qualifiers=utr_quals)
                        utr_feats.append(new_utr)

                # cds_locs should be sorted in forward direction (even if on reverse strand)
                cds_locs = sorted(cds_locs, key=lambda loc: loc.start)

                mrna_id = sf.qualifiers['ID'][0]
                pep_id = re.sub(r'-R([A-Z]+)', r'-P\1', mrna_id)

                mrna_quals = {}
                mrna_quals['locus_tag'] = locus_tag
                mrna_quals['gene'] = locus_tag
                mrna_quals['note'] = ("transcript_id=%s" % mrna_id)

                cds_quals = {}
                cds_quals['transl_table'] = 1
                cds_quals['locus_tag'] = locus_tag
                cds_quals['gene'] = locus_tag
                cds_quals['product'] = product_name
                cds_quals['note'] = ("transcript_id=%s" % mrna_id)
                cds_quals['note'] = ("protein_id=%s" % pep_id)

                for x in all_dbxref:
                    mrna_quals['db_xref'] = x

                # Find and fix peptide sequence
                if pep_id not in prot_seq_dict:
                    raise RuntimeError("Could not find protein sequence for id '%s'" % pep_id)

                pep_seq = str(prot_seq_dict[pep_id].seq)

                # If the protein doesn't start with methionine, it is probably a partial one
                fuzzy_start = False
                fuzzy_end = False
                if pep_seq[0] != 'M':
                    fuzzy_start = True
                    if sf.strand == 1:
                        cds_locs[0] = FeatureLocation(BeforePosition(cds_locs[0].start), cds_locs[0].end, sf.strand)
                        cds_quals['codon_start'] = 1
                    else:
                        cds_locs[-1] = FeatureLocation(cds_locs[-1].start, AfterPosition(cds_locs[-1].end), sf.strand)
                        cds_quals['codon_start'] = 1
                # If the protein doesn't end with stop codon, it is probably a partial one
                if not args.no_stop_codon and pep_seq[-1] not in ('.', '*'):
                    fuzzy_end = True
                    if sf.strand == 1:
                        cds_locs[-1] = FeatureLocation(cds_locs[-1].start, AfterPosition(cds_locs[-1].end), sf.strand)
                        cds_quals['codon_start'] = 1
                    else:
                        cds_locs[0] = FeatureLocation(BeforePosition(cds_locs[0].start), cds_locs[0].end, sf.strand)
                        cds_quals['codon_start'] = 1

                # Check that there is not already an identical mRNA
                if str(cds_locs) in str(seen_cds_locs):
                    print("mRNA %s on %s is identical to a previous one. Skipping." % (mrna_id, rec.name), file=sys.stderr)
                    continue

                seen_cds_locs.append(cds_locs)

                last_char = str(pep_seq[-1:])
                if last_char == '.' or last_char == '*':
                    pep_seq = pep_seq[:-1]

                # Replace stop codons in the middle of sequences too
                if '.' in pep_seq or '*' in pep_seq:
                    # As it occurs for manually edited genes, mark this as an exception (to avoid validation errors)
                    cds_quals['exception'] = "reasons given in citation"
                    pep_seq = pep_seq.replace('.', 'X')
                    pep_seq = pep_seq.replace('*', 'X')
                    print("Gene %s on %s contains stop codon. Marking as exception." % (locus_tag, rec.name), file=sys.stderr)

                # Check that sequence length is a multiple of 3
                cds_length = 0
                for loc in cds_locs:
                    cds_length += len(loc)

                if cds_length % 3 != 0 and not fuzzy_start:
                    # As it occurs for manually edited genes, mark this as an exception (to avoid validation errors)
                    cds_quals['exception'] = "reasons given in citation"
                    print("Gene %s on %s has a length which is not a multiple of 3. Marking as exception." % (locus_tag, rec.name), file=sys.stderr)

                cds_quals['translation'] = pep_seq

                if len(cds_locs) > 1:
                    # join cds locations
                    joined_loc = CompoundLocation(cds_locs)
                else:
                    joined_loc = cds_locs[0]

                mrna_joined_feature = SeqFeature(joined_loc, type="mRNA", qualifiers=mrna_quals)
                cds_joined_feature = SeqFeature(joined_loc, type="CDS", qualifiers=cds_quals)

                new_feats.append(mrna_joined_feature)
                new_feats.append(cds_joined_feature)

                new_feats += utr_feats

    if keep_rec:
        rec.features = new_feats

        rec.description = args.description

        rec.annotations['organism'] = args.species
        rec.annotations['taxonomy'] = lineage
        rec.annotations['data_file_division'] = args.division

        ref.location = [FeatureLocation(0, len(rec))]
        rec.annotations['references'] = [ref]

        rec.dbxrefs = [('Project:%s' % args.project)]

        rec.annotations['keywords'] = ['CON']  # CON is appropriate for scaffolds: https://www.ebi.ac.uk/training/online/course/nucleotide-sequence-data-resources-ebi/what-ena/how-sequence-assembled

        SeqIO.write(rec, args.out, "embl")
