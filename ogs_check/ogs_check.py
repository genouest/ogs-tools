#!/usr/bin/env python

# Check that an OGS GFF file is ready for release

import argparse
import logging
import sys

from BCBio import GFF


logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


class OgsCheck():

    def __init__(self):
        self.gene_ids = []
        self.mRNA_ids = []
        self.exon_ids = []
        self.skipped_types = set()

    def check_valid_mrna(self, mrna):

        if mrna.type == 'transcript':
            mrna.type = "mRNA"

        if mrna.type != 'mRNA':
            self.skipped_types.add(mrna.type)
            return None

        if 'ID' not in mrna.qualifiers or len(mrna.qualifiers['ID']) == 0:
            log.error("Found an mRNA without an ID attribute")
            return None

        if len(mrna.qualifiers['ID']) != 1:
            log.error("Found an mRNA with too many ID attributes")
            return None

        if mrna.qualifiers['ID'][0] in self.mRNA_ids:
            log.error("Duplicate mRNA id: %s" % mrna.qualifiers['ID'][0])
            return None

        self.mRNA_ids.append(mrna.qualifiers['ID'][0])

        exon_coords = {}
        cds_cumul = 0

        # Find positions
        kept_gchild = []
        self.exon_ids = []
        for gchild in mrna.sub_features:  # exons, cds, utr

            if gchild.type == "exon":
                exon_coords[gchild.location.start] = gchild.location.end
            elif gchild.type == "CDS":
                cds_cumul += gchild.location.end - gchild.location.start - 1

            if gchild.type in ['five_prime_utr', "5'UTR"]:
                gchild.type = 'five_prime_UTR'

            if gchild.type in ['three_prime_utr', "3'UTR"]:
                gchild.type = 'three_prime_UTR'

            if gchild.type in ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']:
                kept_gchild.append(gchild)
            else:
                self.skipped_types.add(gchild.type)

            if gchild.type == "exon":
                if 'ID' not in gchild.qualifiers or len(gchild.qualifiers['ID']) == 0:
                    log.error("Found an exon without an ID attribute")
                    return None

                if len(gchild.qualifiers['ID']) != 1:
                    log.error("Found an exon with too many ID attributes")
                    return None

                if gchild.qualifiers['ID'][0] in self.exon_ids:
                    log.error("Duplicate exon id: %s" % mrna.qualifiers['ID'][0])
                    return None

                self.exon_ids.append(gchild.qualifiers['ID'][0])

        mrna.sub_features = kept_gchild

        # Check minimum intron size
        start_sorted = sorted(exon_coords)
        previous_end = None
        for exon_start in start_sorted:
            if previous_end is not None:
                intron_size = exon_start - previous_end
                if intron_size < 9:
                    log.warning("Discarding '%s' because intron size %s < 9" % (mrna.qualifiers['ID'][0], intron_size))
                    return None

            previous_end = exon_coords[exon_start]

        # Check minimum cds size
        if cds_cumul < 15:
            log.warning("Discarding '%s' because CDS size < 15" % mrna.qualifiers['ID'][0])
            return None

        return mrna

    def check(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
        parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        args = parser.parse_args()

        scaffs = []
        for scaff in GFF.parse(args.infile):
            scaff.annotations = {}
            scaff.seq = ""
            new_genes = []

            for gene in scaff.features:

                if gene.type != 'gene':
                    self.skipped_types.add(gene.type)
                    continue

                if 'ID' not in gene.qualifiers or len(gene.qualifiers['ID']) == 0:
                    log.error("Found a gene without an ID attribute")
                    continue

                if len(gene.qualifiers['ID']) != 1:
                    log.error("Found a gene with too many ID attributes")
                    continue

                if gene.qualifiers['ID'][0] in self.gene_ids:
                    log.error("Duplicate gene id: %s" % gene.qualifiers['ID'][0])
                    continue

                self.gene_ids.append(gene.qualifiers['ID'][0])

                new_mrnas = []

                for mrna in gene.sub_features:

                    mrna = self.check_valid_mrna(mrna)

                    if mrna is not None:
                        new_mrnas.append(mrna)

                gene.sub_features = new_mrnas
                new_genes.append(gene)

            scaff.features = new_genes

            if len(new_genes):
                scaffs.append(scaff)

        GFF.write(scaffs, args.outfile)

        if self.skipped_types:
            log.warning("Skipped unknown/misplaced feature types: %s" % (self.skipped_types))


if __name__ == '__main__':
    ogsc = OgsCheck()

    ogsc.check()
