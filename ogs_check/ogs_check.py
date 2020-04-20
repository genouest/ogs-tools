#!/usr/bin/env python

# Check that an OGS GFF file is ready for release

import argparse
import logging
import sys

from BCBio import GFF

from Bio.SeqFeature import FeatureLocation, SeqFeature


logging.basicConfig(level=logging.INFO)
log = logging.getLogger()


def change_parentname(feature, parentKeyName, parentName):

    for child in feature.sub_features:
        child.qualifiers[parentKeyName][0] = parentName

    return


class OgsCheck():

    def __init__(self):
        self.mRNA_ids = []
        self.exon_ids = []
        self.skipped_types = set()
        self.qlistName = ['Name', 'ID']

    def check_valid_mrna(self, mrna, is_complete=True):

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

        if 'Name' not in mrna.qualifiers:
            mrna.qualifiers['Name'] = mrna.qualifiers['ID']

        if is_complete:
            self.mRNA_ids.append(mrna.qualifiers['ID'][0])

        exon_coords = {}
        cds_cumul = 0

        # Find positions
        kept_gchild = []
        self.exon_ids = []

        for gchild in mrna.sub_features:  # exons, cds, utr

            if self.args.source:
                gchild.qualifiers['source'][0] = self.args.source

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

        # Only check CDS/intron sizes when we're sure the mrna is complete
        if is_complete:
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
            if cds_cumul < 18:
                log.warning("Discarding '%s' because CDS size < 18 (%s)" % (mrna.qualifiers['ID'][0], cds_cumul))
                return None

        if self.args.source:
            mrna.qualifiers['source'][0] = self.args.source

        return mrna

    def find_inferred_parents(self, features):
        inferred = {}
        for topfeat in features:
            if topfeat.type == 'inferred_parent':
                inferred[topfeat.qualifiers['ID'][0]] = topfeat

        return inferred

    def create_parent(self, orphan, parent_id, orphan_id, parent_type):
        q = {}
        for key in orphan.qualifiers:
            q[key] = list(orphan.qualifiers[key])
        new_parent = SeqFeature(FeatureLocation(orphan.location.start, orphan.location.end), type=parent_type, strand=orphan.location.strand, qualifiers=q)
        for qn in self.qlistName:
            if qn in new_parent.qualifiers:
                new_parent.qualifiers[qn][0] = parent_id
        for qn in self.qlistName:
            if qn in orphan.qualifiers:
                # The new panret is assigned the id from the orphan, and the orphan id might be modified
                orphan.qualifiers[qn][0] = orphan_id
        new_parent.sub_features = []
        new_parent.sub_features.append(orphan)
        orphan.qualifiers['Parent'] = new_parent.qualifiers['ID']
        if 'Parent' in new_parent.qualifiers:
            del new_parent.qualifiers['Parent']
        change_parentname(orphan, 'Parent', orphan.qualifiers['ID'][0])

        if self.args.source:
            new_parent.qualifiers['source'][0] = self.args.source

        return new_parent

    def adopt_orphan_mrna(self, orphan, is_complete=True):
        # Validate it, create a gene parent, and look if we have a corresponding inferred_parent containing children from this mRNA
        if 'Parent' in orphan.qualifiers and len(orphan.qualifiers['Parent']) == 1:
            parent_id = orphan.qualifiers['Parent'][0]
            orphan_id = orphan.qualifiers['ID'][0]
        else:
            parent_id = orphan.qualifiers['ID'][0]
            orphan_id = parent_id + '-R'

        if 'ID' in orphan.qualifiers and len(orphan.qualifiers['ID']) == 1:
            if len(orphan.sub_features) == 0 and orphan.qualifiers['ID'][0] in self.inferred_parents:
                orphan.sub_features = self.inferred_parents[orphan.qualifiers['ID'][0]].sub_features
                del self.inferred_parents[orphan.qualifiers['ID'][0]]

        orphan = self.check_valid_mrna(orphan, is_complete)

        if orphan is not None:

            if parent_id in self.new_genes:
                potential_parent = self.new_genes[parent_id]

                if potential_parent.location.strand != orphan.location.strand:
                    log.error("Conflict between an orphan %s and its potential parent %s strand: %s != %s" % (orphan.type, parent_id, orphan.location.strand, potential_parent.location.strand))
                    return None

                potential_parent.sub_features.append(orphan)

                if potential_parent.location.start > orphan.location.start:
                    potential_parent.location = FeatureLocation(orphan.location.start, potential_parent.location.end, strand=potential_parent.location.strand)

                if potential_parent.location.end < orphan.location.end:
                    potential_parent.location = FeatureLocation(potential_parent.location.start, orphan.location.end, strand=potential_parent.location.strand)

                self.new_genes[parent_id] = potential_parent

            else:
                new_g = self.create_parent(orphan, parent_id, orphan_id, 'gene')
                self.new_genes[parent_id] = new_g

            self.all_mrnas[orphan.qualifiers['ID'][0]] = orphan

        return orphan

    def adopt_orphan_exoncds(self, orphan, last_one=True):
        # Validate it, create a gene parent, and look if we have a corresponding inferred_parent containing children from this mRNA
        if 'Parent' in orphan.qualifiers and len(orphan.qualifiers['Parent']) == 1:
            parent_id = orphan.qualifiers['Parent'][0]
            orphan_id = orphan.qualifiers['ID'][0]
        else:
            parent_id = orphan.qualifiers['ID'][0]
            orphan_id = '%s-%s' % (parent_id, orphan.type)

        if 'ID' in orphan.qualifiers and len(orphan.qualifiers['ID']) == 1:
            if len(orphan.sub_features) == 0 and orphan.qualifiers['ID'][0] in self.inferred_parents:
                orphan.sub_features = self.inferred_parents[orphan.qualifiers['ID'][0]].sub_features
                del self.inferred_parents[orphan.qualifiers['ID'][0]]

        if parent_id in self.all_mrnas:
            potential_parent = self.all_mrnas[parent_id]

            if potential_parent.location.strand != orphan.location.strand:
                log.error("Conflict between an orphan %s and its potential parent %s strand: %s != %s" % (orphan.type, parent_id, orphan.location.strand, potential_parent.location.strand))
                return None

            del orphan.qualifiers['Parent']  # previous parent is no longer parent
            potential_parent.sub_features.append(orphan)

            if potential_parent.location.start > orphan.location.start:
                potential_parent.location = FeatureLocation(orphan.location.start, potential_parent.location.end, strand=potential_parent.location.strand)

            if potential_parent.location.end < orphan.location.end:
                potential_parent.location = FeatureLocation(potential_parent.location.start, orphan.location.end, strand=potential_parent.location.strand)

            potential_parent = self.check_valid_mrna(potential_parent, last_one)

            if potential_parent is None:
                return None

            self.all_mrnas[parent_id] = potential_parent

            # update its gene parent
            gene_children = []
            for mrna in self.new_genes[potential_parent.qualifiers['Parent'][0]].sub_features:
                if mrna.qualifiers['ID'][0] == parent_id:
                    gene_children.append(potential_parent)
                else:
                    gene_children.append(mrna)
            self.new_genes[potential_parent.qualifiers['Parent'][0]].sub_features = gene_children

            # Fix gene location
            if self.new_genes[potential_parent.qualifiers['Parent'][0]].location.start > potential_parent.location.start:
                self.new_genes[potential_parent.qualifiers['Parent'][0]].location = FeatureLocation(potential_parent.location.start, self.new_genes[potential_parent.qualifiers['Parent'][0]].location.end, strand=potential_parent.location.strand)
            if self.new_genes[potential_parent.qualifiers['Parent'][0]].location.end < potential_parent.location.end:
                self.new_genes[potential_parent.qualifiers['Parent'][0]].location = FeatureLocation(self.new_genes[potential_parent.qualifiers['Parent'][0]].location.start, potential_parent.location.end, strand=potential_parent.location.strand)

        else:
            new_mRNA = self.create_parent(orphan, parent_id, orphan_id, "mRNA")
            self.all_mrnas[parent_id] = new_mRNA

            self.adopt_orphan_mrna(new_mRNA, is_complete=last_one)

        return orphan

    def check(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
        parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
        parser.add_argument('--source', help="Change the source to given value for all features")
        self.args = parser.parse_args()

        scaffs = []
        for scaff in GFF.parse(self.args.infile):
            scaff.annotations = {}
            scaff.seq = ""

            # Genes and mRNA list, reset on each new scaff
            self.new_genes = {}
            self.all_mrnas = {}

            # First check if we have inferred_parent (generated by bcbio-gff)
            self.inferred_parents = self.find_inferred_parents(scaff.features)

            for topfeat in scaff.features:

                if topfeat.type not in ['gene', 'mRNA', 'CDS', 'exon']:
                    if topfeat.type != 'inferred_parent':
                        self.skipped_types.add(topfeat.type)
                    continue

                if 'ID' not in topfeat.qualifiers or len(topfeat.qualifiers['ID']) == 0:
                    log.error("Found a top level %s feature without an ID attribute" % topfeat.type)
                    continue

                if len(topfeat.qualifiers['ID']) != 1:
                    log.error("Found a top level %s feature with too many ID attributes" % topfeat.type)
                    continue

                if topfeat.qualifiers['ID'][0] in self.new_genes.keys():
                    log.error("Duplicate top level %s feature id: %s" % (topfeat.qualifiers['ID'][0], topfeat.type))
                    continue

                if topfeat.type == 'gene':
                    # Simple case: a gene with sub features
                    new_mrnas = []

                    for mrna in topfeat.sub_features:

                        mrna = self.check_valid_mrna(mrna)

                        if mrna is not None:
                            new_mrnas.append(mrna)
                            self.all_mrnas[mrna.qualifiers['ID'][0]] = mrna

                    topfeat.sub_features = new_mrnas

                    if self.args.source:
                        topfeat.qualifiers['source'][0] = self.args.source

                    if 'Name' not in topfeat.qualifiers:
                        topfeat.qualifiers['Name'] = topfeat.qualifiers['ID']

                    self.new_genes[topfeat.qualifiers['ID'][0]] = topfeat

                elif topfeat.type == 'mRNA':
                    # Found an mRNA without gene parent
                    self.adopt_orphan_mrna(topfeat)

                elif topfeat.type in ['exon', 'CDS']:
                    # Found an exon/cds without gene parent
                    self.adopt_orphan_exoncds(topfeat)
                else:
                    log.error('Unexpected feature type %s. There is bug.' % topfeat.type)

            # Now handle the remaining inferred_parents
            for topfeat_name in self.inferred_parents:
                topfeat = self.inferred_parents[topfeat_name]

                if len(topfeat.sub_features) < 1:
                    log.error("Skipping an inferred_parent without children %s" % topfeat)
                    continue

                guessed_type = None
                if topfeat.sub_features[0].type in ['exon', 'CDS', 'start_codon', 'stop_codon', "5'UTR", 'five_prime_UTR', 'five_prime_utr', "3'UTR", 'three_prime_UTR', 'three_prime_utr']:
                    guessed_type = 'mRNA'
                elif topfeat.sub_features[0].type == 'mRNA':
                    guessed_type = 'gene'
                else:
                    log.error("Skipping an inferred_parent: failed to guess type %s" % topfeat)
                    continue

                if guessed_type == 'mRNA':
                    num_seen = 0
                    for sub in topfeat.sub_features:
                        num_seen += 1
                        last_one = num_seen == len(topfeat.sub_features)
                        self.adopt_orphan_exoncds(sub, last_one=last_one)

                elif guessed_type == 'gene':
                    for sub in topfeat.sub_features:
                        self.adopt_orphan_mrna(sub)
                else:
                    log.error('Unexpected feature type %s. There is bug.' % topfeat.type)

            scaff.features = self.new_genes.values()

            if len(self.new_genes):
                scaffs.append(scaff)

        GFF.write(scaffs, self.args.outfile)

        if self.skipped_types:
            log.warning("Skipped unknown/misplaced feature types: %s" % (self.skipped_types))


if __name__ == '__main__':
    ogsc = OgsCheck()

    ogsc.check()
