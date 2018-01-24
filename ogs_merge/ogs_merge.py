#!/usr/bin/python

import argparse
import datetime
import os
import re
import shutil
import string
import sys
import tempfile

from subprocess import call

from BCBio import GFF

from Bio.SeqFeature import FeatureLocation, SeqFeature

# TODO handle (rare) cases where the id prefix vary depending on the feature type (e.g.: GSSPFG for genes, GSSPFT for transcripts, GSSPFP for proteins)


class OgsMerger():

    # Parse a gene feature and create a hash representing its content (structure and attributes)
    # we can not always trust date_last_modified as some fix scripts don't change its value
    # instead, build a hash representing the gene content (exons, attributes, ...)
    def make_content_hash(self, f, check_name=False):
        blacklist_attr = ['ID', 'source', 'Parent', 'owner', 'Alias', 'date_creation', 'date_last_modified']
        if not check_name:
            blacklist_attr += ['Name']

        fhash = str(f.location) + "__"
        for qk in sorted(f.qualifiers):
            if qk not in blacklist_attr:
                fhash += qk + "->" + str(sorted(f.qualifiers[qk], key=str.lower)) + "__"

        # we need to sort subfeatures to have comparable hashes
        subs = sorted(f.sub_features, key=lambda e: (e.type, e.location.start, e.location.end))
        for sub in subs:
            # We don't pass the check_name as we don't want to check names on sub features (they are most of the time automatic)
            fhash += "__SUB[" + self.make_content_hash(sub) + "]"

        return fhash

    # Parse a gene feature and create a hash representing its content (structure only)
    def make_structure_hash(self, f):
        fhash = str(f.location) + "__"

        # we need to sort subfeatures to have comparable hashes
        subs = sorted(f.sub_features, key=lambda e: (e.type, e.location.start, e.location.end))
        for sub in subs:
            if sub.type in ['mRNA', 'CDS']:
                fhash += "__SUB[" + self.make_structure_hash(sub) + "]"

        return fhash

    # Clean a gene (and its children) attributes to output in the final GFF file
    def clean_feature(self, f):
        gene_id = f.qualifiers['ID'][0]
        f.qualifiers['ID'][0] = gene_id

        f.qualifiers['source'][0] = self.source
        if 'filtertag' in f.qualifiers:
            del f.qualifiers['filtertag']
        if 'owner' in f.qualifiers:
            del f.qualifiers['owner']

        # First letter in capital is not valid for GFF3 (no longer a problem in apollo 2.0.7, kept for compatibility with older versions)
        if 'Allele' in f.qualifiers:
            f.qualifiers['allele'] = f.qualifiers['Allele']
            del f.qualifiers['Allele']
        if 'Part' in f.qualifiers:
            f.qualifiers['part'] = f.qualifiers['Part']
            del f.qualifiers['Part']
        if 'Synonym' in f.qualifiers:
            f.qualifiers['synonym'] = f.qualifiers['Synonym']
            del f.qualifiers['Synonym']

        # Some qualifiers are not needed outside apollo
        if 'status' in f.qualifiers:
            del f.qualifiers['status']
        if 'annotGroup' in f.qualifiers:
            del f.qualifiers['annotGroup']

        # Cleanup some attributes
        if 'symbol' in f.qualifiers:
            f.qualifiers['symbol'] = [x.strip() for x in f.qualifiers['symbol']]
        if 'description' in f.qualifiers:
            f.qualifiers['description'] = [x.strip() for x in f.qualifiers['description']]
        if 'Note' in f.qualifiers:
            f.qualifiers['Note'] = [x.strip() for x in f.qualifiers['Note']]
        if 'Dbxref' in f.qualifiers:
            f.qualifiers['Dbxref'] = [x.strip() for x in f.qualifiers['Dbxref']]
        if 'allele' in f.qualifiers:
            f.qualifiers['allele'] = [x.strip() for x in f.qualifiers['allele']]
        if 'part' in f.qualifiers:
            f.qualifiers['part'] = [x.strip() for x in f.qualifiers['part']]
        if 'Name' in f.qualifiers:
            f.qualifiers['full_name'] = [x.strip() for x in f.qualifiers['Name']]
        if 'synonym' in f.qualifiers:
            f.qualifiers['synonym'] = [x.strip() for x in f.qualifiers['synonym']]

        mrna_count = 0
        mrna_count_cycle = 1

        # Add previous versions as aliases
        gene_id_splitted = gene_id.split(".")
        if len(gene_id_splitted) > 1:
            gene_id_no_version = gene_id_splitted[0]
            gene_version = int(gene_id_splitted[1])
            if 'Alias' not in f.qualifiers:
                f.qualifiers['Alias'] = []
            f.qualifiers['Alias'].append(gene_id_no_version)
            for old_version in range(1, gene_version):
                f.qualifiers['Alias'].append(gene_id_no_version + '.' + str(old_version))

        has_multiple_isoforms = len(f.sub_features) > 1

        for child in f.sub_features:  # mRNA
            child.qualifiers['source'][0] = self.source
            if 'filtertag' in child.qualifiers:
                del child.qualifiers['filtertag']
            if 'owner' in child.qualifiers:
                del child.qualifiers['owner']
            mrna_id = gene_id + "-" + "R" + (string.ascii_uppercase[mrna_count] * mrna_count_cycle)
            if 'Name' not in child.qualifiers:
                child.qualifiers['Name'] = [mrna_id]
            else:
                child.qualifiers['Name'][0] = mrna_id

            # put apollo id to aliases and write a good ID
            current_id = child.qualifiers['ID'][0]
            if re.match("^[A-F0-9]{32}$", current_id) or re.match("^[a-f0-9-]{36}$", current_id):
                if 'Alias' not in child.qualifiers:
                    child.qualifiers['Alias'] = []
                child.qualifiers['Alias'].append(current_id)

            if 'ID' not in child.qualifiers:
                child.qualifiers['ID'] = [mrna_id]
            else:
                child.qualifiers['ID'][0] = mrna_id

            child.qualifiers['Parent'][0] = gene_id

            # Some qualifiers are not needed outside apollo
            if 'status' in child.qualifiers:
                del child.qualifiers['status']
            if 'annotGroup' in child.qualifiers:
                del child.qualifiers['annotGroup']

            # Transfer attributes to mRNA level
            if 'symbol' in f.qualifiers:
                child.qualifiers['symbol'] = f.qualifiers['symbol']
            if 'allele' in f.qualifiers:
                child.qualifiers['allele'] = f.qualifiers['allele']
            if 'part' in f.qualifiers:
                child.qualifiers['part'] = f.qualifiers['part']
            if 'synonym' in f.qualifiers:
                if 'synonym' not in child.qualifiers:
                    child.qualifiers['synonym'] = []
                child.qualifiers['synonym'] += f.qualifiers['synonym']
            if 'Note' in f.qualifiers:
                if 'Note' not in child.qualifiers:
                    child.qualifiers['Note'] = []
                child.qualifiers['Note'] += f.qualifiers['Note']
            if 'Dbxref' in f.qualifiers:
                if 'Dbxref' not in child.qualifiers:
                    child.qualifiers['Dbxref'] = []
                child.qualifiers['Dbxref'] += f.qualifiers['Dbxref']

            # Transfer only if not multiple isoforms
            if not has_multiple_isoforms:
                if 'Name' in f.qualifiers:
                    child.qualifiers['full_name'] = f.qualifiers['Name']

            # Transfer without replacing isoform information
            if not has_multiple_isoforms or 'description' not in child.qualifiers or not child.qualifiers['description']:
                if 'description' in f.qualifiers:
                    child.qualifiers['description'] = f.qualifiers['description']

            # Add previous versions as aliases
            if len(gene_id_splitted) > 1:
                if 'Alias' not in child.qualifiers:
                    child.qualifiers['Alias'] = []
                child.qualifiers['Alias'].append(gene_id_no_version + "-" + "R" + (string.ascii_uppercase[mrna_count] * mrna_count_cycle))
                for old_version in range(1, gene_version):
                    child.qualifiers['Alias'].append(gene_id_no_version + '.' + str(old_version) + "-" + "R" + (string.ascii_uppercase[mrna_count] * mrna_count_cycle))

            # Remove uppercase variants if any
            if 'Allele' in child.qualifiers:
                del child.qualifiers['Allele']
            if 'Part' in child.qualifiers:
                del child.qualifiers['Part']
            if 'Synonym' in child.qualifiers:
                del child.qualifiers['Synonym']

            id_count = 1
            for gchild in child.sub_features:  # exons, cds, ...
                gchild.qualifiers['source'][0] = self.source
                if 'filtertag' in gchild.qualifiers:
                    del gchild.qualifiers['filtertag']
                if 'ID' in gchild.qualifiers:
                    gchild.qualifiers['ID'][0] = mrna_id + "-" + gchild.type + "-" + str(id_count)
                if 'owner' in gchild.qualifiers:
                    del gchild.qualifiers['owner']
                gchild.qualifiers['Name'] = [mrna_id + "-" + gchild.type]  # Will create the Name array if not yet present
                gchild.qualifiers['Parent'][0] = mrna_id
                id_count_gg = 1
                for ggchild in gchild.sub_features:  # exotic stuff (non_canonical_five_prime_splice_site non_canonical_three_prime_splice_site stop_codon_read_through)
                    ggchild.qualifiers['source'][0] = self.source
                    if 'filtertag' in ggchild.qualifiers:
                        del ggchild.qualifiers['filtertag']
                    if 'owner' in ggchild.qualifiers:
                        del ggchild.qualifiers['owner']
                    ggchild.qualifiers['ID'][0] = mrna_id + "-" + gchild.type + "-" + str(id_count) + "-" + ggchild.type + "-" + str(id_count_gg)
                    ggchild.qualifiers['Name'][0] = mrna_id + "-" + gchild.type + "-" + ggchild.type
                    ggchild.qualifiers['Parent'][0] = mrna_id + "-" + gchild.type + "-" + str(id_count)
                    id_count_gg += 1
                id_count += 1

            mrna_count += 1
            if mrna_count >= 25:
                mrna_count = 0
                mrna_count_cycle += 1

        return f

    # Add exons subfeatures to a feature, guessing from UTR and CDS subfeatures
    def guess_exons(self, cleaned_f):

        for child in cleaned_f.sub_features:  # mRNA
            cds_coords = {}
            utr_coords = {}
            exon_coords = {}
            seen_exons = False
            for gchild in child.sub_features:  # exons, cds, ...
                if gchild.type == "CDS":
                    cds_coords[gchild.location.start] = gchild.location.end
                elif gchild.type in ["UTR", "five_prime_UTR", "three_prime_UTR"]:
                    utr_coords[gchild.location.start] = gchild.location.end
                elif gchild.type == "exon":
                    # We found some exons, it means we don't need to guess them
                    seen_exons = True

            if seen_exons:
                continue

            cds_sorted = sorted(cds_coords)
            for cds_start in cds_sorted:
                cds_end = cds_coords[cds_start]
                utr_len_init = len(utr_coords)

                if (cds_end) in utr_coords:
                    # 3' UTR
                    exon_coords[cds_start] = utr_coords[cds_end]
                    del utr_coords[cds_end]

                if cds_start in utr_coords.values():
                    # 5' UTR
                    found_start = None
                    for start, end in utr_coords.items():
                        if found_start is None and end == cds_start:
                            found_start = start
                    if cds_start in exon_coords:
                        # we found utr on both sides
                        cds_end = exon_coords[cds_start]
                        del exon_coords[cds_start]

                    exon_coords[found_start] = cds_end

                    del utr_coords[found_start]

                if len(utr_coords) == utr_len_init:
                    exon_coords[cds_start] = cds_end

            # Add all remaining utr = utr that are not adjacent to cds
            for utr_start in utr_coords:
                exon_coords[utr_start] = utr_coords[utr_start]

            count_id = 1
            for exon_start in exon_coords:
                new_subf = SeqFeature(FeatureLocation(exon_start, exon_coords[exon_start]), type="exon", strand=child.location.strand, qualifiers={"source": [self.source], 'ID': [child.qualifiers['Name'][0] + "-exon-" + str(count_id)], 'Name': [child.qualifiers['Name'][0] + "-exon"]})
                child.sub_features.append(new_subf)
                count_id += 1

        return cleaned_f

    # Add UTR subfeatures to a feature, guessing from CDS and exon subfeatures
    def guess_utrs(self, cleaned_f):

        for child in cleaned_f.sub_features:  # mRNA
            cds_coords = {}
            exon_coords = {}
            utr_coords = {}
            for gchild in child.sub_features:  # exons, cds, ...
                if gchild.type == "CDS":
                    cds_coords[gchild.location.start] = gchild.location.end
                elif gchild.type == "exon":
                    exon_coords[gchild.location.start] = gchild.location.end

            exon_sorted = sorted(exon_coords)
            for exon_start in exon_sorted:
                exon_end = exon_coords[exon_start]

                found_cds = False
                for cds_start in cds_coords:
                    cds_end = cds_coords[cds_start]

                    if cds_start == exon_start and cds_end == exon_end:
                        # whole exon is a CDS
                        found_cds = True

                    else:
                        if cds_start > exon_start and cds_start < exon_end:
                            # 5'UTR
                            found_cds = True
                            utr_coords[exon_start] = cds_start

                        if cds_end > exon_start and cds_end < exon_end:
                            # 3' UTR
                            found_cds = True
                            utr_coords[cds_end] = exon_end

                if not found_cds:
                    # No CDS on the exon, create a UTR on the whole exon
                    utr_coords[exon_start] = exon_end

            count_id = 1
            for utr_start in utr_coords:
                new_subf = SeqFeature(FeatureLocation(utr_start, utr_coords[utr_start]), type="UTR", strand=child.location.strand, qualifiers={"source": [self.source], 'ID': [child.qualifiers['Name'][0] + "-UTR-" + str(count_id)], 'Name': [child.qualifiers['Name'][0] + "-UTR"]})
                new_subf.sub_features = [] # See https://github.com/biopython/biopython/issues/928
                child.sub_features.append(new_subf)
                count_id += 1

        return cleaned_f

    def renumber_exons(self, cleaned_f):

        for child in cleaned_f.sub_features:  # mRNA
            exon_coords = {}
            cds_coords = {}
            utr_coords = {}
            strand = child.location.strand

            # Find positions
            for gchild in child.sub_features:  # exons, cds, ...
                if gchild.type == "exon":
                    exon_id = gchild.qualifiers['ID'][0]
                    start = gchild.location.start
                    exon_coords[start] = exon_id
                elif gchild.type == "CDS":
                    if 'ID' in gchild.qualifiers:
                        cds_id = gchild.qualifiers['ID'][0]
                    else:
                        cds_id = gchild.qualifiers['Name'][0]
                    start = gchild.location.start
                    cds_coords[start] = cds_id
                elif gchild.type in ["UTR", "five_prime_UTR", "three_prime_UTR"]:
                    if 'ID' in gchild.qualifiers:
                        utr_id = gchild.qualifiers['ID'][0]
                    else:
                        utr_id = gchild.qualifiers['Name'][0]
                    start = gchild.location.start
                    utr_coords[start] = utr_id

            # Sort by position
            if strand > 0:
                exons_sorted = sorted(exon_coords)
                cds_sorted = sorted(cds_coords)
                utr_sorted = sorted(utr_coords)
            else:
                exons_sorted = sorted(exon_coords, reverse=True)
                cds_sorted = sorted(cds_coords, reverse=True)
                utr_sorted = sorted(utr_coords, reverse=True)

            # Reassign ids
            exon_num = 1
            exon_subs = {}
            for exon_s in exons_sorted:
                exon_subs[exon_coords[exon_s]] = re.sub('-exon(-[0-9]+)?$', '-exon-' + str(exon_num), exon_coords[exon_s])
                exon_num += 1

            cds_num = 1
            cds_subs = {}
            for cds_s in cds_sorted:
                cds_subs[cds_coords[cds_s] + str(cds_s)] = re.sub('-CDS(-[0-9]+)?$', '-CDS-' + str(cds_num), cds_coords[cds_s])
                cds_num += 1

            utr_num = 1
            utr_subs = {}
            for utr_s in utr_sorted:
                utr_subs[utr_coords[utr_s] + str(utr_s)] = re.sub('-UTR(-[0-9]+)?$', '-UTR-' + str(utr_num), utr_coords[utr_s])
                utr_num += 1

            for gchild in child.sub_features:  # exons, cds, ...
                if gchild.type == "exon":
                    if gchild.qualifiers['ID'][0] in exon_subs:
                        # The gchild id can be absent from exon_subs if we are treating an exon shared between multiple isoforms,
                        # and it was already renamed while treating previous isoform
                        gchild.qualifiers['ID'][0] = exon_subs[gchild.qualifiers['ID'][0]]
                elif gchild.type == "CDS":
                    if 'ID' in gchild.qualifiers:
                        cds_id = gchild.qualifiers['ID'][0]
                    else:
                        cds_id = gchild.qualifiers['Name'][0]

                    if cds_id + str(gchild.location.start) in cds_subs:
                        # The gchild id can be absent from cds_subs if we are treating an exon shared between multiple isoforms,
                        # and it was already renamed while treating previous isoform
                        gchild.qualifiers['ID'] = [cds_subs[cds_id + str(gchild.location.start)]]
                elif gchild.type in ["UTR", "five_prime_UTR", "three_prime_UTR"]:
                    if 'ID' in gchild.qualifiers:
                        utr_id = gchild.qualifiers['ID'][0]
                    else:
                        utr_id = gchild.qualifiers['Name'][0]

                    if utr_id + str(gchild.location.start) in utr_subs:
                        # The gchild id can be absent from utr_subs if we are treating an exon shared between multiple isoforms,
                        # and it was already renamed while treating previous isoform
                        gchild.qualifiers['ID'] = [utr_subs[utr_id + str(gchild.location.start)]]

        return cleaned_f

    def match_apollo_against_base(self):

        print("Converting GFF to BED...")

        base_gff_in = open(self.filtered_base_gff, 'r')
        base_gff_out = open(self.tmpdir + '/base_cds.gff', 'w+')
        for l in base_gff_in:
            cols = l.strip().split()
            # FIXME CDS could be more appropriate (or maybe not...)
            if not l.startswith("#") and cols[2] == 'exon':
                cols[8] = re.sub(r'ID=([a-zA-Z0-9]+)', r'exID=\1', cols[8])  # remove already set id
                cols[8] = re.sub(r'Parent=([a-zA-Z0-9]+)([\.0-9]+)?([-_]R[A-Z]+)?', r'ID=\1', cols[8])  # generate a fake id based on Parent
                cols[8] = cols[8].rstrip(";")  # gff2bed doesn't like trailing ;
                print('\t'.join(cols), file=base_gff_out)
        base_gff_out.close()

        try:
            retcode = call("awk '{if ($3 == \"gene\") {print}}' FS='\t' OFS='\t' " + self.apollo_gff + " | gff2bed --do-not-sort | awk '{ if ($2 <= 0) {$2++}; print }' FS='\t' OFS='\t' | sort-bed - > " + self.tmpdir + "/apollo.bed", shell=True)
            if retcode < 0:
                print("Child was terminated by signal " + str(-retcode), file=sys.stderr)

            retcode = call("cat " + self.tmpdir + "/base_cds.gff | gff2bed > " + self.tmpdir + "/base_cds.bed", shell=True)
            if retcode < 0:
                print("Child was terminated by signal " + str(-retcode), file=sys.stderr)
        except OSError as e:
            print("Execution failed:" + e, file=sys.stderr)
            sys.exit(1)

        print("Computing intersect between new Apollo annotation and base annotation...")
        in_map = self.tmpdir + "/wa_to_all_base.tsv"
        try:
            retcode = call("intersectBed -s -wao -a " + self.tmpdir + "/apollo.bed -b " + self.tmpdir + "/base_cds.bed | cut -f4,14,21 | awk '{ if ($2 != \".\") {print} }'  | awk '{a[$1\"\t\"$2]+=$3;} END {for(i in a) print i\"\t\"a[i];}' | sort > " + in_map, shell=True)
            if retcode < 0:
                print("Child was terminated by signal " + str(-retcode), file=sys.stderr)

        except OSError as e:
            print("Execution failed:" + e, file=sys.stderr)
            sys.exit(1)

        # Load the bedtools intersect result
        in_map_h = open(in_map)

        for gene in in_map_h:
            cols = gene.strip().split("\t")
            if len(cols) < 2:
                print("Failed loading bedtools result: " + gene)
                sys.exit(1)

            cur_waid = cols[0]
            cur_gid = cols[1]
            cur_len = int(cols[2])

            if cur_waid not in self.primary_matches:
                self.primary_matches[cur_waid] = {'gid': cur_gid, 'len': cur_len}  # key is wa id, gid is base, len is exon length coverage
            else:
                # found another match for current apollo gene
                if cur_len > self.primary_matches[cur_waid]['len']:
                    # found a better match
                    if cur_waid not in self.secondary_matches:
                        self.secondary_matches[cur_waid] = []
                    self.secondary_matches[cur_waid].append(self.primary_matches[cur_waid])
                    self.primary_matches[cur_waid] = {'gid': cur_gid, 'len': cur_len}
                else:
                    # found a worse match
                    if cur_waid not in self.secondary_matches:
                        self.secondary_matches[cur_waid] = []
                    self.secondary_matches[cur_waid].append({'gid': cur_gid, 'len': cur_len})

        in_map_h.close()

    # Parse the original base annotation file
    def parse_base_annotation(self):

        agff_handle = open(self.filtered_base_gff)

        for rec in GFF.parse(agff_handle):
            rec.annotations = {}
            for f in rec.features:  # gene
                if f.type == "gene":
                    gid = f.qualifiers['ID'][0]

                    search = re.search(self.id_regex, gid)
                    if search:
                        pid_number = int(search.group(1))
                        if pid_number > self.highest_id:
                            self.highest_id = pid_number
                            self.padding_length = len(search.group(1))
                    else:
                        print("ERROR: could not parse ID from base annotation: " + gid)

                    gid = self.id_syntax.replace('{id}', search.group(1))  # Remove version suffix
                    if search.group(2):
                        self.base_genes_version[gid] = int(search.group(2).strip("."))
                    self.base_genes_structure[gid] = self.make_structure_hash(f)
                    self.base_genes_positions[gid] = {'start': f.location.start, 'end': f.location.end, 'strand': f.location.strand}

        agff_handle.close()

    # Load the previous annotation
    def parse_previous_annotation(self):

        ogff_handle = open(self.previous_gff)

        for rec in GFF.parse(ogff_handle):
            rec.annotations = {}
            for f in rec.features:  # gene
                if f.type == "gene":
                    self.num_total_genes_in_previous += 1
                    if 'Alias' in f.qualifiers:  # suspect it was already from apollo in previous annotation

                        pid = f.qualifiers['ID'][0]
                        search = re.search(self.id_regex, pid)
                        if search:
                            pid_number = int(search.group(1))
                            if pid_number > self.highest_id:
                                self.highest_id = pid_number
                        else:
                            print("ERROR: could not parse previous ID: " + pid)

                        previous_assigned_id = self.id_syntax.replace('{id}', str(pid_number).zfill(self.padding_length))

                        for alias in f.qualifiers['Alias']:
                            # Alias format is different in apollo 1.x and apollo 2.x, try to support both
                            if re.match("^[A-F0-9]{32}$", alias) or re.match("^[a-f0-9-]{36}$", alias):

                                # Getting there means that the gene was in the previous annotation and that the gene is no longer present in the new version,
                                # OR that it is there but with no bedtools result
                                # OR that it is there with the same or a different bedtools result

                                # Store the information just in case the same alias is in the new annotation
                                self.name_map[alias] = previous_assigned_id
                                self.previous_genes_content[alias] = self.make_content_hash(f, check_name=True)
                                self.previous_genes_name[alias] = pid
                                break  # we use the first alias matching a apollo id, ignoring any subsequent

        ogff_handle.close()

    # When a gene has been splitted on apollo, assign the base id to the most covered
    def handle_splitted(self):
        assigned_gids = {}
        removed_wa = []
        for wa in self.primary_matches:
            gid = self.primary_matches[wa]['gid']
            if gid in assigned_gids:
                self.num_splitted += 1
                if assigned_gids[gid]['len'] > self.primary_matches[wa]['len']:
                    # a previous apollo gene is more covered by the base gene than the current one
                    # move the current match from the primary to the secondary list
                    if wa not in self.secondary_matches:
                        self.secondary_matches[wa] = []
                    self.secondary_matches[wa].append(self.primary_matches[wa])
                    removed_wa.append(wa)
                else:
                    # a new best match
                    # move the previous best from the primary to the secondary list
                    if assigned_gids[gid]['waid'] not in self.secondary_matches:
                        self.secondary_matches[assigned_gids[gid]['waid']] = []
                    self.secondary_matches[assigned_gids[gid]['waid']].append(self.primary_matches[assigned_gids[gid]['waid']])
                    removed_wa.append(assigned_gids[gid]['waid'])

                    # keep current in the primary list
                    assigned_gids[gid] = {'waid': wa, 'len': self.primary_matches[wa]['len']}
            else:
                assigned_gids[gid] = {'waid': wa, 'len': self.primary_matches[wa]['len']}

        for wa in removed_wa:
            del self.primary_matches[wa]

    # When an apollo gene covers several base genes completely, make sure we remove all 100% covered base genes
    # => 3 genes completely covered by 1 apollo gene wll give just 1 gene in output
    def handle_merged(self):

        for wa in self.primary_matches:
            if wa in self.secondary_matches:
                apollo_pos = self.apollo_genes_positions[wa]
                for secs in self.secondary_matches[wa]:
                    base_pos = self.base_genes_positions[secs['gid']]

                    search = re.search(self.id_regex, secs['gid'])
                    raw_merged_gene_id = self.id_syntax.replace('{id}', search.group(1))

                    # Check that the gene is not already blacklisted (ie, the primary match of another gene?)
                    if raw_merged_gene_id not in self.blacklist_base:
                        if base_pos['strand'] == apollo_pos['strand'] and (base_pos['start'] >= apollo_pos['start'] and base_pos['start'] <= apollo_pos['end']) or (base_pos['end'] >= apollo_pos['start'] and base_pos['end'] <= apollo_pos['end']):
                            # base gene is somewhat included in the apollo gene, delete it
                            # we could trust secondary_matches and don't check the position, but you never know...
                            if wa not in self.blacklist_merged:
                                self.blacklist_merged[wa] = []
                            self.blacklist_merged[wa].append(raw_merged_gene_id)
                        else:
                            print('ERROR: Strange merge case for ' + str(wa) + ' -> ' + str(raw_merged_gene_id))

    def check_id_compatibility(self):
        for wa in self.primary_matches:

            if wa in self.name_map:
                if self.primary_matches[wa]['gid'] != self.name_map[wa]:
                    self.kept_compatibility += 1
                    print("WARNING: In the previous annotation, gene '" + wa + "' was assigned id '" + self.name_map[wa] + "' but in the new one it would be assigned id '" + self.primary_matches[wa]['gid'] + "'. Keeping '" + self.name_map[wa] + "' for compatibility.")
                # else: nothing to change, the previous id and the one found with bedtools is the same
            else:
                if self.primary_matches[wa]['gid'] not in self.name_map.values():
                    self.name_map[wa] = self.primary_matches[wa]['gid']
                else:
                    already_assigned = ""
                    for w, g in self.name_map.items():
                        if g == self.primary_matches[wa]['gid']:
                            already_assigned = w
                    # The id was already used for another gene, don't store any id mapping for this gene, we will generate a new id later
                    # This happens when a gene was splitted by annotators (but not only)
                    print("WARNING: Gene '" + wa + "' should be assigned id '" + self.primary_matches[wa]['gid'] + "' but it is already used by gene '" + already_assigned + "'. A new id will be created.")

    def parse_apollo_annotation(self):
        # Load the new WA annotation
        wgff_handle = open(self.apollo_gff)
        for rec in GFF.parse(wgff_handle):
            rec.annotations = {}
            for f in rec.features:  # gene
                if 'ID' in f.qualifiers and f.type == "gene":
                    wa_id = f.qualifiers['ID'][0]
                    self.apollo_ids_in_latest.append(wa_id)

                    base_id = ""

                    if wa_id not in self.name_map:
                        # A new gene that have no bedtools result and that wasn't found in the previous annotation version
                        # Give it a completely new id
                        base_id = self.id_syntax.replace('{id}', str(self.highest_id).zfill(self.padding_length))
                        self.name_map[wa_id] = base_id + ".1"
                        self.highest_id += 1
                        self.num_new_id += 1
                    elif wa_id in self.previous_genes_content:
                        # Found in the previous and the new annotation => we have to compare the previous and new version to choose the good suffix
                        f = self.guess_utrs(f)  # apollo doesn't generate UTRs lines in the gff output, guess them before comparing
                        new_hash = self.make_content_hash(f, check_name=True)
                        search = re.search(self.id_regex, self.previous_genes_name[wa_id])
                        if self.previous_genes_content[wa_id] == new_hash:
                            # Gene was the same in the previous annotation, keeping its previous ID
                            self.num_unchanged_from_previous += 1
                            self.name_map[wa_id] = self.previous_genes_name[wa_id]
                            if search:
                                base_id = self.id_syntax.replace("{id}", search.group(1))
                            else:
                                print("ERROR: could not parse previous ID: " + self.previous_genes_name[wa_id])
                        else:
                            # Gene changed, updating its suffix
                            if search:
                                # The previous gene is sure to have a version suffix as it was already incremented when creating the previous OGS
                                suffix = int(search.group(2).strip(".")) + 1
                            else:
                                print("ERROR: could not parse previous ID: " + self.previous_genes_name[wa_id])
                            base_id = self.name_map[wa_id]
                            self.name_map[wa_id] = self.name_map[wa_id] + "." + str(suffix)
                            self.num_version_bumped += 1
                    else:
                        # The gene was not seen in previous annotation and we have a bedtools result (=found in base annotation)
                        base_id = self.name_map[wa_id]
                        search = re.search(self.id_regex, self.name_map[wa_id])
                        if base_id not in self.base_genes_version:
                            # No version suffix in base annotation, add one
                            self.name_map[wa_id] = self.name_map[wa_id] + ".1"
                        else:
                            # There was already a version suffix, increment it
                            suffix = self.base_genes_version[base_id] + 1
                            self.name_map[wa_id] = self.id_syntax.replace('{id}', search.group(1)) + "." + str(suffix)
                        self.num_new_id += 1

                    # See if the gene had exactly the same structure in the base annotation
                    if base_id in self.base_genes_structure and self.base_genes_structure[base_id] == self.make_structure_hash(f):
                        self.num_unchanged_from_base += 1

                    self.apollo_genes_positions[wa_id] = {'start': f.location.start, 'end': f.location.end, 'strand': f.location.strand}

        wgff_handle.close()

    # Generate a list of genes to delete from the base annotation
    def prepare_base_blacklist(self):
        # It also contains ids that were generated for new genes (will not be found in base annotation file)
        self.blacklist_base = []
        for waid in self.name_map:
            # Only keep in the blacklist the genes that are still in the latest apollo annotation
            # We do this to avoid removing a gene that was replaced in the previous annotation, but is no longer replaced in the new
            # The version suffix is not present in blastlist_base
            if waid in self.apollo_ids_in_latest:
                search = re.search(self.id_regex, self.name_map[waid])
                if search:
                    self.blacklist_base.append(self.id_syntax.replace('{id}', search.group(1)))
                else:
                    print("ERROR: could not parse name mapping ID: " + self.name_map[waid])

    # Remove all isoforms marked as deleted by annotators
    # Write a filtered gff into the tmpdir
    def filter_deleted_isoforms(self):

        print("Filtering deleted isoforms...")

        mrnas_to_delete = []

        # Check the mRNAs that were specifically deleted by annotators
        if self.deleted:
            with open(self.deleted, 'rU') as f:
                for line in f:
                    mrnas_to_delete.append(line.strip())

        in_handle = open(self.base_gff)

        recs = []
        for rec in GFF.parse(in_handle):
            rec.annotations = {}
            rec.seq = ""
            new_feats = []

            for gene in rec.features:  # gene

                new_subs = []

                for mRNA in gene.sub_features:
                    if mRNA.qualifiers['ID'][0] not in mrnas_to_delete:
                        new_subs += [mRNA]
                    else:
                        self.num_deleted_isoforms += 1

                if len(new_subs) == 0:
                    # No more isoforms in the gene, delete the whole gene
                    self.whole_genes_deleted.append(gene.qualifiers['ID'][0])
                    continue

                gene.sub_features = new_subs
                new_feats.append(gene)

            rec.features = new_feats

            if len(rec.features):
                recs.append(rec)

        in_handle.close()

        out_handle = open(self.filtered_base_gff, "w")
        GFF.write(recs, out_handle)
        out_handle.close()

    # Write final GFF
    def write_gff(self):
        # First write genes kept from the base annotation...
        in_handle = open(self.filtered_base_gff)
        recs = []

        all_merged = [j for i in self.blacklist_merged.values() for j in i]
        self.num_merged = len(all_merged)

        for rec in GFF.parse(in_handle):
            rec.annotations = {}
            rec.seq = ""
            new_feats = []

            for f in rec.features:  # gene

                gene_id = f.qualifiers['ID'][0]

                search = re.search(self.id_regex, gene_id)
                raw_gene_id = self.id_syntax.replace('{id}', search.group(1))

                if raw_gene_id not in self.blacklist_base and raw_gene_id not in all_merged:
                    cleaned_f = self.clean_feature(f)
                    cleaned_f = self.guess_exons(cleaned_f)
                    cleaned_f = self.renumber_exons(cleaned_f)
                    new_feats.append(cleaned_f)
                    self.num_kept_base += 1
                    self.num_total_genes += 1
                else:
                    self.num_replaced_base += 1

            rec.features = new_feats

            if len(rec.features):
                recs.append(rec)

        in_handle.close()

        # ... then write genes models from apollo
        in_handle = open(self.apollo_gff)
        for rec in GFF.parse(in_handle):
            rec.annotations = {}
            rec.seq = ""
            new_feats = []

            for f in rec.features:  # gene

                gene_id = f.qualifiers['ID'][0]

                if gene_id in self.name_map:
                    f.qualifiers['Alias'] = [gene_id]
                    f.qualifiers['ID'][0] = self.name_map[gene_id]
                    nfeat = self.clean_feature(f)
                    nfeat = self.guess_utrs(nfeat)
                    nfeat = self.renumber_exons(nfeat)
                    new_feats.append(nfeat)
                    self.num_total_genes += 1

            rec.features = new_feats

            if len(rec.features):
                recs.append(rec)

        in_handle.close()

        out_handle = open(self.out_gff, "w")
        GFF.write(recs, out_handle)
        out_handle.close()

    def print_statistics(self, file):
        num_apollo_genes = self.num_new_id + self.num_unchanged_from_previous + self.num_version_bumped

        print("Starting numbering of new genes from: {}".format(self.highest_id), file=file)

        print("    {} genes in the base annotation".format(len(self.base_genes_structure) + len(self.whole_genes_deleted)), file=file)
        if self.previous_gff:
            print("    {} genes in the previous annotation, including {} Apollo genes".format(self.num_total_genes_in_previous, len(self.previous_genes_name)), file=file)
        print("    {} genes in the new Apollo annotation".format(num_apollo_genes), file=file)
        print("", file=file)
        print("    {} genes in the new final gff file".format(self.num_total_genes), file=file)
        print("        {} genes from the base annotation were kept untouched".format(self.num_kept_base), file=file)
        print("        {} genes from apollo have no equivalent in base annotation".format(self.num_total_genes - len(self.base_genes_structure) + self.num_merged), file=file)
        print("        {} genes from the base annotation were replaced by an apollo gene".format(self.num_replaced_base), file=file)
        print("            {} base genes were merged into {} Apollo genes".format(self.num_merged, len(self.blacklist_merged)), file=file)
        print("        {} isoforms ({} whole genes) were deleted from the base annotation (by annotator request)".format(self.num_deleted_isoforms, len(self.whole_genes_deleted)), file=file)
        print("", file=file)
        if self.previous_gff:
            print("    {} ids compared to previous annotation: kept {} unchanged, updated {} ids, and {} are new.".format(num_apollo_genes, self.num_unchanged_from_previous, self.num_version_bumped, self.num_new_id), file=file)
            print("    {} genes have an id that was kept from previous annotation for compatibility reason".format(self.kept_compatibility), file=file)
            print("", file=file)
        print("    {} apollo genes have the exact same structure as in the base annotation".format(self.num_unchanged_from_base), file=file)
        print("    {} genes were splitted by annotators".format(self.num_splitted), file=file)
        print("", file=file)

    def print_convertion(self, file):

        print("# Id conversion table from {} to {}".format(os.path.basename(self.base_gff), self.source), file=file)
        if self.previous_gff:
            print("# {}\t{}\t{}\tapollo_gene_id".format(self.source, os.path.basename(self.base_gff), os.path.basename(self.base_gff)), file=file)
        else:
            print("# {}\t{}\tapollo_gene_id".format(self.source, os.path.basename(self.base_gff)), file=file)

        for waid in self.name_map:
            old_id = '-'
            search = re.search(self.id_regex, self.name_map[waid])
            if search:
                old_id = self.id_syntax.replace('{id}', search.group(1))
                if old_id not in self.base_genes_structure:
                    # Ahah, the id in name_map is a new one
                    old_id = '-'

            if self.previous_gff:
                previous_id = '-'
                if waid in self.previous_genes_name:
                    previous_id = self.previous_genes_name[waid]

                print("{}\t{}\t{}\t{}".format(self.name_map[waid], old_id, previous_id, waid), file=file)
            else:
                print("{}\t{}\t{}".format(self.name_map[waid], old_id, waid), file=file)

        for d in self.whole_genes_deleted:
            if self.previous_gff:
                print("deleted\t{}\t-\t-".format(d), file=file)
            else:
                print("deleted\t{}\t-".format(d), file=file)

    # Create fasta output using gffread
    def write_fasta(self):
        print("Generating fasta files using gffread...")
        try:
            retcode = call("gffread " + self.out_gff + " -g " + self.genome + " -w " + self.out_transcript + " -x " + self.out_cds + " -y " + self.tmpdir + '/proteins.fa', shell=True)
            if retcode < 0:
                print("Child was terminated by signal " + str(-retcode), file=sys.stderr)
        except OSError as e:
            print("Execution failed:" + e, file=sys.stderr)
            sys.exit(1)

        # Protein fasta file need to have modified id
        prot_in = open(self.tmpdir + '/proteins.fa', 'r')
        prot_out = open(self.out_protein, 'w+')
        for prot in prot_in:
            prot = prot.strip()
            if prot.startswith(">"):
                prot = re.sub(r'-R([A-Z]+)', r'-P\1', prot)
            print(prot, file=prot_out)
        prot_out.close()

    def parse_args(self):

        parser = argparse.ArgumentParser()
        parser.add_argument("genome", help="Genome file (fasta)")
        parser.add_argument("ogs_name", help="Name of the new OGS")
        parser.add_argument("id_regex", help="Regex with a capturing group around the incremental part of gene ids, and a second one around the version suffix (e.g. 'GSSPF[GPT]([0-9]{8})[0-9]{3}(\.[0-9]+)?')")
        parser.add_argument("id_syntax", help="String representing a gene id, with {id} where the incremental part of the id should be placed (e.g. 'GSSPFG{id}001')")
        parser.add_argument("base_gff", help="The gff from the base annotation (usually automatic annotation)")
        parser.add_argument("apollo_gff", help="The gff from the new Apollo valid annotation")
        # In --previous_gff, the genes should have an Alias attribute containing the Apollo ID
        parser.add_argument("-p", "--previous_gff", help="The gff from the previous annotation version (if different than <base_gff>)")
        parser.add_argument("-d", "--deleted", help="File containing a list of mRNAs to remove")
        parser.add_argument("-o", "--out_prefix", help="Prefix for output files (default=<ogs_name>_<today's date>)")
        args = parser.parse_args()

        self.base_gff = args.base_gff
        self.filtered_base_gff = self.tmpdir + '/filtered_base.gff'  # base_gff without the deleted isoforms/genes
        self.previous_gff = args.previous_gff
        self.apollo_gff = args.apollo_gff
        self.genome = args.genome
        self.source = args.ogs_name
        self.deleted = args.deleted

        self.id_regex = args.id_regex
        self.id_syntax = args.id_syntax

        self.out_prefix = args.out_prefix
        if not self.out_prefix:
            self.out_prefix = self.source + '_' + datetime.date.today().strftime('%Y%m%d')

        self.out_gff = self.out_prefix + '.gff3'
        self.out_transcript = self.out_prefix + '_transcript.fasta'
        self.out_cds = self.out_prefix + '_cds.fasta'
        self.out_protein = self.out_prefix + '_protein.fasta'
        self.out_report = self.out_prefix + '_merge_report.txt'

        self.highest_id = 0
        self.padding_length = 0

        # Intersect results
        self.primary_matches = {}
        self.secondary_matches = {}

        # Information from base annotation
        self.base_genes_structure = {}
        self.base_genes_version = {}  # Store the original version of the genes
        self.base_genes_positions = {}  # Store positions of all base genes

        # Information from previous annotation
        self.previous_genes_content = {}
        self.previous_genes_name = {}
        self.num_total_genes_in_previous = 0
        self.kept_compatibility = 0  # Number of genes which id was kept from previous annotation for compatibility

        # Information from apollo annotation
        self.apollo_ids_in_latest = []
        self.apollo_genes_positions = {}  # Store positions of all apollo genes

        # Name mapping
        self.name_map = {}  # keys = Apollo id, values = best corresponding base annotation id
        self.blacklist_merged = {}

        # Statistics
        self.num_unchanged_from_base = 0
        self.num_unchanged_from_previous = 0
        self.num_version_bumped = 0
        self.num_new_id = 0
        self.num_kept_base = 0
        self.num_total_genes = 0
        self.num_replaced_base = 0
        self.num_splitted = 0
        self.num_merged = 0
        self.num_deleted_isoforms = 0
        self.whole_genes_deleted = []

    def main(self):

        # Create a temp directory to launch bedtools
        self.tmpdir = tempfile.mkdtemp()

        self.parse_args()

        # Delete isoforms (when asked by annotators)
        self.filter_deleted_isoforms()

        # Launch external tools to compute intersection between apollo genes and genes from base annotation
        self.match_apollo_against_base()

        # Parse the base annotation file
        self.parse_base_annotation()

        # Load the previous annotation (if any)
        if self.previous_gff:
            self.parse_previous_annotation()

        # Assign IDs to genes that got splitted in apollo
        self.handle_splitted()

        # Ensure we assign the same ids as in previous annotation
        self.check_id_compatibility()

        # Increase a little to be sure we don't reuse the id of a gene that was removed
        self.highest_id += 100

        # Parse latest apollo annotation
        self.parse_apollo_annotation()

        # Blacklist genes that were replaced
        self.prepare_base_blacklist()

        # Find genes that got merged in apollo
        self.handle_merged()

        # Print output gff
        self.write_gff()

        # Print output gff
        self.write_fasta()

        # Print some general statistics
        self.print_statistics(sys.stdout)
        self.print_statistics(open(self.out_prefix + '_merge_stats.txt', 'w'))

        # Print ID convertion table from base annotation to new one
        self.print_convertion(open(self.out_prefix + '_merge_ids_conversion.tsv', 'w'))

        # Remove temp dir
        shutil.rmtree(self.tmpdir)


if __name__ == '__main__':
    ogsmerger = OgsMerger()

    ogsmerger.main()
