from Bio import SeqFeature
from Bio.SeqIO.InsdcIO import EmblWriter
from Bio.SeqIO.InsdcIO import _insdc_location_string_ignoring_strand_and_subfeatures

"""
This version of the EMBL writer will write a template EMBL file ready for submission to EBI.
Some fields are intentionnally left empty as they are automatically filled by EBI affter submission.
"""


def _insdc_location_string(location, rec_length, deep=False):
    """Build a GenBank/EMBL location from a (Compound) FeatureLocation (PRIVATE).

    There is a choice of how to show joins on the reverse complement strand,
    GenBank used "complement(join(1,10),(20,100))" while EMBL used to use
    "join(complement(20,100),complement(1,10))" instead (but appears to have
    now adopted the GenBank convention). Notice that the order of the entries
    is reversed! This function therefore uses the first form. In this situation
    we expect the CompoundFeatureLocation and its parts to all be marked as
    strand == -1, and to be in the order 19:100 then 0:10.
    """
    try:
        parts = location.parts
        # CompoundFeatureLocation
        if location.strand == -1:
            # Special case, put complement outside the join/order/... and reverse order
            return "complement(%s(%s))" % (location.operator,
                                           ",".join(_insdc_location_string(p, rec_length, True) for p in parts))
        else:
            return "%s(%s)" % (location.operator,
                               ",".join(_insdc_location_string(p, rec_length, True) for p in parts))
    except AttributeError:
        # Simple FeatureLocation
        loc = _insdc_location_string_ignoring_strand_and_subfeatures(location, rec_length)
        if (location.strand == -1) and (not deep):
            return "complement(%s)" % loc
        else:
            return loc


class BipaaEmblSubmitWriter(EmblWriter):
    def _write_feature(self, feature, record_length):
        """Write a single SeqFeature object to features table."""
        assert feature.type, feature
        location = _insdc_location_string(feature.location, record_length)
        f_type = feature.type.replace(" ", "_")
        line = (self.QUALIFIER_INDENT_TMP % f_type)[:self.QUALIFIER_INDENT] + self._wrap_location(location) + "\n"
        self.handle.write(line)
        # Now the qualifiers...
        for key in sorted(feature.qualifiers.keys()):
            values = feature.qualifiers[key]
            if isinstance(values, list) or isinstance(values, tuple):
                for value in values:
                    self._write_feature_qualifier(key, value)
            else:
                # String, int, etc - or None for a /pseudo tpy entry
                self._write_feature_qualifier(key, values)

    def _write_the_first_lines(self, record):
        """Write the ID and AC lines."""
        if "." in record.id and record.id.rsplit(".", 1)[1].isdigit():
            accession = self._get_annotation_str(record, "accession",
                                                 record.id.rsplit(".", 1)[0],
                                                 just_first=True)
        else:
            accession = self._get_annotation_str(record, "accession",
                                                 record.id,
                                                 just_first=True)

        if ";" in accession:
            raise ValueError("Cannot have semi-colon in EMBL accession, %s"
                             % repr(str(accession)))
        if " " in accession:
            # This is out of practicallity... might it be allowed?
            raise ValueError("Cannot have spaces in EMBL accession, %s"
                             % repr(str(accession)))

        # TODO - Full ID line
        handle = self.handle
        # ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
        # 1. Primary accession number
        # 2. Sequence version number
        # 3. Topology: 'circular' or 'linear'
        # 4. Molecule type
        # 5. Data class
        # 6. Taxonomic division
        # 7. Sequence length
        self._write_single_line("ID", "%s; %s; %s; %s; XXX; %s; %s %s."
                                % ('XXX', 'SV XXX', 'linear', 'genomic DNA',
                                   'XXX', 'XXX', 'XXX'))
        handle.write("XX\n")
        self.handle.write("AC   ;\n")
        handle.write("XX\n")
        self.handle.write("AC * _%s\n" % (accession))
        handle.write("XX\n")

    def _write_references(self, record):
        # The order should be RN, RC, RP, RX, RG, RA, RT, RL
        number = 0
        for ref in record.annotations["references"]:
            if not isinstance(ref, SeqFeature.Reference):
                continue
            number += 1
            self._write_single_line("RN", "[%i]" % number)
            # TODO - support for RC line (needed in parser too)
            # TODO - support more complex record reference locations?
            if ref.location and len(ref.location) == 1:
                self._write_single_line(
                    "RP", "%i-%i" % (ref.location[0].nofuzzy_start + 1,
                                     ref.location[0].nofuzzy_end))
            # TODO - record any DOI or AGRICOLA identifier in the reference object?
            if ref.pubmed_id:
                self._write_single_line("RX", "PUBMED; %s." % ref.pubmed_id)
            if ref.consrtm:
                self._write_single_line("RG", "%s" % ref.consrtm)
            if ref.authors:
                # We store the AUTHORS data as a single string
                self._write_multi_line("RA", "XXX;")
            if ref.title and ref.title != " ":
                # We store the title as a single string
                self._write_multi_line("RT", '"%s";' % ref.title)
            else:
                self._write_multi_line("RT", ';')
            if ref.journal:
                # We store this as a single string - holds the journal name,
                # volume, year, and page numbers of the citation
                self._write_multi_line("RL", ref.journal)
            self.handle.write("XX\n")
