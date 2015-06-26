from re import split
from warnings import warn
from collections import defaultdict
from HTSeq import GenomicInterval
from Bio.Seq import Seq

class WrongCDBoxOrderExeption(Exception): pass
class ToManyBoxesException(Exception): pass
class IncompatibleStrandAndCoordsException(Exception): pass
class WrongCDBoxPlacementException(Exception): pass
class NoBoxException(Exception): pass
class WrongBoxException(Exception): pass
class InteractionRegionTooShortException(Exception): pass

class snoRNA(object):

    """A class representing snoRNA"""

    def __init__(self, snor_id, organism, chrom, start, end, strand,
                 sequence, snor_type, **kwargs):
        """@todo: to be defined1.

        Args:
            snor_id (str): a unique id for snoRNA
            organism (str): species in which snoRNA can be found
            chrom (str): @todo
            start (int): @todo
            end (int): @todo
            strand (str): @todo
            sequence (str): @todo
            snor_type (str): @todo

        Kwargs:
            alias (str): anlternative name for snoRNA
            gene_name (str): snoRNA gene name
            accession (str): accession for the snoRNA (eg. NCBI)
            modified_sites (str): string of modified sites eg. 28S:U46,18S:G52 separated
                                  by the coma.
                                  It will be transformed to dictionary of the form
                                  {rna: [(position, nucleotide)]}
            host_gene (str): host gene
            host_id (str): id for host locus
            organization (str): organization of the locus
            note (str): additional information about snoRNA


        """
        #
        # args
        #
        self.snor_id = snor_id
        self.organism = organism
        self.position = GenomicInterval(chrom, start, end, strand)
        self.sequence = Seq(sequence.upper())
        self.snor_type = snor_type
        #
        # kwargs
        #
        self.alias = kwargs.get("alias", None)
        self.gene_name = kwargs.get("gene_name", None)
        self.accession = kwargs.get("accession", None)
        self.__assign_modification_sites(kwargs.get("modified_sites", None))
        self.host_gene = kwargs.get("host_gene", None)
        self.host_id = kwargs.get("host_id", None)
        self.organization = kwargs.get("organization", None)
        self.note = kwargs.get("note", None)
        self.__validate()


    def __validate(self):
        """Validate if CD box snoRNA is correct and can be
        used in calculations"""
        if self.position.length != len(self.sequence):
            m = "Length of sequence is different than length deduced from coordinates of %s: " % self.snor_id
            m += "coords - %i vs sequence - %i. Sequence: \n" % (self.position.length, len(self.sequence))
            m += "%s" % (self.sequence)
            # raise IncompatibleStrandAndCoordsException(m)
            warn(m, Warning)

    def __assign_modification_sites(self, sites):
        if sites is not None:
            try:
                modsites = filter(None, split(r",\s*", sites))
            except Exception, e:
                warn("Exception raised in __assign_modification_sites: %s" % str(e), Warning)
                modsites = []
            if len(modsites) == 0:
                self.modified_sites = None
            else:
                self.modified_sites = defaultdict(list)
                for count, site in enumerate(modsites):
                    try:
                        rna, pos = site.split(":")
                        self.modified_sites[rna].append((int(pos[1:]), pos[0]))
                    except Exception, e:
                        warn("Exception raised in __assign_modification_sites: %s, %s" % (str(e),
                                                                                          str(modsites)), Warning)
                        if count + 1 == len(modsites):
                            self.modified_sites = None
                            break
        else:
            self.modified_sites = None


    def get_bed_string(self):
        bed_str = "%s\t%s\t%s\tsnoRNA:%s\t%s\t%s\n" % (self.position.chrom,
                                                self.position.start,
                                                self.position.end,
                                                self.snor_id,
                                                self.alias,
                                                self.position.strand)
        return bed_str

    def get_fasta_string(self):
        return ">%s\n%s\n" % (self.snor_id, self.sequence)

    def is_orphan(self):
        if not self.modified_sites:
            return True
        else:
            return False

    def get_gff_ensembl_string(self):
        attr_tup = (self.gene_name,
                    self.snor_id)
        chrom = self.position.chrom
        attr = 'gene_id "%s"; gene_version "1"; gene_name "%s"; gene_source "snOPY"; gene_biotype "snoRNA";' % attr_tup
        gff_str = "%s\tsnOPY\tgene\t%i\t%i\t.\t%s\t.\t%s\n" % (chrom[3:] if chrom.startswith("chr") else chrom,
                                                               self.position.start,
                                                               self.position.end,
                                                               self.position.strand,
                                                               attr)
        return gff_str

    def __repr__(self):
        return "<%s %s>" % (self.snor_type, self.snor_id)

class CD_snoRNA(snoRNA):

    """A class representing C/D box snoRNA"""

    def __init__(self, snor_id, organism, chrom, start, end, strand,
                 sequence, snor_type, d_boxes, **kwargs):
        """@todo: to be defined1.

        Args:
            snor_id (str): a unique id for snoRNA
            organism (str): species in which snoRNA can be found
            chrom (str): @todo
            start (int): @todo
            end (int): @todo
            strand (str): @todo
            sequence (str): @todo
            snor_type (str): @todo
            d_boxes (str): coma separated list of d boxes positions: "1..2,5..6".
                             Start and end should be separated by .. and maximum length is 2.

        Kwargs:
            c_boxes (str): coma separated list of d boxes positions: "1..2,5..6".
                           Start and end should be separated by .. and maximum length is 2.
            switch_boxes (bool): if the boxes are in wrong order switch the order, default False
            alias (str): anlternative name for snoRNA
            gene_name (str): snoRNA gene name
            accession (str): accession for the snoRNA (eg. NCBI)
            modified_sites (dict): dictionary of modified sites: {rna: {position, nucleotide}}
            host_gene (str): host gene
            host_id (str): id for host locus
            organization (str): organization of the locus
            note (str): additional information about snoRNA


        """
        snoRNA.__init__(self, snor_id, organism, chrom, start, end, strand, sequence, snor_type, **kwargs)
        self.switch_boxes = kwargs.get("switch_boxes", False)
        self.__asign_d_boxes(d_boxes)
        self.__asign_c_boxes(kwargs.get("c_boxes", None))
        self.__validate_cd()


    def __asign_d_boxes(self, d_boxes):
        """
        Assign D-boxes to snoRNA based on the string provided. D-boxes should be in the
        form of "D-box_start..D-box_end" or "D'-box_start..D'-box_end,D-box_start..D-box_end"
        """
        if d_boxes != None:
            try:
                boxes = filter(None, split(r",\s*", d_boxes))
            except TypeError:
                boxes = []
            if len(boxes) == 0 or d_boxes == None:
                warn("No D-boxes were assigned for %s" % self.snor_id, Warning)
                self.d_box = None
                self.dprime_box = None
            elif len(boxes) == 1:
                dstart, dend = [int(i) for i in boxes[0].split("..")]
                if dstart < (len(self.sequence)/2):
                    raise WrongCDBoxPlacementException("D-box in first half of the snoRNA sequence for %s. It looks like D' box!" % self.snor_id)
                self.d_box = (dstart, dend)
                self.dprime_box = None
            elif len(boxes) == 2:
                dpstart, dpend = [int(i) for i in boxes[0].split("..")]
                dstart, dend = [int(i) for i in boxes[1].split("..")]
                if not self.switch_boxes:
                    if dpstart > dstart:
                        raise WrongCDBoxOrderExeption("D'-box start is higher than D-box start for %s" % self.snor_id)
                    else:
                        self.d_box = (dstart, dend)
                        self.dprime_box = (dpstart, dpend)
                else:
                    if dpstart > dstart:
                        warn("Switching D boxes for %s" % self.snor_id)
                        self.d_box = (dpstart, dpend)
                        self.dprime_box = (dstart, dend)
                    else:
                        self.d_box = (dstart, dend)
                        self.dprime_box = (dpstart, dpend)
            else:
                raise ToManyBoxesException("To many D-boxes for %s" % self.snor_id)
        else:
            raise NoBoxException("No D-box found in this snoRNA!")

    def __asign_c_boxes(self, c_boxes):
        """
        Assign C-boxes to snoRNA based on the string provided. C-boxes should be in the
        form of "C-box_start..C-box_end" or "C'-box_start..C'-box_end,C-box_start..C-box_end"
        """
        if c_boxes != None:
            try:
                boxes = filter(None, split(r",\s*", c_boxes))
            except TypeError:
                boxes = []
            if len(boxes) == 0:
                warn("No C-boxes were assigned for %s" % self.snor_id)
                self.c_box = None
                self.cprime_box = None
            elif len(boxes) == 1:
                cstart, cend = [int(i) for i in boxes[0].split("..")]
                if cstart > (len(self.sequence)/2):
                    warn("C-box in the second half of the snoRNA sequence for %s. It looks like C' box so I am not assigning any C box!" % self.snor_id,
                         Warning)
                self.c_box = (cstart, cend)
                self.cprime_box = None
            elif len(boxes) == 2:
                cpstart, cpend = [int(i) for i in boxes[0].split("..")]
                cstart, cend = [int(i) for i in boxes[1].split("..")]
                if not self.switch_boxes:
                    if cpstart < cstart:
                        raise WrongCDBoxOrderExeption("C'-box start is lower than C-box start for %s" % self.snor_id)
                    else:
                        self.c_box = (cstart, cend)
                        self.cprime_box = (cpstart, cpend)
                else:
                    if cpstart < cstart:
                        warn("Switching C boxes for %s" % self.snor_id)
                        self.c_box = (cpstart, cpend)
                        self.cprime_box = (cstart, cend)
                    else:
                        self.c_box = (cstart, cend)
                        self.cprime_box = (cpstart, cpend)
            else:
                raise ToManyBoxesException("To many C-boxes for %s" % self.snor_id)
        else:
            warn("No C-boxes were assigned for %s" % self.snor_id)
            self.c_box = None
            self.cprime_box = None

    def __validate_cd(self):
        """Validate if CD box snoRNA is correct and can be
        used in calculations"""
        if self.d_box:
            if (self.d_box[1] > len(self.sequence)) or (self.d_box[0] < 0):
                raise WrongBoxException("Box outside snoRNA")
            if len(self.get_dbox_sequence()) != 4:
                raise WrongBoxException("D-box should be four nucleotide long")

    def get_dbox_sequence(self):
        if self.d_box:
            return self.sequence[self.d_box[0] - 1: self.d_box[1]]
        else:
            raise NoBoxException("No D-box in this snoRNA")

    def get_cbox_sequence(self):
        if self.c_box:
            return self.sequence[self.c_box[0] - 1: self.c_box[1]]
        else:
            return None

    def get_dprime_box_sequence(self):
        if self.dprime_box:
            return self.sequence[self.dprime_box[0] - 1: self.dprime_box[1]]
        else:
            return None

    def get_cprime_box_sequence(self):
        if self.cprime_box:
            return self.sequence[self.cprime_box[0] - 1: self.cprime_box[1]]
        else:
            return None

    def __get_box_types(self):
        if self.d_box and not self.dprime_box and not self.c_box and not self.cprime_box:
            return "d"
        elif self.d_box and self.dprime_box and not self.c_box and not self.cprime_box:
            return "d_dp"
        elif self.d_box and not self.dprime_box and self.c_box and not self.cprime_box:
            return "d_c"
        elif self.d_box and not self.dprime_box and not self.c_box and self.cprime_box:
            return "d_cp"
        elif self.d_box and self.dprime_box and self.c_box and not self.cprime_box:
            return "d_dp_c"
        elif self.d_box and not self.dprime_box and self.c_box and self.cprime_box:
            return "d_c_cp"
        elif self.d_box and self.dprime_box and self.c_box and self.cprime_box:
            return "d_dp_c_cp"
        else:
            raise Exception("Wrong composition of boxes for %s" % self.snor_id)

    def get_plexy_string(self):
        if self.__get_box_types() == "d":
            plexy = ">%s-_NNNNNNN_%i_%s_%i\n%s\n" % (self.snor_id,
                                                     10,
                                                     self.get_dbox_sequence(),
                                                     self.d_box[0],
                                                     self.sequence)

            return plexy
        elif self.__get_box_types() == "d_dp":
            plexy = ">%s-_NNNNNNN_%i_%s_%i_NNNNNNN_%i_%s_%i\n%s\n" % (self.snor_id,
                                                                      self.d_box[0] - 10,
                                                                      self.get_dprime_box_sequence(),
                                                                      self.dprime_box[0],
                                                                      10,
                                                                      self.get_dbox_sequence(),
                                                                      self.d_box[0],
                                                                      self.sequence)

            return plexy
        elif self.__get_box_types() == "d_c" or self.__get_box_types() == "d_c_cp":
            plexy = ">%s-_%s_%i_%s_%i\n%s\n" % (self.snor_id,
                                                self.get_cbox_sequence(),
                                                self.c_box[0],
                                                self.get_dbox_sequence(),
                                                self.d_box[0],
                                                self.sequence)
            return plexy
        elif self.__get_box_types() == "d_dp_c":
            plexy = ">%s-_NNNNNNN_%i_%s_%i_%s_%i_%s_%i\n%s\n" % (self.snor_id,
                                                                 self.d_box[0] - 10,
                                                                 self.get_dprime_box_sequence(),
                                                                 self.dprime_box[0],
                                                                 self.get_cbox_sequence(),
                                                                 self.c_box[0],
                                                                 self.get_dbox_sequence(),
                                                                 self.d_box[0],
                                                                 self.sequence)

            return plexy
        elif self.__get_box_types() == "d_dp_c_cp":
            plexy = ">%s-_%s_%i_%s_%i_%s_%i_%s_%i\n%s\n" % (self.snor_id,
                                                                 self.get_cprime_box_sequence(),
                                                                 self.cprime_box[0],
                                                                 self.get_dprime_box_sequence(),
                                                                 self.dprime_box[0],
                                                                 self.get_cbox_sequence(),
                                                                 self.c_box[0],
                                                                 self.get_dbox_sequence(),
                                                                 self.d_box[0],
                                                                 self.sequence)

            return plexy
        else:
            raise Exception("Wrong composition of boxes for %s" % self.snor_id)


    def get_d_interaction_region(self, length):
        region = self.sequence[self.d_box[0] - length - 1: self.d_box[0] - 1]
        if len(region) == length:
            return region
        else:
            raise InteractionRegionTooShortException("Cannot extract D-box interaction region" +
                                                     " for %s (too short)" % self.snor_id)

    def get_dprime_interaction_region(self, length):
        if self.dprime_box:
            region = self.sequence[self.dprime_box[0] - length - 1: self.dprime_box[0] - 1]
            if len(region) == length:
                return region
            else:
                raise InteractionRegionTooShortException("Cannot extract D'-box interaction region" +
                                                         " for %s (too short)" % self.snor_id)
        else:
            return None


class HACA_snoRNA(snoRNA):

    """A class representing H/ACA box snoRNA"""

    def __init__(self, snor_id, organism, chrom, start, end, strand,
                 sequence, snor_type, h_box, aca_box, **kwargs):
        """@todo: to be defined1.

        Args:
            snor_id (str): a unique id for snoRNA
            organism (str): species in which snoRNA can be found
            chrom (str): @todo
            start (int): @todo
            end (int): @todo
            strand (str): @todo
            sequence (str): @todo
            snor_type (str): @todo
            h_box (str): H box position: "1..5".
                             Start and end should be separated by "..".
            aca_box (str): ACA box position: "1..5".
                             Start and end should be separated by "..".

        Kwargs:
            alias (str): anlternative name for snoRNA
            gene_name (str): snoRNA gene name
            accession (str): accession for the snoRNA (eg. NCBI)
            modified_sites (dict): dictionary of modified sites: {rna: {position, nucleotide}}
            host_gene (str): host gene
            host_id (str): id for host locus
            organization (str): organization of the locus
            note (str): additional information about snoRNA


        """
        snoRNA.__init__(self, snor_id, organism, chrom, start, end, strand, sequence, snor_type, **kwargs)
        self.__asign_h_box(h_box)
        self.__asign_aca_box(aca_box)
        self.__validate_haca()

    def __asign_h_box(self, h_box):
        """
        Assign H-box to snoRNA based on the string provided. H-box should be in the
        form of "H-box_start..H-box_end"
        """
        if h_box != None:
            hstart, hend = [int(i) for i in h_box.split("..")]
            self.h_box = (hstart, hend)
        else:
            warn("No H-box was assigned for %s" % self.snor_id)
            self.h_box = None

    def __asign_aca_box(self, aca_box):
        """
        Assign H-box to snoRNA based on the string provided. H-box should be in the
        form of "H-box_start..H-box_end"
        """
        if aca_box != None:
            acastart, acaend = [int(i) for i in aca_box.split("..")]
            self.aca_box = (acastart, acaend)
        else:
            warn("No ACA-box was assigned for %s" % self.snor_id)
            self.aca_box = None

    def __validate_haca(self):
        """Validate if H/ACA box snoRNA is correct and can be
        used in calculations"""
        if self.h_box:
            if (self.h_box[1] > len(self.sequence)) or (self.h_box[0] < 0):
                raise WrongBoxException("Box outside snoRNA")

    def get_hbox_sequence(self):
        if self.h_box:
            return self.sequence[self.d_box[0] - 1: self.d_box[1]]
        else:
            raise NoBoxException("No H-box in this H/ACA snoRNA")

    def get_acabox_sequence(self):
        if self.aca_box:
            return self.sequence[self.aca_box[0] - 1: self.aca_box[1]]
        else:
            raise NoBoxException("No ACA-box in this H/ACA snoRNA")

    def get_left_stem_sequence(self):
        """Get the sequence of the left stem of the snoRNA"""
        return self.sequence[:self.h_box[0]]

    def get_right_stem_sequence(self):
        """Get the sequence of the right stem of the snoRNA"""
        return self.sequence[self.h_box[1]:]

    def get_left_stem_fasta(self):
        """Get the fasta string of the left stem of the snoRNA"""
        return ">%s_stem1\n%s" % (self.snor_id, self.get_left_stem_sequence())

    def get_right_stem_fasta(self):
        """Get the fasta string of the right stem of the snoRNA"""
        return ">%s_stem2\n%s" % (self.snor_id, self.get_right_stem_sequence())
