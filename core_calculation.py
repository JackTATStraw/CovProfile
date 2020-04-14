from BCBio import GFF
import pysam
from scipy.optimize import minimize
import numpy as np
import math

tsv = "tsv"
out = "out"
bam = "bam"

class Calculation:
    def __init__(self, gff_path, prefix, pcr_pool, offset, ref_name):
        self.gff_path = gff_path
        self.prefix = prefix
        self.pcr_pool = pcr_pool
        self.offset = offset
        self.ref_name = ref_name
        self.gene_locations = dict()
        self.pos_primer = list()
        self.neg_primer = list()
        self.options = self.infer_options()
        self.read_gff()
        self.read_primer()
        self.do_calculation()
        return

    def out_primer_efficiencies(self, _primer_efficiencies: list):
        n_pool = len(self.pcr_pool)
        for idx in range(0, n_pool):
            primers_pairs = self.pos_primer[idx]
            primer_efficiencies = _primer_efficiencies[idx]
            out_name = self.options[out][idx] + ".primer_eff"
            sample = self.options[out][idx]
            outf = open(out_name, "w")
            outf.write("sample")
            for i in range(0, len(primers_pairs)):
                primer_name = primers_pairs[i][0]
                outf.write("\t%s" % primer_name)
            outf.write("\n%s" % sample)
            for i in range(0, len(primers_pairs)):
                primer_efficiency = primer_efficiencies[i]
                outf.write("\t%f" % primer_efficiency)
            outf.close()

    def out_gene_exp(self, gene_exp: dict):
        n_pool = len(self.pcr_pool)
        for idx in range(0, n_pool):
            out_name = self.options[out][idx] + ".gene_exp"
            sample = self.options[out][idx]
            outf = open(out_name, "w")
            outf.write("sample")
            for gene_id in gene_exp.keys():
                outf.write("\t%s" % gene_id)
            outf.write("\n%s" % sample)
            for gene_id in gene_exp.keys():
                exp = gene_exp[gene_id]
                if exp == "NA":
                    outf.write("\tNA")
                else:
                    outf.write("\t%f" % exp)
            outf.close()

    def do_calculation(self):
        gene_primer_map, primer_link_two_gene, primer_cover_entire_gene, primer_pairs_on_gene = self.map_gene_primer()
        primer_pairs_on_gene_reverse = self.reverse_gene_primer(primer_pairs_on_gene)
        primer_pairs_reads_count = self.map_reads_primer()
        gene_primer_reads_count = self.map_reads_gene_primer(gene_primer_map)
        primer_off_target = self.get_primer_off_target(primer_pairs_reads_count)
        primer_on_copy = self.get_primer_on_copy(primer_pairs_reads_count, primer_cover_entire_gene, primer_pairs_on_gene,
                                            gene_primer_reads_count)

        primer_efficiencies, gene_exp, K, Copy, a = self.optimize_all(primer_pairs_reads_count, gene_primer_reads_count,
                                                                 primer_on_copy, primer_pairs_on_gene_reverse,
                                                                 primer_off_target)
        self.out_primer_efficiencies(primer_efficiencies)
        self.out_gene_exp(gene_exp)

        n_pool = len(self.pcr_pool)
        for idx in range(0, n_pool):
            o_name = self.options[out][idx] + ".stats"
            f = open(o_name, "w")
            f.write("K\tCopy\ta\n")
            f.write("%f\t%f\t%f\n" % (K[idx], Copy, a))
            f.close()
        return

    def infer_options(self):
        options = dict()
        options[tsv] = list()
        options[out] = list()
        options[bam] = list()

        for pool in self.pcr_pool:
            options[tsv].append("../primer_align/" + self.prefix + "_highconf." + pool + ".tsv")
            options[out].append(self.prefix + "." + pool)
            options[bam].append("../virus_align/" +self.prefix + "." + pool + ".sorted.bam")
        return options

    def read_gff(self):
        inf = open(self.gff_path, "r")
        e = GFF.GFFExaminer()
        # tmp = e.available_limits(inf)
        # pprint.pprint(tmp)

        for r in GFF.parse(inf):
            for record in r.features:
                if len(record.sub_features) >= 1:
                    self.gene_locations[record.id] = (record.location.nofuzzy_start, record.location.nofuzzy_end)
        inf.close()

    def read_primer(self):
        pos_primer = list()
        neg_primer = list()
        for inn in self.options[tsv]:
            inf = open(inn, "r")
            lines = inf.readlines()
            pos_primer.append([])
            neg_primer.append([])
            for i in range(1, len(lines), 2):
                line = lines[i]
                buffer = line.strip().split()
                p_name = buffer[0]
                p_start = int(buffer[5]) - 1
                p_end = int(buffer[6]) - 1
                line = lines[i + 1]
                buffer = line.strip().split()
                n_name = buffer[0]
                n_start = int(buffer[5]) - 1
                n_end = int(buffer[6]) - 1
                pos_primer[-1].append((p_name, p_start, p_end))
                neg_primer[-1].append((n_name, n_start, n_end))
            inf.close()
        self.pos_primer = list(pos_primer)
        self.neg_primer = list(neg_primer)

    def map_gene_primer(self):
        n_pool = len(self.pcr_pool)
        _gene_primer_map = list()
        _primer_link_two_gene = list()
        _primer_cover_entire_gene = list()
        _primer_on_gene = list()
        _already_set_primers = list()

        for idx in range(0, n_pool):
            pos_primers = self.pos_primer[idx]
            neg_primers = self.neg_primer[idx]
            gene_primer_map = dict()
            primer_link_two_gene = dict()
            primer_cover_entire_gene = dict()
            primer_on_gene = dict()
            already_set_primers = dict()
            for gene_id in self.gene_locations.keys():
                gene_start, gene_end = self.gene_locations[gene_id]
                if gene_id not in gene_primer_map.keys():
                    gene_primer_map[gene_id] = []
                    primer_on_gene[gene_id] = []
                for i in range(0, len(pos_primers)):
                    p_name, p_start, p_end = pos_primers[i]
                    n_name, n_start, n_end = neg_primers[i]
                    if p_start >= gene_start - 5 and n_end <= gene_end + 5:
                        primer_on_gene[gene_id].append(i)
                already_set_primers_gene = set()
                for i in range(0, len(pos_primers)):
                    p_name, p_start, p_end = pos_primers[i]
                    n_name, n_start, n_end = neg_primers[i]

                    if p_end + 30 < gene_start and gene_start + 30 < n_start:
                        # start of the gene to negative primer on gene
                        gene_primer_map[gene_id].append((i, 's', 'n'))
                        already_set_primers_gene.add(i)
                        if i in already_set_primers.keys():
                            if i in primer_link_two_gene.keys():
                                primer_link_two_gene[i].add(already_set_primers[i])
                            else:
                                primer_link_two_gene[i] = set()
                                primer_link_two_gene[i].add(already_set_primers[i])
                                primer_link_two_gene[i].add(gene_id)

                        already_set_primers[i] = gene_id
                    if p_end + 30 < gene_end and gene_end + 30 < n_start:
                        # positive primer to end of gene
                        gene_primer_map[gene_id].append((i, 'e', 'p'))
                        if i not in already_set_primers_gene:
                            already_set_primers[i] = gene_id
                        else:
                            primer_cover_entire_gene[i] = gene_id
            _gene_primer_map.append(gene_primer_map)
            _primer_link_two_gene.append(primer_link_two_gene)
            _primer_cover_entire_gene.append(primer_cover_entire_gene)
            _primer_on_gene.append(primer_on_gene)
            _already_set_primers.append(already_set_primers)
        return _gene_primer_map, _primer_link_two_gene, _primer_cover_entire_gene, _primer_on_gene

    def reverse_gene_primer(self, _primer_pairs_on_gene: list):
        _primer_pairs_on_gene_reverse = list()
        n_pool = len(_primer_pairs_on_gene)
        for idx in range(0, n_pool):
            primer_pairs_on_gene = _primer_pairs_on_gene[idx]
            primer_pairs_on_gene_reverse = dict()
            for gene in primer_pairs_on_gene.keys():
                primers = primer_pairs_on_gene[gene]
                for primer in primers:
                    primer_pairs_on_gene_reverse[primer] = gene
            _primer_pairs_on_gene_reverse.append(primer_pairs_on_gene_reverse)
        return _primer_pairs_on_gene_reverse

    def map_reads_primer(self):
        n_pool = len(self.pcr_pool)
        _primer_pair_reads_count = list()
        for idx in range(0, n_pool):
            pos_primers = self.pos_primer[idx]
            neg_primers = self.neg_primer[idx]
            bam_name = self.options[bam][idx]
            inf = pysam.AlignmentFile(bam_name, "rb")
            primer_pair_reads_count = dict()
            for i in range(0, len(pos_primers)):
                p_name, p_start, p_end = pos_primers[i]
                n_name, n_start, n_end = neg_primers[i]
                reads_count = 0
                record: pysam.AlignedSegment
                for record in inf.fetch(self.ref_name, p_start - self.offset, n_end + self.offset):
                    # TODO: consider secondary alignment
                    if record.is_unmapped or record.is_secondary:
                        continue
                    # TODO: don't count a single alignment record for multiple times
                    if record.reference_start > p_start - self.offset and record.reference_start < p_start + self.offset:
                        if record.reference_end > n_end - self.offset and record.reference_end < n_end + self.offset:
                            reads_count += 1
                primer_pair_reads_count[i] = reads_count
            inf.close()
            _primer_pair_reads_count.append(primer_pair_reads_count)
        return _primer_pair_reads_count  # = primer_efficiency * virus_genome_copy * 2^k = M

    def map_reads_gene_primer(self, _gene_primer_map: list):
        n_pool = len(self.pcr_pool)
        _gene_primer_reads_count = list()
        for idx in range(0, n_pool):
            pos_primers = self.pos_primer[idx]
            neg_primers = self.neg_primer[idx]
            gene_primer_map = _gene_primer_map[idx]
            bam_name = self.options[bam][idx]
            inf = pysam.AlignmentFile(bam_name, "rb")
            gene_primer_reads_count = dict()
            for gene_id in gene_primer_map:
                gene_start, gene_end = self.gene_locations[gene_id]
                gene_primer_reads_count[gene_id] = dict()
                for primer in gene_primer_map[gene_id]:
                    primer_id, gene_loc, primer_direction = primer
                    start: int
                    end: int
                    if gene_loc == 's':
                        # equivalent to negative primer
                        start = gene_start - self.offset
                        end = neg_primers[primer_id][2] + self.offset
                    else:
                        start = pos_primers[primer_id][1] - self.offset
                        end = gene_end + self.offset
                    reads_count = 0
                    for record in inf.fetch(self.ref_name, start, end):
                        if record.reference_start > start - self.offset and record.reference_start < start + self.offset:
                            if record.reference_end > end - self.offset and record.reference_end < end + self.offset:
                                reads_count += 1
                    if primer_id not in gene_primer_reads_count[gene_id].keys():
                        gene_primer_reads_count[gene_id][primer_id] = reads_count
                    else:
                        gene_primer_reads_count[gene_id][primer_id] += reads_count
            inf.close()
            _gene_primer_reads_count.append(gene_primer_reads_count)
        return _gene_primer_reads_count  # = primer_efficiency/2 * gene_expression * K = N

    def get_primer_off_target(self, _primer_pair_reads_count: list):
        _primer_off_target = list()
        for idx in range(0, len(_primer_pair_reads_count)):
            primer_off_target = set()
            primer_pair_reads_count = _primer_pair_reads_count[idx]
            for i in primer_pair_reads_count.keys():
                if primer_pair_reads_count[i] < 15:
                    primer_off_target.add(i)
            _primer_off_target.append(primer_off_target)

        return _primer_off_target

    def get_primer_on_copy(self, _primer_pair_reads_count: list, _primer_cover_entire_gene: list, _primer_pairs_on_gene: list,
                           _gene_primer_reads_count: list):
        n_pool = len(self.pcr_pool)
        _primer_on_copy = list()
        for idx in range(0, n_pool):
            primer_pair_reads_count = _primer_pair_reads_count[idx]
            primer_cover_entire_gene = _primer_cover_entire_gene[idx]
            primer_pairs_on_gene = _primer_pairs_on_gene[idx]
            gene_primer_reads_count = _gene_primer_reads_count[idx]
            primer_on_copy = set()
            primer_on_gene = set()
            for i in primer_cover_entire_gene.keys():
                primer_on_copy.add(i)
            for gene_id in primer_pairs_on_gene:
                primers1 = primer_pairs_on_gene[gene_id]
                primers2 = gene_primer_reads_count[gene_id]
                for i in primers1:
                    primer_on_gene.add(i)
                for i in primers2.keys():
                    primer_on_gene.add(i)

            for i in primer_pair_reads_count.keys():
                if i not in primer_on_gene:
                    primer_on_copy.add(i)
            _primer_on_copy.append(primer_on_copy)
        return _primer_on_copy

    def optimize_all(self, _primer_pairs_reads_count: list, _gene_primer_reads_count: list, _primer_on_copy: list,
                     _primer_on_gene: list, _primer_off_target: list):
        n_pool = len(_primer_pairs_reads_count)
        Ks = list()
        for idx in range(0, n_pool):
            Ks.append(20)
        K = 20
        Copy = 2
        a = 2
        _primer_eff = list()
        gene_exp = dict()
        x0 = list()
        for idx in range(0, n_pool):
            primer_eff = dict()
            for i in _primer_pairs_reads_count[idx].keys():
                primer_eff[i] = 1.5
            _primer_eff.append(primer_eff)
            x0 += list(primer_eff.values())
        n_primer = len(x0)
        for i in _gene_primer_reads_count[0].keys():
            gene_exp[i] = 2
        primer_accross_gene = self.get_primer_accross_gene(_gene_primer_reads_count)
        x0 += list(gene_exp.values())
        x0.append(a)
        x0.append(Copy)
        x0 += Ks
        # x0.append(K)
        x0 = np.asarray(x0)
        arg = list()
        arg.append(_primer_pairs_reads_count)
        arg.append(_gene_primer_reads_count)
        arg.append(primer_accross_gene)
        arg.append(_primer_on_copy)
        arg.append(_primer_on_gene)
        arg.append(_primer_off_target)

        arg2 = (n_primer,)
        arg3 = (len(gene_exp.keys()),)
        bond_exp = self.con_exp(arg3)
        bond_eff = self.con_eff(arg2)
        bond = bond_eff + bond_exp
        bond.append((2, None))
        bond.append((1, None))  # Copy
        for idx in range(0, n_pool):
            bond.append((20, 30))  # K
        r = minimize(self.optimize_sub, x0=x0, args=arg, method="L-BFGS-B", bounds=bond)
        Ks = list()
        for idx in range(0, n_pool):
            Ks.append(r.x[-1 - idx])
        Ks.reverse()
        Copy = r.x[-1 - n_pool]
        a = r.x[-2 - n_pool]
        n_primer = 0
        _primer_eff = list()
        for idx in range(0, n_pool):
            primer_eff = dict()

            for i in _primer_pairs_reads_count[idx].keys():
                primer_eff[i] = r.x[i + n_primer]
            n_primer += len(_primer_pairs_reads_count[idx].keys())
            _primer_eff.append(primer_eff)

        n_gene = len(_gene_primer_reads_count[0].keys())

        for i in _gene_primer_reads_count[0].keys():
            gene_exp[i] = r.x[n_primer]
            if gene_exp[i] == 2.0:
                gene_exp[i] = "NA"
            n_primer += 1

        count = 0
        # for i in primer_pairs_reads_count.keys():
        #    primer_eff[i] = r.x[count]
        #    count+= 1
        # for i in gene_primer_reads_count.keys():
        #    gene_exp[i] = r.x[count]
        #    count += 1

        return _primer_eff, gene_exp, Ks, Copy, a,

    def get_primer_accross_gene(self, _gene_primer_reads_count: list):
        n_pool = len(_gene_primer_reads_count)
        _primer_accross_gene = list()
        for idx in range(0, n_pool):
            gene_primer_reads_count = _gene_primer_reads_count[idx]
            primer_accross_gene = dict()
            for gene in gene_primer_reads_count.keys():
                primers = gene_primer_reads_count[gene]
                for primer in primers.keys():
                    if primer not in primer_accross_gene.keys():
                        primer_accross_gene[primer] = list()
                    if primers[primer] > 3:
                        primer_accross_gene[primer].append(gene)
            tmp = dict(primer_accross_gene)
            for primer in primer_accross_gene.keys():
                if len(primer_accross_gene[primer]) == 0:
                    tmp.pop(primer)
            _primer_accross_gene.append(tmp)
        return _primer_accross_gene

    def con_eff(self, args):
        p_count = args[0]
        cons = list()
        for i in range(0, p_count):
            c = (1, 2)
            cons.append(c)
        return cons

    def con_exp(self, args):
        p_count = args[0]
        cons = list()
        for i in range(0, p_count):
            c = (1, None)
            cons.append(c)
        return cons

    @staticmethod
    def optimize_sub(x0, args):
        _primer_pair_reads_count, _gene_primer_reads_count, _primer_accross_gene, _primer_on_copy, _primer_on_gene, _primer_off_target = \
        args[0], args[1], args[2], args[3], args[4], args[5]
        n_pool = len(_primer_pair_reads_count)
        _p_eff = list()
        gene_exp = list()
        gene_name2id = dict()
        count = 0
        for gene in _gene_primer_reads_count[0].keys():
            gene_name2id[gene] = count
            count += 1
        count = 0
        Ks = list()
        for idx in range(0, n_pool):
            Ks.append(x0[-1 - idx])
        Ks.reverse()
        Copy = x0[-1 - n_pool]
        a = x0[-2 - n_pool]
        n_primer = 0
        for idx in range(0, n_pool):
            _p_eff.append([])
            p_eff = _p_eff[-1]
            for i in _primer_pair_reads_count[idx].keys():
                p_eff.append(x0[i + n_primer])
            n_primer += len(_primer_pair_reads_count[idx].keys())

        n_gene = len(_gene_primer_reads_count[0].keys())

        for i in range(n_primer, n_gene + n_primer):
            gene_exp.append(x0[i])
        f = 0
        for idx in range(0, n_pool):
            primer_pair_reads_count = _primer_pair_reads_count[idx]
            primer_accross_gene = _primer_accross_gene[idx]
            primer_off_target = _primer_off_target[idx]
            primer_on_gene = _primer_on_gene[idx]
            primer_on_copy = _primer_on_copy[idx]
            gene_primer_reads_count = _gene_primer_reads_count[idx]
            p_eff = _p_eff[idx]
            K = Ks[idx]
            special_primer = set()
            for primer_id in primer_pair_reads_count.keys():
                if primer_id in primer_accross_gene.keys():
                    gene_names = primer_accross_gene[primer_id]
                    if primer_id not in primer_off_target:
                        f += (K * math.log2(p_eff[primer_id]) + math.log2(Copy) - math.log2(
                            primer_pair_reads_count[primer_id]) - math.log2(a)) ** 2
                    for gene_name in gene_names:
                        special_primer.add(primer_id)
                        gene_id = gene_name2id[gene_name]
                        if gene_primer_reads_count[gene_name][primer_id] < 10:
                            continue
                        f += (math.log2(K) + math.log2(p_eff[primer_id]) + math.log2(gene_exp[gene_id]) - 1 - math.log2(
                            gene_primer_reads_count[gene_name][primer_id]) - math.log2(a)) ** 2
                elif primer_id in primer_on_gene.keys():
                    gene_id = gene_name2id[primer_on_gene[primer_id]]
                    if primer_id not in primer_off_target:
                        f += (K * math.log2(p_eff[primer_id]) + math.log2(Copy + gene_exp[gene_id]) - math.log2(
                            primer_pair_reads_count[primer_id]) - math.log2(a)) ** 2
                elif primer_id in primer_on_copy:
                    if primer_id not in primer_off_target:
                        f += (K * math.log2(p_eff[primer_id]) + math.log2(Copy) - math.log2(
                            primer_pair_reads_count[primer_id]) - math.log2(a)) ** 2
                else:
                    if primer_id not in primer_off_target:
                        f += (K * math.log2(p_eff[primer_id]) + math.log2(Copy) - math.log2(
                            primer_pair_reads_count[primer_id]) - math.log2(a)) ** 2

            for j in range(0, len(p_eff)):
                eff = p_eff[j]
                if j in special_primer:
                    f += math.log2(eff) ** 2
                else:
                    f += math.log2(eff) ** 2
            # for exp in gene_exp:
            #    f += math.log2(exp) ** 2

            f += 0.5 * K
        f += math.log2(Copy) ** 2
        f += math.log2(a) ** 2
        return f
