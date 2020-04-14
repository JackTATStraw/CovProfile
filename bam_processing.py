import subprocess
import pysam

class BamProcess:
    def __init__(self, Minimap_path, samtools_path, fastq_path, ref_path, prefix, pcr_pools, primer_info, offset):
        self.minimap_path = Minimap_path
        self.fastq_path = fastq_path
        self.ref_path = ref_path
        self.prefix = prefix
        self.samtools_path = samtools_path
        self.pcr_pools = pcr_pools
        self.primer_info = primer_info
        self.offset = offset
        self.create_ref_index()
        self.align_fastq()
        self.filter_bam()
        self.sort_bam()
        pos_primer, neg_primer = self.read_primer_alignment()
        self.split_bam(pos_primer, neg_primer)
        return

    def create_ref_index(self):
        cmd = self.minimap_path + " -x map-ont -d ref.mmi " + self.ref_path
        p = subprocess.call(cmd, shell=True, stderr=subprocess.PIPE)

    def align_fastq(self):
        cmd = self.minimap_path + " -ax map-ont -o " + self.prefix+".sam" + " ref.mmi " + self.fastq_path
        p = subprocess.call(cmd, shell=True, stderr=subprocess.PIPE)

    def filter_bam(self):
        cmd = self.samtools_path + " view -F 0x104 -b -h -o " + self.prefix+".aln.bam " + self.prefix+".sam"
        p = subprocess.call(cmd, shell=True, stderr=subprocess.PIPE)

    def sort_bam(self):
        cmd = self.samtools_path + " sort -o " + self.prefix+".aln.sorted.bam " + self.prefix+".aln.bam"
        p = subprocess.call(cmd, shell=True, stderr=subprocess.PIPE)

    def index_bam(self, bam_path):
        cmd = self.samtools_path + " index " + bam_path
        p = subprocess.call(cmd, shell=True, stderr=subprocess.PIPE)

    def read_primer_alignment(self):
        pos_primer = dict()
        neg_primer = dict()
        for pool in self.pcr_pools:
            in_name = "../primer_align/" + self.prefix + "_highconf." + pool + ".tsv"
            inf = open(in_name, "r")
            lines = inf.readlines()
            for i in range(1, len(lines), 2):
                line = lines[i]
                buffer = line.strip().split()
                p_name = buffer[0]
                pool, l = self.primer_info[p_name]

                p_start = int(buffer[5]) - 1
                p_end = int(buffer[6]) - 1
                line = lines[i + 1]
                buffer = line.strip().split()
                n_name = buffer[0]
                n_start = int(buffer[5]) - 1
                n_end = int(buffer[6]) - 1
                pos_primer[p_name] = (p_start, p_end)
                neg_primer[n_name] = (n_start, n_end)
            inf.close()
        return pos_primer, neg_primer

    def split_bam(self, pos_primers, neg_primers):
        outs = dict()
        inf = pysam.AlignmentFile(self.prefix + ".aln.sorted.bam", "rb")
        tmp_name = self.prefix + ".extra.bam"
        out_not_on_pcr = pysam.AlignmentFile(tmp_name, "wb", template=inf)
        for pool in self.pcr_pools:
            out_name = self.prefix + "." + pool + ".sorted.bam"
            outs[pool] = pysam.AlignmentFile(out_name, "wb", template=inf)
        record: pysam.AlignedSegment
        for record in inf:
            record_written = False
            for pos_primer in pos_primers.keys():
                p_start, p_end = pos_primers[pos_primer]
                pool, l = self.primer_info[pos_primer]
                if record.reference_start > p_start - self.offset and record.reference_start < p_start + self.offset:
                    outs[pool].write(record)
                    record_written = True
                    break
            if not record_written:
                for neg_primer in neg_primers.keys():
                    n_start, n_end = neg_primers[neg_primer]
                    pool, l = self.primer_info[neg_primer]
                    if record.reference_end > n_end - self.offset and record.reference_end < n_end + self.offset:
                        outs[pool].write(record)
                        record_written = True
                        break
            if not record_written:
                out_not_on_pcr.write(record)

        for pool in self.pcr_pools:
            outs[pool].close()
            out_name = self.prefix + "." + pool + ".sorted.bam"
            self.index_bam(out_name)
        inf.close()
        out_not_on_pcr.close()
        return