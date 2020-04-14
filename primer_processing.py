import os
import subprocess

class PrimerProcess:
    def __init__(self, blat_path, primer_path, ref_path, prefix):
        self.blat_path = blat_path
        self.in_primer_file = primer_path
        self.ref_path = ref_path
        self.prefix = prefix
        self.pools = set()
        self.primer_info = dict()
        self.generate_fasta_for_primers()
        self.align_primer()
        self.separate_m8_file()
        self.screen_m8_file()

    def generate_fasta_for_primers(self):
        inf = open(self.in_primer_file, "r")  # name pool seq length
        out_fa = open("primer.fa", "w")
        header = inf.readline()
        for line in inf:
            buffer = line.strip().split()
            name, pool, seq, length = buffer[0], buffer[1], buffer[2], int(buffer[3])
            self.primer_info[name] = (pool, length)
            self.pools.add(pool)
            out_fa.write(">%s\n" % name)
            out_fa.write("%s\n" % seq)
        out_fa.close()
        return

    def align_primer(self):
        cmd = self.blat_path + " -t=dna -q=dna -out=blast8 -tileSize=6  -stepSize=3 -minMatch=1 -repMatch=1000000 -minScore=0 -minIdentity=80 " \
                             + self.ref_path\
                             + " ./primer.fa ./primer2virus.m8"
        p = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
        return

    def separate_m8_file(self):
        outs = dict()
        for pcr_pool in self.pools:
            name = self.prefix + "." + pcr_pool + ".m8"
            outs[pcr_pool] = open(name, "w")
        inf = open("primer2virus.m8", "r")
        for line in inf:
            buffer = line.strip().split()
            p_name = buffer[0]
            pool, l = self.primer_info[p_name]
            out = outs[pool]
            out.write(line)
        for pcr_pool in self.pools:
            outs[pcr_pool].close()
        inf.close()

        return

    def screen_m8_file(self):
        for pcr_pool in self.pools:
            in_name = self.prefix + "." + pcr_pool + ".m8"
            out_name = self.prefix+ "_highconf." + pcr_pool + ".tsv"
            infile = open(in_name, "r")
            outfile = open(out_name, "w")
            outfile.write("Primer_name\tTarget\tAlignment_len\tQuery_start\tQuery_end\tTarget_start\tTarget_end\n")
            for line in infile:
                buffer = line.split()
                name = buffer[0]
                target = buffer[1]
                percentage = float(buffer[2])
                aligned_len = int(buffer[3])
                qs = int(buffer[6])
                qe = int(buffer[7])
                ts = int(buffer[8])
                te = int(buffer[9])
                if ts > te:
                    tmp = ts
                    ts = te
                    te = tmp
                pool, query_len = self.primer_info[name]

                if (aligned_len + 5) > query_len and percentage > 70.0:
                    outfile.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % (name, target, aligned_len, qs, qe, ts, te))
            infile.close()
            outfile.close()
        return
