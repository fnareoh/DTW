import pysam
from collections import defaultdict


class Read_Badread:
    """ stores all information needed about a read """

    def __init__(self, fastq_text, sequence, quality, sam_id, sam_alignement, G):
        self.sequence = "".join(sequence.split())
        self.quality = quality
        info = fastq_text[1:].split(" ")
        self.id = info[0]
        self.strand = info[1].split(",")[1]
        # position from which the read was generated as given by the fastq header
        self.position = int(info[1].split(",")[2].split("-")[0])
        # the ending position from which the read was generated as given by the fastq header
        self.end_position = int(info[1].split(",")[2].split("-")[1])
        if self.strand == "-strand":
            # Correction due to badread header on - strand
            self.position, self.end_position = max(0, len(G) - self.end_position), max(
                0, len(G) - self.position
            )
        self.length = int(info[2].split("=")[-1])
        self.error_free_length = int(info[3].split("=")[-1])
        self.read_identity = float(info[4].split("=")[-1].split("%")[0])
        self.sam_alignement = sam_alignement
        assert self.id == sam_id
        assert self.length == len(self.sequence)

    def __str__(self):
        res = f"Read Id = {self.id}"
        return res

    def sequence_str(self):
        return f"\nSequence = {self.sequence}"

    def detail_str(self):
        res = f"Read Id = {self.id}"
        res += f"\nStrand = {self.strand}, Start position={self.position}, End position={self.end_position}"
        res += f"\nLength = {self.length}, Error-free_length={self.error_free_length}, Sam alignment={self.sam_alignement}"
        res += f", Read identity = {self.read_identity}"
        return res


class Read_Nanosim:
    def __init__(
        self, fastq_text, sequence, quality, sam_id, sam_alignement, error_length
    ):
        self.sequence = "".join(sequence.split())
        self.quality = quality
        info = fastq_text[1:].split("_")
        self.id = fastq_text[1:-1]
        self.position = int(info[1])
        self.strand = info[4]
        self.head = int(info[-3])
        self.middle = int(info[-2])
        self.tail = int(info[-1])
        self.length = len(self.sequence)
        self.read_identity = float(100 * (self.length - error_length) / self.length)
        self.sam_alignement = sam_alignement
        assert self.id == sam_id
        assert self.length == len(self.sequence)

    def detail_str(self):
        res = f"Read Id = {self.id}"
        res += f"\nStrand = {self.strand}, Start position={self.position}"
        res += f"\nLength = {self.length}, Sam alignment={self.sam_alignement}"
        res += f", head, middle, tail = {self.head,self.middle,self.tail}"
        return res


def parse_input(genome_file, read_pref, simulator="badread"):
    """
    Takes as input a fasta file for the genome, a fastq file and a sam file
    for the reads and output the genome (a string) G and a list of Read objects.
    """
    file_G = open(genome_file, "r")
    G = "".join(["".join(l.split()) for l in file_G.readlines()[1:]])
    file_G.close()
    fastq_R = open(read_pref + ".fastq", "r")
    sam_R = pysam.AlignmentFile(read_pref + ".sam", "r")
    list_R_al = []
    list_R_not_al = []
    nb_read_not_al = 0
    nb_read = 0
    if simulator == "nanosim":
        error_file = open(f"{read_pref}.error_length")
        error_file.readline()
        err_length_dict = defaultdict(int)
        for line in error_file:
            seq_name, error_length = line.split(",")
            err_length_dict[seq_name] = int(error_length)
        error_file.close()
    for read in sam_R.fetch():
        nb_read += 1
        text = fastq_R.readline()
        assert text[0] == "@"
        seq = fastq_R.readline()
        optional = fastq_R.readline()
        assert optional[0] == "+"
        qualities = fastq_R.readline()
        if simulator == "badread":
            R = Read_Badread(
                text, seq, qualities, read.query_name, read.reference_start, G
            )
        elif simulator == "nanosim":
            id = text[1:-1]
            R = Read_Nanosim(
                text,
                seq,
                qualities,
                read.query_name,
                read.reference_start,
                int(err_length_dict[id]),
            )
        else:
            raise Error("Unkwoned simulator!")
        if R.sam_alignement == -1:
            nb_read_not_al += 1
            list_R_not_al.append(R)
        else:
            list_R_al.append(R)
    print(f"Genome size: {len(G)}")
    print(f"Number of read not alligned in the sam file = {nb_read_not_al}")
    print(f"Total number of read = {nb_read}")
    return G, list_R_al, list_R_not_al
