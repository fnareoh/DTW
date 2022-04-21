import pysam


class Read:
    """ stores all information needed about a read """

    def __init__(
        self, fastq_text, sequence, quality, sam_id, sam_alignement, simulator="badread"
    ):
        self.sequence = sequence
        self.quality = quality
        if simulator == "nanosim":
            info = fastq_text[1:].split("_")
            self.id = fastq_text[1:-1]
            self.position = int(info[1])
            self.end_position = -1
            self.length = len(self.sequence) - 1
        else:
            info = fastq_text[1:].split(" ")
            self.id = info[0]
            self.strand = info[1].split(",")[1]
            # position from which the read was generated as given by the fastq header
            self.position = int(info[1].split(",")[2].split("-")[0])
            # the ending position from which the read was generated as given by the fastq header
            self.end_position = int(info[1].split(",")[2].split("-")[1])
            self.length = int(info[2].split("=")[-1])
            self.error_free_length = int(info[3].split("=")[-1])
            self.read_identity = info[4].split("=")[-1]
        self.sam_alignement = sam_alignement
        assert self.id == sam_id
        assert self.length == len(self.sequence) - 1

    def __str__(self):
        res = f"Read Id = {self.id}"
        return res

    def detail_str(self):
        res = f"Read Id = {self.id}"
        res += f"\nStrand = {self.strand}, Position={self.position}"
        res += f"\nLength = {self.length}, Error-free_length={self.error_free_length}"
        res += f", Read identity = {self.read_identity}"
        res += f"\nSequence = {self.sequence}"
        return res


def parse_input(genome_file, read_fastq_file, read_sam_file, simulator="badread"):
    """
    Takes as input a fasta file for the genome, a fastq file and a sam file
    for the reads and output the genome (a string) G and a list of Read objects.
    """
    file_G = open(genome_file, "r")
    G = "".join(file_G.readlines()[1:])
    file_G.close()
    fastq_R = open(read_fastq_file, "r")
    sam_R = pysam.AlignmentFile(read_sam_file, "r")
    list_R_al = []
    list_R_not_al = []
    nb_read_not_al = 0
    nb_read = 0
    for read in sam_R.fetch():
        nb_read += 1
        text = fastq_R.readline()
        assert text[0] == "@"
        seq = fastq_R.readline()
        optional = fastq_R.readline()
        assert optional[0] == "+"
        qualities = fastq_R.readline()
        R = Read(text, seq, qualities, read.query_name, read.reference_start, simulator)
        if R.sam_alignement == -1:
            nb_read_not_al += 1
            list_R_not_al.append(R)
        else:
            list_R_al.append(R)
    print(f"Genome size: {len(G)}")
    print(f"Read size: {len(list_R_al[0].sequence)}")
    print(f"Number of read not alligned in the sam file = {nb_read_not_al}")
    print(f"Total number of read = {nb_read}")
    return G, list_R_al, list_R_not_al
