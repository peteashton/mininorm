""" Module to implement quit and simple handling of FASTQ sequences in FASTQ files (and possibly gzipped FASTQ files).  
Implements the classes FastqReader, FastqSeq and FastqFormatError (the last used just for raising exceptions for 
invalid FASTQ file).  Format and error checking is minimal, and quality strings are not decoded, to maintain the 
speed of the module.
"""

import gzip

complementer = str.maketrans("ACGTYRacgtyr", "TGCARYtgcary")

class FastqFormatError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class FastqReader:
    """ Class that allows iteration over FASTQ files.  Can read compressed (gzipped) files, and handles decoding 
    from utf-8 to native python strings.  Performs rudimentary format checking, but is designed to be much faster 
    than the Bio.SeqIO module.  A FastqReader object is a fully compliant interator, and can be used multiple times
    to read through the file, without having to close and re-open it.
    """
    def __init__(self, filename):
        self.filehandle = False
        self.compressed = False
        self.filename = filename
        if filename.endswith(".gz"):
            self.filehandle = gzip.open(filename, 'r')
            self.compressed = True
        else:
            self.filehandle = open(filename, 'r')

        if not self.filehandle:
            raise NameError("Can't open %s for FastQ access" % (filename))

    def __iter__(self):
        if self.filehandle:
            self.filehandle.seek(0,0)
        return self

    def __next__(self):

# Read the first line which should be a sequence header

        header_line = self.filehandle.readline()
        if header_line == "" or header_line == b'':
            raise StopIteration
        header_line = header_line.strip()

# Deal with the bytestrings that come from gzip compressed files, and convert them to strings

        if self.compressed:
            header_line = header_line.decode('utf-8')

# Make sure the file is in the right format

        if not header_line.startswith("@"):
            raise FastqFormatError("Invalid header line")

# Get the sequence ID (without the @) and description, if there is one

        if " " in header_line:
            seqid, description = header_line.split(" ", 1)
        else:
            seqid = header_line
            description = None
        seqid = seqid[1:]

# Read the next line, which should contain the sequence itself

        sequence = self.filehandle.readline().strip()
        if self.compressed:
            sequence = sequence.decode('utf-8')

# Skip the line with the "+" that separates sequence and quality

        self.filehandle.readline()

# Get the quality string from the next line

        quality_string = self.filehandle.readline().strip()
        if self.compressed:
            quality_string = quality_string.decode('utf-8')
        return FastqSeq(seqid, description, sequence, quality_string)

    def __repr__(self):
        return "FastaReader(\"%s\")" % (self.filename)

class FastqSeq:
    def __init__(self, id=None, description=None, seq=None, quality_string=None):
        """ Class to hold the representation of a DNA sequence from a FASTQ file, containing the sequence itself plus the
        encoded quality string representing the quality score of each base in the sequence
        """
        self.id = id
        self.description = description
        self.seq = seq
        self.quality_string = quality_string

    def __repr__(self):
        return "FastqSeq(id=\"%s\", description=\"%s\", seq=\"%s\", quality_string=\"%s\")" % (self.id, self.description, self.seq, self.quality_string)

    def __str__(self):
        return self.seq

    def reverse_complement(self):
        """ Generates the reverse complement sequence of a DNA sequence 
        """
        revcomp = self.seq.translate(complementer)
        return revcomp[::-1]
    
    def to_fastq(self):
        """ Generates the text representation of the FASTQ sequence in a format that is suitable to be written to a FASTQ format text file
        """
        out = "@" + self.id
        if self.description:
            out += " " + self.description

        out += "\n" + self.seq + "\n+\n"
        out += self.quality_string

        return out

    def to_fasta(self):
        """ Generates the text representation of the FASTQ seqeunce in a format that can be written to a FASTA format file (which omits the quality information)
        """
        out = ">" + self.id
        if self.description:
            out += " " + self.description

        seq_lines = []
        for i in range(len(self.seq) // 60):
            seq_lines.append(self.seq[i*60:(i+1)*60])
        if len(self.seq) % 60 > 0:
            seq_lines.append(self.seq[(len(self.seq) // 60)*60:])

        out += "\n" + "\n".join(seq_lines)

        return out

if __name__ == "__main__":
    fr = FastqReader("short.fastq")
    for s in fr:
        print(s.to_fasta())
        # print(s.id, s.seq, s.quality_string)
