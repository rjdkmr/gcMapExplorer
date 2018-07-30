#!/usr/bin/env python3
#
# This file is part of gcMapExplorer
# Copyright (C) 2016-2018  Rajendra Kumar, Ludvig Lizana, Per Stenberg
#
# gcMapExplorer is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gcMapExplorer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gcMapExplorer.  If not, see <http://www.gnu.org/licenses/>.
#
# =============================================================================
"""
Fast parser for Hi-C contact map format described in:
Juicer Provides a One-Click System for Analyzing Loop-Resolution Hi-C Experiments
Durand, Neva C. et al.
Cell Systems, Volume 3, Issue 1, p. 95-98
"""
import zlib
from argparse import ArgumentParser, ArgumentTypeError, FileType
from array import array
from collections import namedtuple
from enum import Enum
from functools import lru_cache
from struct import unpack, unpack_from, Struct


def _read_string(file):
    buf = b""
    for b in iter(lambda: file.read(1), b"\0"):
        if b == "":
            raise EOFError("Buffer unexpectedly empty while trying to read null-terminated string")
        buf += b
    return buf.decode("utf-8")


_struct_int = Struct("<i")


def _read_int(file):
    return _struct_int.unpack(file.read(4))[0]


Chromosome = namedtuple("Chromosome", "length index")
Record = namedtuple("Record", "position size")
ScaleFactor = namedtuple("ScaleFactor", "index scale_factor")
ExpectedValue = namedtuple("ExpectedValue", "bin_size unit values scale_factors")
NormExpectedValue = namedtuple("NormExpectedValue", "type bin_size unit values scale_factors")
NormVector = namedtuple("NormVector", "type index unit bin_size position n_bytes")


class Unit(Enum):
    BP = "BP"
    FRAG = "FRAG"


class NormType(Enum):
    VC = "VC"
    VC_SQRT = "VC_SQRT"
    KR = "KR"


class BlockReader:
    """The BlockReader is an iterable that decompresses block data on the fly
    which allows for huge files to be read.

    This class is not supposed to be used directly, it is supposed to be instantiated
    by the HicParser.

    For hic version 7 and greater.
    """
    __struct_hh = Struct("<hh")

    def __init__(self, file, unit, bin_size, bin_count, column_count, blocks):
        """Initialize a BlockReader

        Parameters
        ----------
        file
            hic file object (must be open)
        unit : Unit
            the resolution unit
        bin_size : int
            the bin size
        bin_count : int
            the size of each block in bins
        column_count : int
            the number of columns for the block grid
        blocks : iterable tuple
            id, position, size
        """
        self.file = file
        self.unit = unit
        self.bin_size = bin_size
        self.bin_count = bin_count
        self.column_count = column_count
        self.blocks = blocks

    def __iter__(self):
        """Generator which yields all block data.

        Yields
        ------
        int, int, int
            bin_x, bin_y, count
        """
        for __, position, size in self.blocks:
            self.file.seek(position)
            compressed = self.file.read(size)
            data = zlib.decompress(compressed)

            bin_x_offset, bin_y_offset, use_float, block_type = unpack_from("<ii?b", data, 4)
            offset = 14

            if block_type == 1:
                row_count = unpack_from("<h", data, offset)[0]
                offset += 2
                for _ in range(row_count):
                    y, col_count = self.__struct_hh.unpack_from(data, offset)
                    offset += 4
                    bin_y = y + bin_y_offset
                    if use_float:
                        fmt = "<{}".format("hf" * col_count)
                        col = unpack_from(fmt, data, offset)  # x, count, x, count, ...
                        offset += 6 * col_count
                    else:
                        col = array("h", data[offset:offset + 4 * col_count])  # x, count, x, count, ...
                        offset += 4 * col_count
                    for x, count in zip(*[iter(col)] * 2):
                        yield x + bin_x_offset, bin_y, float(count)

            elif block_type == 2:
                n_pts, w = unpack_from("<ih", data, offset)
                offset += 6
                for i in range(n_pts):
                    row = int(i / w)
                    col = i - row * w
                    if use_float:
                        count = unpack_from("<f", data, offset)[0]
                        offset += 4
                    else:
                        count = float(unpack_from("<h", data, offset)[0])
                        offset += 2
                    yield bin_x_offset + col, bin_y_offset + row, count


class BlockReaderV6(BlockReader):
    """Specialized version of BlockReader for hic version 6.

    """

    def __iter__(self):
        """Generator which yields all block data.

        Yields
        ------
        int, int, float
            bin_x, bin_y, count
        """
        for __, position, size in self.blocks:
            self.file.seek(position)
            compressed = self.file.read(size)
            data = zlib.decompress(compressed)
            n_values = unpack_from("<i", data)[0]
            values = unpack_from("<{}".format("iif" * n_values), data, 4)  # x, y, count, x, y, count, ...
            for bin_x, bin_y, count in zip(*[iter(values)] * 3):
                yield bin_x, bin_y, count


class HicParser:
    """Fast access to hic data.

    Supports version 6 and greater.

    Examples
    --------
    >>> with open("somefile.hic", "rb") as f:
    >>>     hic = HicParser(f)
    >>>     hic.bp_resolutions
    [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000]
    >>>     record = hic.record("X", "X")
    >>>     blocks = hic.blocks(record, 5000)
    >>>     for x, y, count in blocks:
    >>>         print(x * blocks.bin_size, y * blocks.bin_size, count)
    0 10000 1.0
    10000 10000 1.0
    ...

    Attributes
    ----------
    name: str
        the file name
    version : int
        the format version
    genome_id : str
        Genome identifier or file that contains list of genome identifiers
    attributes : dict of str to str
        Dictionary of metadata that describes the experiment
    chromosomes : dict of str -> Chromosome
        Chromosome dictionary
    bp_resolutions : array of int
        Base pair resolutions
    frag_resolutions : array of int
        Fragment resolutions
    sites : array of int
        Restriction sites
    records : dict of str to Record
        Records of chromosome pairs
    expected_value_vectors : list of ExpectedValue
    norm_expected_value_vectors : list of NormExpectedValue
    norm_vectors : list of NormVector
    """

    def __init__(self, file):
        """Initialize a HicParser

        Parameters
        ----------
        file
            hic file object (must be open as "rb")
        """
        self.file = file
        self.name = file.name

        if file.mode != "rb":
            raise ValueError("File mode must be 'rb'")

        # Read header
        file.seek(0)

        magic_string = unpack("<3sx", file.read(4))[0]

        if magic_string != b"HIC":
            raise Exception("No magic string found!")

        self.version = _read_int(file)

        if self.version < 6:
            raise Exception("Version {} no longer supported".format(self.version))

        master_index_pos = unpack("<q", file.read(8))[0]

        self.genome_id = _read_string(file)

        self.attributes = dict()
        n_attributes = _read_int(file)
        for _ in range(n_attributes):
            key = _read_string(file)
            value = _read_string(file)
            self.attributes[key] = value

        self.chromosomes = dict()
        n_chromosomes = _read_int(file)
        for index in range(n_chromosomes):
            chr_name = _read_string(file)
            chr_length = _read_int(file)
            self.chromosomes[chr_name] = Chromosome(chr_length, index)

        n_bp_resolutions = _read_int(file)
        self.bp_resolutions = array("i", file.read(4 * n_bp_resolutions))

        n_frag_resolutions = _read_int(file)
        self.frag_resolutions = array("i", file.read(4 * n_frag_resolutions))

        self.sites = []
        if n_frag_resolutions > 0:
            n_sites = _read_int(file)
            self.sites = array("i", file.read(4 * n_sites))

        # Read footer
        file.seek(master_index_pos)

        n_bytes_v5 = _read_int(file)

        self.records = dict()
        n_entries = _read_int(file)
        for _ in range(n_entries):
            key = _read_string(file)  # E.g. 1_1, 2_2, ...
            position, size = unpack("<qi", file.read(12))
            self.records[key] = Record(position, size)

        self.expected_value_vectors = []
        n_expected_value_vectors = _read_int(file)
        for _ in range(n_expected_value_vectors):
            unit = Unit(_read_string(file))
            bin_size, n_values = unpack("<2i", file.read(8))
            values = array("d", file.read(8 * n_values))
            n_chr_scale_factors = _read_int(file)
            scale_factors = [ScaleFactor(*unpack("<id", file.read(12))) for _ in range(n_chr_scale_factors)]
            self.expected_value_vectors.append(ExpectedValue(bin_size, unit, values, scale_factors))

        self.norm_expected_value_vectors = []
        n_norm_expected_value_vectors = _read_int(file)
        for _ in range(n_norm_expected_value_vectors):
            norm_type = NormType(_read_string(file))
            unit = Unit(_read_string(file))
            bin_size, n_values = unpack("<2i", file.read(8))
            values = array("d", file.read(8 * n_values))
            n_chr_scale_factors = _read_int(file)
            scale_factors = [ScaleFactor(*unpack("<id", file.read(12))) for _ in range(n_chr_scale_factors)]
            self.norm_expected_value_vectors.append(NormExpectedValue(norm_type, bin_size, unit, values, scale_factors))

        self.norm_vectors = []
        n_norm_vectors = _read_int(file)
        for _ in range(n_norm_vectors):
            norm_type = NormType(_read_string(file))
            chr_index = _read_int(file)
            unit = Unit(_read_string(file))
            bin_size, position, n_bytes = unpack("<iqi", file.read(16))
            self.norm_vectors.append(NormVector(norm_type, chr_index, unit, bin_size, position, n_bytes))

    @lru_cache(maxsize=32)
    def read_header(self, record):
        """Read all block headers for chromosome pairs defined by the record.

        Parameters
        ----------
        record : Record
            the chromosome pair record

        Returns
        -------
        list of BlockReader
        """
        file = self.file
        file.seek(record.position)

        chr1_idx, chr2_idx, n_resolutions = unpack("<3i", file.read(12))

        block_readers = []
        block_reader = BlockReaderV6 if self.version < 7 else BlockReader

        for _ in range(n_resolutions):
            unit = Unit(_read_string(file))

            resolution_index = _read_int(file)

            # currently unused
            sum_counts, occupied_cell_count, std_dev, percent95 = unpack("<4f", file.read(16))

            bin_size, block_bin_count, block_column_count, block_count = unpack("<4i", file.read(16))

            # id, position, size
            blocks = [unpack("<iqi", file.read(16)) for _ in range(block_count)]

            block_readers.append(block_reader(file, unit, bin_size, block_bin_count, block_column_count, blocks))

        return block_readers

    @lru_cache(maxsize=32)
    def norm_vector(self, chromosome, norm_type, bin_size, unit=Unit.BP):
        """Get norm vector for chromosome.

        Parameters
        ----------
        chromosome : str
            the chromosome name
        norm_type : NormType
            the norm type
        bin_size : int
            the bin size
        unit : Unit, optional
            the resolution unit

        Returns
        -------
        array of floats

        Raises
        ------
        LookupError
            If no norm vector was found matching the arguments.

        Examples
        --------
        >>> hic = HicParser(f)  # f: file object
        >>> c1, c2, res = "X", "Y", 5000
        >>> c1_norm = hic.norm_vector(c1, NormType.VC, res)
        >>> c2_norm = hic.norm_vector(c2, NormType.VC, res)
        >>> record = hic.record(c1, c2)
        >>> blocks = hic.blocks(record, res)
        >>> for bin_x, bin_y, count in blocks:
        >>>     print(count / (c1_norm[bin_x] * c2_norm[bin_y]))  # print normalized count
        """
        index = self.chromosome_index[chromosome]
        norm = self._get_norm_vector(index, norm_type, bin_size, unit)
        return self._read_norm_vector(norm.position)

    def _get_norm_vector(self, index, norm_type, bin_size, unit):
        try:
            return next(n for n in self.norm_vectors
                        if n.index == index and n.type == norm_type and n.unit == unit and n.bin_size == bin_size)
        except StopIteration:
            raise LookupError(
                "No norm vector found for chromosome {}, norm type {}, unit {}, bin size {}".format(
                    self.chromosome_name[index],
                    norm_type.name,
                    unit.name,
                    bin_size))

    def _read_norm_vector(self, position):
        self.file.seek(position)
        n_values = _read_int(self.file)
        return array("d", self.file.read(8 * n_values))

    @property
    def chromosome_index(self):
        """Chromosome index by name

        Returns
        -------
        int
            chromosome index
        """
        return {name: index for name, (__, index) in self.chromosomes.items()}

    @property
    def chromosome_name(self):
        """Chromosome name by index

        Returns
        -------
        str
            chromosome name
        """
        return {index: name for name, (__, index) in self.chromosomes.items()}

    def chromosome_names(self, record_key):
        """Get chromosome names from record key

        Parameters
        ----------
        record_key : str
            The key in the record dictionary corresponding to the chromosome pair.

        Returns
        -------
        tuple of str, str
            Two chromosome names

        Raises
        ------
        LookupError
            If chromosome names were not found.
        """
        index1, index2 = map(int, record_key.split("_"))

        try:
            return self.chromosome_name[index1], self.chromosome_name[index2]
        except KeyError:
            raise LookupError("No chromosomes found for record {}".format(record_key))

    def record(self, chromosome1, chromosome2):
        """Get record for chromosome pair.

        Parameters
        ----------
        chromosome1 : str
            chromosome name
        chromosome2 : str
            chromosome name

        Returns
        -------
        Record
            The record for the chromosome pair

        Raises
        ------
        LookupError
            If no record was found for the chromosome pair.
        """
        index1, index2 = sorted((self.chromosome_index[chromosome1], self.chromosome_index[chromosome2]))
        key = "{}_{}".format(index1, index2)

        try:
            return self.records[key]
        except KeyError:
            raise LookupError("No record found for chromosomes {} {}".format(chromosome1, chromosome2))

    def blocks(self, record, bin_size, unit=Unit.BP):
        """Get BlockReader for chromosome pair.

        Parameters
        ----------
        record : Record
            record for the chromosome pair
        bin_size : int
            bin size
        unit : Unit, optional
            resolution unit

        Returns
        -------
        BlockReader

        Raises
        ------
        LookupError
            If no block header was found matching the arguments.

        Examples
        --------
        >>> hic = HicParser(f)  # f: file object
        >>> record = hic.record("X", "Y")
        >>> blocks = hic.blocks(record, 5000)
        >>> for bin_x, bin_y, count in blocks:
        >>>     x = bin_x * blocks.bin_size
        >>>     y = bin_y * blocks.bin_size
        >>>     print(x, y, count)
        """
        try:
            return next(h for h in self.read_header(record) if h.unit == unit and h.bin_size == bin_size)
        except StopIteration:
            raise LookupError(
                "No block header found for unit {}, bin size {}".format(unit.name, bin_size))


class HicFileType(FileType):
    """Factory for creating hic file types for use with ArgumentParser

    See Also
    --------
    argparse.FileType
    """

    def __init__(self):
        super().__init__(mode="rb")

    def __call__(self, string):
        f = super().__call__(string)
        try:
            return HicParser(f)
        except Exception as e:
            raise ArgumentTypeError(e)


if __name__ == "__main__":
    parser = ArgumentParser(description="Fast hic file parser")
    parser.add_argument("input", type=HicFileType(), help="hic input file")
    args = parser.parse_args()

    hic = args.input
    print("Hic version {}\n".format(hic.version))

    if len(hic.attributes):
        print("Attributes:")
        for key, value in hic.attributes.items():
            print("{}: {}".format(key, value))

    if len(hic.chromosomes):
        print("Chromosomes")
        print(", ".join(hic.chromosomes), end="\n\n")

    if len(hic.bp_resolutions):
        print("Base pair resolutions:")
        print(", ".join(map(str, hic.bp_resolutions)), end="\n\n")

    if len(hic.frag_resolutions):
        print("Fragment resolutions:")
        print(", ".join(map(str, hic.frag_resolutions)), end="\n\n")

    if len(hic.sites):
        print("Sites:")
        print(", ".join(map(str, hic.sites)), end="\n\n")
