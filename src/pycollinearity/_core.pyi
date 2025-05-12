from __future__ import annotations
import typing
__all__ = ['Alignment', 'DynIndex', 'Index', 'Request', 'Response', 'ResponseGenerator']


class Alignment:
    """
    An alignment object that contains the reference name, strand, start and end positions, and fraction of k-mer matches.
    """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, header: str, fwd: bool, start: int, pres_frac: float, qry_len: int) -> None:
        """
        Creates an alignment object
        :param header: name of the reference to which the query was aligned
        :param fwd: forward or reverse strand of the reference
        :param start: start position in the reference
        :param qry_len: length of the query
        :param pres_frac: fraction of k-mers in the query that are also present in reference[start: start + qry_len]
        """
        ...
    @property
    def ctg(self) -> str:
        """
        :return: name of the reference to which the query was aligned. A '*' means that an alignment couldn't be found
        """
        ...
    @property
    def pres_frac(self) -> float:
        """
        :return: fraction of k-mers in the query that are also present in reference[start: start + qry_len]
        """
        ...
    @property
    def r_en(self) -> int:
        """
        :return: end position of the alignment in the reference
        """
        ...
    @property
    def r_st(self) -> int:
        """
        :return: end position of the alignment in the reference
        """
        ...
    @property
    def strand(self) -> int:
        """
        :return: forward or reverse strand of the reference
        """
        ...


class DynIndex:
    def __init__(self, *args, **kwargs) -> None:
        ...
    def add(self, arg0: str, arg1: str) -> None:
        ...
    def add_batch(self, arg0: list, arg1: list) -> None:
        ...
    def merge(self) -> None:
        ...
    def query(self, arg0: str) -> Alignment:
        ...
    def query_batch(self, arg0: list) -> list[Alignment]:
        ...


class Index:
    def __init__(self, input: str, *args, **kwargs) -> None:
        """
        Create an index
        :param input: input fasta file
        :param args: arguments
        :param kwargs: keyword arguments
        :keyword jaccard : Use jaccard similarity. [implicit: "true", default: false]
        :keyword fr : Index both forward and reverse strands of the reference. [implicit: "true", default: false]
        :keyword k: k-mer length [default: 15]
        :keyword pf : Fraction of k-mers that must be present in an alignment. [default: 0.1]
        :keyword bw : Width of the band in which kmers contained will be considered collinear [default: 15]
        :keyword jc-frag-len : If jaccard is set, the sequence are indexed and queried in overlapping fragments of this length. [default: 180]
        :keyword jc-frag-ovlp-len : If jaccard is set, the sequence are indexed and queried in fragments which overlap this much. [default: 120]
        """
        ...
    def dump(self, basename: str) -> None:
        """
        Dump index to file
        :param basename: basename of the output file. A .cidx extension will be added to the basename
        :return: None
        """
        ...
    def load(self, basename: str) -> None:
        """
        Load index from file
        :param basename: basename of the input file. A .cidx extension will be added to the basename
        :return: None
        """
        ...
    def query(self, sequence: str) -> Alignment:
        """
        Query a sequence in the index
        :param sequence: query sequence
        :return: an alignment of the query
        """
        ...
    def query_batch(self, sequences: list[str]) -> list[Alignment]:
        """
        Query a batch of sequences in the index using multiple threads
        :param sequences: list of query sequences
        :return: a list of alignments
        """
        ...
    def query_stream(self, requests: typing.Iterator) -> ResponseGenerator:
        """
        Query a stream of requests and return a stream of responses (Readfish compatible)
        :param requests: an iterator of query request objects
        :return: A generator of query response objects
        """
        ...


class Request:
    """
    A query request object with channel, query id, and query sequence
    """
    channel: int
    id: str
    seq: str
    def __init__(self, channel: int, id: str, seq: str) -> None:
        ...


class Response:
    """
    A query response object with channel, query id, and its alignment in the index
    """
    alignment: Alignment
    channel: int
    id: str
    def __init__(self, channel: int, id: str, alignment: Alignment) -> None:
        ...


class ResponseGenerator:
    """
    A generator for query responses
    """
    def __iter__(self) -> ResponseGenerator:
        ...
    def __next__(self) -> Response:
        ...
