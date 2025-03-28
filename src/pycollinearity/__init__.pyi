
class Alignment:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    def __init(self, header: str, fwd: bool, start: int, pres_frac: float, qry_len: int) -> None:
        ...
    @property
    def ctg(self) -> str:
        ...
    @property
    def r_st(self) -> int:
        ...
    @property
    def r_en(self) -> int:
        ...
    @property
    def strand(self) -> int:
        ...
    @property
    def pres_frac(self) -> float:
        ...

class Index:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
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
        ...
    def query(self, sequence: str) -> Alignment:
        ...
    def query_batch(self, sequences: list[str]) -> list[Alignment]:
        ...
