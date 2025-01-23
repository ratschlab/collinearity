from __future__ import annotations
__all__ = ['alignment_t', 'rt_index_t']

class alignment_t:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    @property
    def presence(self) -> float:
        ...
    @property
    def trg_name(self) -> str:
        ...
    @property
    def trg_pos(self) -> int:
        ...

class rt_index_t:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, batch_size: int, filename: str) -> None:
        """
        Create an index and the associated query client
        :param batch_size:
        :param filename:
        """
        ...
    def dump(self, basename: str) -> None:
        ...
    def query_batch(self, sequences: list[str]) -> list[alignment_t]:
        ...