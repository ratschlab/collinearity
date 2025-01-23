# Adapted from ReadFish - TODO - insert citation + ack

import os
from collections import namedtuple
from pathlib import Path

from packaging.version import Version, parse
from minknow_api.manager import Manager, FlowCellPosition
from minknow_api import Connection
from basecaller import BaseCaller
from config import Config
from prelude import *
from _core import *
from read_cache import AccumulatingCache
from readuntil_client import ReadUntilClient

__all__ = ['Pipeline']

from utils import Result


def _get_minknow_version(host: str = "127.0.0.1", port: int = None) -> Version:
    """
    Get the version of MinKNOW

    :param host: The host the RPC is listening on, defaults to "127.0.0.1"
    :param port: The port the RPC is listening on, defaults to None

    :return: The version of MinKNOW readfish is connected to
    """
    manager = Manager(host=host, port=port)
    minknow_version = parse(manager.core_version)
    return minknow_version

def get_device(
    device: str, host: str = "127.0.0.1", port: int = None
) -> FlowCellPosition:
    """Get a position for a specific device over the minknow API

    :param device: The device name - example X1 or MS00000
    :param host: The host the RPC is listening on, defaults to "127.0.0.1"
    :param port: The port the RPC is listening on, defaults to None
    :raises ValueError: If their is no match on any of the positions for the given device name
    :return: The position representation from the MinkKNOW API
    """
    manager = Manager(host=host, port=port)
    for position in manager.flow_cell_positions():
        if position.name == device:
            return position
    raise ValueError(f"Could not find device {device!r}")


class Pipeline:
    def __init__(self):
        self.config = Config()
        mk_addr = self.config.minknow_addr.split(':')
        host, port = mk_addr[0], int(mk_addr[1])

        minknow_version = _get_minknow_version(host=host, port=port)
        if str(minknow_version) == 'fake_server': minknow_version = Version("6.0.0")
        print("Minknow version: ", minknow_version)

        # Fetch sequencing device
        position = get_device(self.config.device, host=host, port=port)

        # Create a read until client
        self.read_until_client = ReadUntilClient(
            mk_host=position.host,
            mk_port=position.description.rpc_ports.secure,
            filter_strands=True,
            cache_type=AccumulatingCache,
            prefilter_classes={
                "strand",
                "adapter",
            },
        )

        info("Building index..")
        self.index = rt_index_t(self.read_until_client.channel_count, self.config.ref)

        info("Connecting to basecaller..")
        self.run_information = self.read_until_client.connection.protocol.get_run_info()
        self.sample_rate = self.read_until_client.connection.device.get_sample_rate().sample_rate

        self.basecaller = BaseCaller(run_information=self.run_information, sample_rate=self.sample_rate)


    def run(self):
        # start the client running
        self.read_until_client.run(
            first_channel=1,
            last_channel=self.read_until_client.channel_count,
            max_unblock_read_length_seconds=5,
            accepted_first_chunk_classifications=[
                "strand",
                "adapter",
            ],
        )

        # begin main loop
        try:
            while self.read_until_client.is_running:
                chunks = self.read_until_client.get_read_chunks(self.read_until_client.channel_count, last=True)
                calls = self.basecaller.basecall(
                    chunks, self.read_until_client.signal_dtype, self.read_until_client.calibration_values
                )
                alignments = self.index.query_batch([x.seq for x in calls])
                self.act(calls, alignments)

        except KeyboardInterrupt:
            info("Keyboard interrupt received, stopping readfish.")
            pass
        finally:
            self.read_until_client.reset()

    def act(self, sequences: list[Result], alignments: list[alignment_t]):
        # todo - this is very crude at the moment. we need to take into account the number of chunks consumed before acting
        if len(sequences) != len(alignments):
            error("Number of sequences and alignments do not match")
        present = []
        absent = []
        for sequence, alignment in zip(sequences, alignments):
            if alignment.trg_pos < 0:
                absent.append((sequence.channel, sequence.read_id))
            else:
                present.append((sequence.channel, sequence.read_id))
        if self.config.background:  # these sequences are present in the background and will be ejected
            self.read_until_client.unblock_read_batch(present)
            self.read_until_client.stop_receiving_batch(absent)
        else:   # these sequences are present in the foreground so
            self.read_until_client.stop_receiving_batch(present)
            self.read_until_client.unblock_read_batch(absent)