# Copied from ReadFish - TODO - insert citation + ack

from __future__ import annotations
from enum import Enum, unique
import attrs
from typing import Optional, Any, Sequence

@unique
class Decision(Enum):
    """Decision readfish has made about a read after Alignment"""

    #: The read aligned to a single location that is within a target region
    single_on: str = "single_on"
    #: The read aligned to a single location that is not in a target region
    single_off = "single_off"
    #: The read aligned to multiple locations, where at least one alignment is within a target region
    multi_on = "multi_on"
    #: The read aligned to multiple locations, none of which were in a target region
    multi_off = "multi_off"
    #: The read was basecalled but did not align
    no_map = "no_map"
    #: The read did not basecall
    no_seq = "no_seq"
    #: Too many signal chunks have been collected for this read
    above_max_chunks = "above_max_chunks"
    #: Fewer signal chunks for this read collected than required
    below_min_chunks = "below_min_chunks"
    #: Potential second half of a duplex read
    duplex_override = "duplex_override"
    #: Read sequenced as translocated portion was of unknown length at start of readfish
    first_read_override = "first_read_override"


@unique
class Action(Enum):
    """
    Action to take for a read.

    This enum class represents different actions that can be taken for a read during sequencing.
    Each action has a corresponding string value used for logging.

    :cvar unblock: Send an unblock command to the sequencer.
    :cvar stop_receiving: Allow the read to finish sequencing.
    :cvar proceed: Sample another chunk of data.

    :Example:

    Define an Action:

    >>> action = Action.unblock

    Access the string value of an Action:

    >>> action.value
    'unblock'
    """

    #: Send an unblock command to the sequencer
    unblock = "unblock"
    #: Allow the read to finish sequencing
    stop_receiving = "stop_receiving"
    #: Sample another chunk of data
    proceed = "proceed"


@attrs.define
class Result:
    """Result holder

    This should be progressively filled with data from the basecaller,
    barcoder, and then the aligner.

    :param channel: The channel that this read is being sequenced on
    :param read_id: The read ID assigned to this read by MinKNOW
    :param seq: The basecalled sequence for this read
    :param decision: The ``Decision`` that has been made, this will by used to determine the ``Action``
    :param barcode: The barcode that has been assigned to this read
    :param basecall_data: Any extra data that the basecaller may want to send to the aligner
    :param alignment_data: Any extra alignment data
    """

    channel: int
    read_id: str
    seq: str
    decision: Decision = attrs.field(default=Decision.no_seq)
    barcode: Optional[str] = attrs.field(default=None)
    basecall_data: Optional[Any] = attrs.field(default=None)
    # alignment_data: Optional[list[alignment_t]] = attrs.field(default=None)


def nice_join(seq: Sequence[Any], sep: str = ", ", conjunction: str = "or") -> str:
    """Join lists nicely

    :param seq: A sequence of objects that have a __str__ method.
    :param sep: The separator for the join, defaults to ", "
    :param conjunction: A conjunction between the joined list and the last element, defaults to "or"
    :return: The nicely joined string
    """
    seq = [str(x) for x in seq]

    if len(seq) <= 1 or conjunction is None:
        return sep.join(seq)
    else:
        return f"{sep.join(seq[:-1])} {conjunction} {seq[-1]}"
