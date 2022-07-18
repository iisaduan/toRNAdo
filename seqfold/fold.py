"""Predict nucleic acid secondary structure"""

import math
from .types import Energies


def _pair(s: str, i: int, i1: int, j: int, j1: int) -> str:
    """Return a stack representation, a key for the NN maps

    Args:
        s: Sequence being folded
        i: leftmost index
        i1: index to right of i
        j: rightmost index
        j1: index to left of j
    
    Returns:
        str: string representation of the pair
    """

    return (
        (s[i] if i >= 0 else ".")
        + (s[i1] if i1 >= 0 else ".")
        + "/"
        + (s[j] if j >= 0 else ".")
        + (s[j1] if j1 >= 0 else ".")
    )



def _d_g(d_h: float, d_s: float, temp: float) -> float:
    """Find the free energy given delta h, s and temp

    Args:
        d_h: The enthalpy increment in kcal / mol
        d_s: The entropy increment in cal / mol
        temp: The temperature in Kelvin

    Returns:
        The free energy increment in kcal / (mol x K)
    """

    return d_h - temp * (d_s / 1000.0)


def _j_s(query_len: int, known_len: int, d_g_x: float, temp: float) -> float:
    """Estimate the free energy of length query_len based on one of length known_len.

    The Jacobson-Stockmayer entry extrapolation formula is used
    for bulges, hairpins, etc that fall outside the 30nt upper limit
    for pre-calculated free-energies. See SantaLucia and Hicks (2004).

    Args:
        query_len: Length of element without known free energy value
        known_len: Length of element with known free energy value (d_g_x)
        d_g_x: The free energy of the element known_len
        temp: Temperature in Kelvin

    Returns:
        float: The free energy for a structure of length query_len
    """

    gas_constant = float(0.0019872)
    return d_g_x + 2.44 * gas_constant * temp * math.log(query_len / known_len)


def _stack(
    seq: str, i: int, i1: int, j: int, j1: int, temp: float, emap: Energies
) -> float:
    """Get the free energy for a stack.

    Using the indexes i and j, check whether it's at the end of
    the sequence or internal. Then check whether it's a match
    or mismatch, and return.

    Two edge-cases are terminal mismatches and dangling ends.
    The energy of a dangling end is added to the energy of a pair
    where i XOR j is at the sequence's end.

    Args:
        seq: The full folding sequence
        i: The start index on left side of the pair/stack
        i1: The index to the right of i
        j: The end index on right side of the pair/stack
        j1: The index to the left of j
        temp: Temperature in Kelvin

    Returns:
        float: The free energy of the NN pairing
    """

    if any(x >= len(seq) for x in [i, i1, j, j1]):
        return 0.0

    pair = _pair(seq, i, i1, j, j1)
    if any(x == -1 for x in [i, i1, j, j1]):
        # it's a dangling end
        d_h, d_s = emap.DE[pair]
        return _d_g(d_h, d_s, temp)

    if i > 0 and j < len(seq) - 1:
        # it's internal
        d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.INTERNAL_MM[pair]
        return _d_g(d_h, d_s, temp)

    if i == 0 and j == len(seq) - 1:
        # it's terminal
        d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.TERMINAL_MM[pair]
        return _d_g(d_h, d_s, temp)

    if i > 0 and j == len(seq) - 1:
        # it's dangling on left
        d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.TERMINAL_MM[pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = seq[i - 1] + seq[i] + "/." + seq[j]
        if pair_de in emap.DE:
            d_h, d_s = emap.DE[pair_de]
            d_g += _d_g(d_h, d_s, temp)
        return d_g

    if i == 0 and j < len(seq) - 1:
        # it's dangling on right
        d_h, d_s = emap.NN[pair] if pair in emap.NN else emap.TERMINAL_MM[pair]
        d_g = _d_g(d_h, d_s, temp)

        pair_de = "." + seq[i] + "/" + seq[j + 1] + seq[j]
        if pair_de in emap.DE:
            d_h, d_s = emap.DE[pair_de]
            d_g += _d_g(d_h, d_s, temp)
        return d_g

    return 0


def _hairpin(seq: str, i: int, j: int, temp: float, emap: Energies) -> float:
    """Calculate the free energy of a hairpin.

    Args:
        seq: The sequence we're folding
        i: The index of start of hairpin
        j: The index of end of hairpin
        temp: Temperature in Kelvin
        emap: Map of energies

    Returns:
        float: The free energy increment from the hairpin structure
    """

    if j - i < 4:
        return math.inf

    hairpin = seq[i : j + 1]
    hairpin_len = len(hairpin) - 2
    pair = _pair(seq, i, i + 1, j, j - 1)

    if emap.COMPLEMENT[hairpin[0]] != hairpin[-1]:
        # not known terminal pair, nothing to close "hairpin"
        raise RuntimeError

    d_g = 0.0
    if emap.TRI_TETRA_LOOPS and hairpin in emap.TRI_TETRA_LOOPS:
        # it's a pre-known hairpin with known value
        d_h, d_s = emap.TRI_TETRA_LOOPS[hairpin]
        d_g = _d_g(d_h, d_s, temp)

    # add penalty based on size
    if hairpin_len in emap.HAIRPIN_LOOPS:
        d_h, d_s = emap.HAIRPIN_LOOPS[hairpin_len]
        d_g += _d_g(d_h, d_s, temp)
    else:
        # it's too large, extrapolate
        d_h, d_s = emap.HAIRPIN_LOOPS[30]
        d_g_inc = _d_g(d_h, d_s, temp)
        d_g += _j_s(hairpin_len, 30, d_g_inc, temp)

    # add penalty for a terminal mismatch
    if hairpin_len > 3 and pair in emap.TERMINAL_MM:
        if pair in emap.TERMINAL_MM:
            d_h, d_s = emap.TERMINAL_MM[pair]
            d_g += _d_g(d_h, d_s, temp)

    # add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
    if hairpin_len == 3 and (hairpin[0] == "A" or hairpin[-1] == "A"):
        d_g += 0.5  # convert to entropy

    return d_g


def _bulge(
    seq: str, i: int, i1: int, j: int, j1: int, temp: float, emap: Energies
) -> float:
    """Calculate the free energy associated with a bulge.

        seq: The full folding DNA sequence
        i: The start index of the bulge
        i1: The index to the right of i
        j: The end index of the bulge
        j1: The index to the left of j
        loop: The sequence of the bulge
        temp: Temperature in Kelvin
        emap: Map to DNA/RNA energies

    Returns:
        float: The increment in free energy from the bulge
    """

    loop_len = max(i1 - i - 1, j - j1 - 1)
    if loop_len <= 0:
        raise RuntimeError

    # add penalty based on size
    if loop_len in emap.BULGE_LOOPS:
        d_h, d_s = emap.BULGE_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large for pre-calculated list, extrapolate
        d_h, d_s = emap.BULGE_LOOPS[30]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, 30, d_g, temp)

    if loop_len == 1:
        # if len 1, include the delta G of intervening NN (SantaLucia 2004)
        pair = _pair(seq, i, i1, j, j1)
        assert pair in emap.NN
        d_g += _stack(seq, i, i1, j, j1, temp, emap)

    # penalize AT terminal bonds
    if any(seq[k] == "A" for k in [i, i1, j, j1]):
        d_g += 0.5

    return d_g


def _internal_loop(
    seq: str, i: int, i1: int, j: int, j1: int, temp: float, emap: Energies
) -> float:
    """Calculate the free energy of an internal loop.

    The first and last bp of both left and right sequences
    are not themselves parts of the loop, but are the terminal
    bp on either side of it. They are needed for when there's
    a single internal looping bp (where just the mismatching
    free energies are used)

    Note that both left and right sequences are in 5' to 3' direction

    This is adapted from the "Internal Loops" section of SantaLucia/Hicks, 2004

    Args:
        seq: The sequence we're folding
        i: The index of the start of structure on left side
        i1: The index to the right of i
        j: The index of the end of structure on right side
        j1: The index to the left of j
        temp: Temperature in Kelvin
        emap: Dictionary mapping to energies for DNA/RNA

    Returns:
        float: The free energy associated with the internal loop
    """

    loop_left = i1 - i - 1
    loop_right = j - j1 - 1
    loop_len = loop_left + loop_right

    if loop_left < 0 or loop_right < 0:
        raise RuntimeError
    elif loop_left == 0 and loop_right == 0:
        return _stack(seq, i, i1, j, j1, temp, emap)
    elif loop_left == 0 or loop_right == 0:
        return _bulge(seq, i, i1, j, j1, temp, emap)

    # single bp mismatch, sum up the two single mismatch pairs
    if loop_left == 1 and loop_right == 1:
        mm_left = _stack(seq, i, i1, j, j1, temp, emap)
        mm_right = _stack(seq, i1 - 1, i1, j1 + 1, j1, temp, emap)
        return mm_left + mm_right

    # apply a penalty based on loop size
    if loop_len in emap.INTERNAL_LOOPS:
        d_h, d_s = emap.INTERNAL_LOOPS[loop_len]
        d_g = _d_g(d_h, d_s, temp)
    else:
        # it's too large an internal loop, extrapolate
        d_h, d_s = emap.INTERNAL_LOOPS[30]
        d_g = _d_g(d_h, d_s, temp)
        d_g = _j_s(loop_len, 30, d_g, temp)

    # apply an asymmetry penalty
    loop_asymmetry = abs(loop_left - loop_right)
    d_g += 0.3 * loop_asymmetry

    # apply penalty based on the mismatching pairs on either side of the loop
    pair_left_mm = _pair(seq, i, i + 1, j, j - 1)
    d_h, d_s = emap.TERMINAL_MM[pair_left_mm]
    d_g += _d_g(d_h, d_s, temp)

    pair_right_mm = _pair(seq, i1 - 1, i1, j1 + 1, j1)
    d_h, d_s = emap.TERMINAL_MM[pair_right_mm]
    d_g += _d_g(d_h, d_s, temp)

    return d_g