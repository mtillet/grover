"""
Code by Raphael Lambert
"""

import numpy as np
from matplotlib import pyplot as plt


def cnx_o_paquet(c, cont, target, losts):
    """
    Adds a len(cont)-qbits-Control(X) gate to c, on r, controlled by the List cont at qbit target,
    losing one qbit in the process
    """
    # c.barrier(r)

    if len(cont) == 1:
        c.cx(cont[0], target)
    elif len(cont) == 2:
        c.ccx(cont[0], cont[1], target)
    elif len(cont) == 3:
        c.ccx(cont[2], cont[1], losts[0])
        c.ccx(cont[0], losts[0], target)
        c.ccx(cont[2], cont[1], losts[0])
        c.ccx(cont[0], losts[0], target)

    else:
        m = len(cont)
        # c.barrier(r)
        c.ccx(cont[0], losts[0], target)
        # c.barrier(r)

        for i in range(1, m - 2):
            c.ccx(cont[i], losts[i], losts[i - 1])
        # c.barrier(r)
        c.ccx(cont[-1], cont[-2], losts[m - 3])
        # c.barrier(r)
        for i in range(m - 3, 0, -1):
            c.ccx(cont[i], losts[i], losts[i - 1])
        # c.barrier(r)
        c.ccx(cont[0], losts[0], target)
        # c.barrier(r)
        for i in range(1, m - 2):
            c.ccx(cont[i], losts[i], losts[i - 1])
        # c.barrier(r)
        c.ccx(cont[-1], cont[-2], losts[m - 3])
        # c.barrier(r)
        for i in range(m - 3, 0, -1):
            c.ccx(cont[i], losts[i], losts[i - 1])
    return 1


def cnx_o(c, r, cont, target, lost):
    """
    Adds a len(cont)-qbits-Control(X) gate to c, on r, controlled by the List cont at qbit target,
    losing one qbit in the process

    :param c: Quantum circuit
    :param r: Quantum Registers
    :param cont: control qbits list
    :param target: target qbit for the gate
    :param lost: ancillary qubit used to accelerate computation

    :return:
    """

    if len(cont) == 1:
        c.cx(cont[0], target)
    elif len(cont) == 2:
        c.ccx(cont[0], cont[1], target)
    elif len(cont) == 3:
        c.ccx(cont[2], cont[1], lost)
        c.ccx(cont[0], lost, target)
        c.ccx(cont[2], cont[1], lost)
        c.ccx(cont[0], lost, target)
    else:
        m = int(np.ceil(len(cont) / 2 + 1))
        m1 = len(cont) - m

        # A more efficient way to do this would be defining a new circuit
        # and just multiply it instead of doing the same thing twice
        cnx_o_paquet(c, cont[m1:], lost, [target] + cont[:m1])
        cnx_o_paquet(c, [lost] + cont[:m1], target, cont[m1:])
        cnx_o_paquet(c, cont[m1:], lost, [target] + cont[:m1])
        cnx_o_paquet(c, [lost] + cont[:m1], target, cont[m1:])

    return 1
