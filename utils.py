from IBMQuantumExperience import IBMQuantumExperience
from qiskit import IBMQ, BasicAer, execute
from qiskit.providers.ibmq import least_busy
from qiskit.tools.monitor import job_monitor, backend_overview

import numpy as np
from matplotlib import pyplot as plt


IBMQ.load_accounts()
print("Account loaded")


def launch(circ, shots=1024, backend_type="local", state_vectors=False, max_credits=3, verbose=False):
    if verbose:
        print("Depth : ", circ.depth())
        print("Qbits : ", circ.width())
        print("Number of gates : ", circ.count_ops())

    if backend_type == "local":
        if state_vectors:
            simulator = 'statevector_simulator'
        else:
            simulator = 'qasm_simulator'
        backend = BasicAer.get_backend(simulator)
    elif backend_type in ["q", "quantum"]:
        min_qbits = circ.width()
        print("Available backends:")
        IBMQ.backends()

        large_enough_devices = IBMQ.backends(
            filters=lambda x: x.configuration().n_qubits > min_qbits and not x.configuration().simulator)
        backend = least_busy(large_enough_devices)
        print("The best backend is " + backend.name())
    elif backend_type in ["hpc", "ibm_sim"]:
        backend = IBMQ.backends(filters=lambda x: x.configuration().simulator)[0]
    else:
        print("Invalid backend name. Switching to local simulator")
        backend = BasicAer.get_backend('qasm_simulator')

    job = execute(circ, backend=backend, shots=shots, max_credits=max_credits)

    if backend_type in ["q", "quantum"]:
        job_monitor(job)
    result = job.result()

    counts = result.get_counts(circ)
    # print(counts)
    return counts


def backends_info():
    backend_overview()


def show_credits():
    api = IBMQuantumExperience(IBMQ.stored_accounts()[0]['token'])
    print(api.get_my_credits())


def plot_hist(data, crop=False):
    names = [""]*len(data)
    val = [0]*len(data)
    maxi = max(data.values())
    for ind, key in enumerate(data.keys()):
        if not crop or data[key] > maxi/10:
            names[ind] = key

    names.sort()
    for ind, key in enumerate(names):
        val[ind] = data[key]

    ind = np.arange(len(data))  # the x locations for the groups
    width = 0.9  # the width of the bars

    fig, ax = plt.subplots()
    rects = ax.bar(ind, val, width, color='b')

    # add some text for labels, title and axes ticks
    ax.set_xticks(ind)
    ax.set_xticklabels(names)

    def autolabel(rectangles):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rectangles:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height,
                    '%d' % int(height),
                    ha='center', va='bottom')

    autolabel(rects)

    plt.show()


def cnx_o_paquet(c, r, cont, target, losts):
    """Adds a len(cont)-qbits-Control(X) gate to c, on r, controlled by the List cont at qbit target,
    losing one qbit in the process"""
    # c.barrier(r)
    # print(len(cont), len(losts))
    if len(cont) == 1:
        c.cx(cont[0], target)
    elif len(cont) == 2:
        c.ccx(cont[0], cont[1], target)
    elif len(cont) == 3:  # and len(losts) == 1:
        c.ccx(cont[2], cont[1], losts[0])
        c.ccx(cont[0], losts[0], target)
        c.ccx(cont[2], cont[1], losts[0])
        c.ccx(cont[0], losts[0], target)
    # elif len(cont) == 4 and len(losts) == 1:
    #    cnx_o(c, r, cont, target, losts[0])
    # elif len(cont) == 5:
    #    cnx_o(c, r, cont, target, losts[0])
    else:
        m = len(cont)
        # print(m)
        # c.barrier(r)
        c.ccx(cont[0], losts[0], target)
        # c.barrier(r)
        # c.h(r)

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
    """Adds a len(cont)-qbits-Control(X) gate to c, on r, controlled by the List cont at qbit target,
    losing one qbit in the process"""
    #c.barrier(r)
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
        # print(m, len(cont)+2, len(cont), m1)
        # A more efficient way to do this would be defining a new circuit
        # and just multiply it istead of doing the same thing twice
        cnx_o_paquet(c, r, cont[m1:], lost, [target] + cont[:m1])
        cnx_o_paquet(c, r, [lost] + cont[:m1], target, cont[m1:])
        cnx_o_paquet(c, r, cont[m1:], lost, [target] + cont[:m1])
        cnx_o_paquet(c, r, [lost] + cont[:m1], target, cont[m1:])
    return 1