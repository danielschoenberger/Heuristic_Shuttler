from get_sequence import parse_qasm

n = 6
filename = "QASM_files/qft_%squbits.qasm" % n


def get_gatelist(filename):
    parsed_file = parse_qasm(filename=filename)
    gatelist = {}
    for gate in parsed_file:
        for qubit in gate:
            if qubit not in gatelist:
                gatelist[qubit] = []
            gatelist[qubit].append(gate)
    return gatelist


def get_candidates(gatelist):
    candidates = {}
    for qubit in gatelist:
        candidates[qubit] = []
        for i, gate in enumerate(gatelist[qubit]):
            if all(commute(gatelist[qubit][j], gate) for j in range(i)):
                candidates[qubit].append(gate)
    return candidates


def commute(gate0, gate1):
    # if name is equal
    return bool(any((len(gate0) == 1, len(gate1) == 1)))


def get_frontlayer():
    pass


gatelist = get_gatelist(filename)
print(gatelist)
candidates = get_candidates(gatelist)
print("\n", candidates, "\n")
