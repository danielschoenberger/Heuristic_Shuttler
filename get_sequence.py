import math
import re
from pathlib import Path


def generate_qft_qasm(n):
    """Generate QASM code for Quantum Fourier Transform for n qubits."""
    qasm_code = []

    # Header
    qasm_code.append("OPENQASM 2.0;")
    qasm_code.append('include "qelib1.inc";')
    qasm_code.append(f"qreg q[{n}];")
    qasm_code.append(f"creg c[{n}];")
    qasm_code.append("")

    # QFT circuit
    for target_qubit in range(n):
        qasm_code.append(f"h q[{target_qubit}];")
        for control_qubit in range(target_qubit + 1, n):
            angle = 2 * math.pi / (2 ** (control_qubit - target_qubit + 1))
            qasm_code.append(f"cu1({angle}) q[{control_qubit}],q[{target_qubit}];")

    return "\n".join(qasm_code)


def write_qft_to_file(n, filename="qft.qasm"):
    """Write QFT QASM code for n qubits to a file."""
    qasm_code = generate_qft_qasm(n)
    with Path(filename).open(mode="w") as file:
        file.write(qasm_code)


def extract_qubits_from_gate(gate_line):
    """Extract qubit numbers from a gate operation line."""
    # Regular expression to match qubits (assuming they are in the format q[<number>])
    pattern = re.compile(r"q\[(\d+)\]")
    matches = pattern.findall(gate_line)

    # Convert matched qubit numbers to integers
    return [int(match) for match in matches]


def parse_qasm(filename):
    """Parse a QASM file and return qubits used for each gate, preserving their order."""
    gates_and_qubits = []

    with Path(filename).open() as file:
        for line_ in file:
            line = line_.strip()

            # Check if line represents a gate operation
            if not line.startswith(("OPENQASM", "include", "qreg", "creg", "gate")):
                qubits = extract_qubits_from_gate(line)
                if qubits:
                    gates_and_qubits.append(tuple(qubits))

    return gates_and_qubits


if __name__ == "__main__":
    # Test
    # Example: Generate QFT for n qubits
    n = 6
    filename = "QASM_files/qft_%squbits.qasm" % n
    # /Users/danielschonberger/Desktop/Heuristic_Github/
    write_qft_to_file(n, filename)
    print(f"QFT for {n} qubits written to {filename}")

    seq = parse_qasm(filename)
    flat_seq = [item for sublist in seq for item in sublist]
    print("sequence: ", seq)

    two_qubit_seq = []
    i = 0
    s = 0
    while i < len(seq):
        if len(seq[i]) == 2:
            two_qubit_seq.append(i + s)
            s += 1
        i += 1

    print("two_qubit_seq: ", two_qubit_seq)
    # print('two_qubit_seq: ', [i for i, gate in enumerate(seq) if len(gate) == 2])
    print("flat sequence: ", flat_seq)
    print(len(flat_seq))
