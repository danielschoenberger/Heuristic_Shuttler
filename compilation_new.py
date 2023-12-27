from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dagdependency


def get_dag_dependency(filename):
    """Get the DAGDependency from a QASM file."""
    qc = QuantumCircuit.from_qasm_file(filename)
    return circuit_to_dagdependency(qc)


def get_front_layer(dag):
    """Get the front layer of the DAG."""
    front_layer = []
    for node in dag.get_nodes():
        # If a node has no predecessors, it's in the front layer
        if not dag.direct_predecessors(node.node_id):
            front_layer.append(node)
    return front_layer


def remove_node(dag, node):
    """Execute a node and update the DAG (remove the node and its edges)."""
    qubits = [qubit.index for qubit in node.qargs]
    # clbits = [clbit.index for clbit in node.cargs]
    print(f"Remove gate: {node.name}, Qubits: {qubits}")  # , Classical bits: {clbits}")

    dag._multi_graph.remove_node(node.node_id)


def find_best_gate(front_layer):
    """Find the best gate to execute based on a cost function."""
    # Example: simple cost function (just choose the first gate)
    # Implement your own cost function here
    return front_layer[0] if front_layer else None


if __name__ == "__main__":
    n = 6
    filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/QASM_files/qft_nativegates_quantinuum_qiskit_opt3_6.qasm"  # "QASM_files/qft_%squbits.qasm" % n
    dag_dep = get_dag_dependency(filename)
    while True:
        front_layer = get_front_layer(dag_dep)
        f = []
        for node in front_layer:
            print(node.qargs)
            f.append(tuple(qubit.index for qubit in node.qargs))
        print(f)
        if not front_layer:
            break  # Exit if no more gates to execute

        best_gate = find_best_gate(front_layer)
        if best_gate:
            remove_node(dag_dep, best_gate)
            # Here, you would actually execute the gate in your quantum circuit
            print(f"Executed gate: {best_gate.name}")
