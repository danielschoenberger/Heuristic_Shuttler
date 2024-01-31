import math

from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dagdependency

# def get_dag_dependency(filename):
#     """Get the DAGDependency from a QASM file."""
#     qc = QuantumCircuit.from_qasm_file(filename)
#     return circuit_to_dagdependency(qc)


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
    qubits_ind = node.qindices
    print(f"Remove gate: {node.name}, Qubits: {qubits_ind}")  # , Classical bits: {clbits}")

    dag._multi_graph.remove_node(node.node_id)


def find_best_gate(front_layer, dist_map):
    """Find the best gate to execute based on distance."""
    min_gate_cost = math.inf
    for _i, gate_node in enumerate(front_layer):
        qubit_indices = gate_node.qindices
        gate_cost = sum([dist_map[qs] for qs in qubit_indices])
        # if both ions of 2-qubit gate are in pz execute 2-qubit gate
        if len(qubit_indices) == 2 and gate_cost == 0:
            return gate_node
        if gate_cost < min_gate_cost:
            min_gate_cost = gate_cost
            best_gate = gate_node
    return best_gate

    # return front_layer[0] if front_layer else None


if __name__ == "__main__":
    n = 6
    filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/QASM_files/qft_nativegates_quantinuum_qiskit_opt3_6.qasm"  # "QASM_files/qft_%squbits.qasm" % n
    qc = QuantumCircuit.from_qasm_file(filename)
    print(qc.to_gate().definition)
    dag_dep = circuit_to_dagdependency(qc)
    i = 0
    while True:
        # plot the dag
        filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/dags_test/dag_dep_%s.png" % i
        dag_dep.draw(filename=filename)

        # create qubit map to access the qubit indices of the DAGDependency
        front_layer = get_front_layer(dag_dep)
        for node in front_layer:
            qbit_indices = node.qindices
            print("qbit_indices of front layer gates: ", qbit_indices)
        if not front_layer:
            break  # Exit if no more gates to execute
        distance_map = {3: 2 + i, 4: 1 + i, 5: 3 - i}
        best_gate = find_best_gate(front_layer, distance_map)
        remove_node(dag_dep, best_gate)
        i += 1
        print(f"Executed gate: {best_gate.name}\n")
