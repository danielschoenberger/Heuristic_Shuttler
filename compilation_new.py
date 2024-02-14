import math

from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dagdependency
from qiskit.dagcircuit import DAGDependency

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
    # if dag.direct_successors(node.node_id):
    #    for successor in dag.direct_successors(node.node_id):
    #        dag._multi_graph.remove_edge(node.node_id, successor)
    dag._multi_graph.remove_node(node.node_id)


def find_best_gate(front_layer, dist_map):
    """Find the best gate to execute based on distance."""
    min_gate_cost = math.inf
    for _i, gate_node in enumerate(front_layer):
        qubit_indices = gate_node.qindices
        gate_cost = max([dist_map[qs] for qs in qubit_indices])
        # if both ions of 2-qubit gate are in pz execute 2-qubit gate
        if len(qubit_indices) == 2 and gate_cost == 0:
            return gate_node
        if gate_cost < min_gate_cost:
            min_gate_cost = gate_cost
            best_gate = gate_node
    return best_gate


def manual_copy_dag(dag):
    new_dag = DAGDependency()

    # Recreate quantum registers in the new DAG
    for qreg in dag.qregs.values():
        new_dag.add_qreg(qreg)

    # Iterate over all operation nodes in the original DAG and copy them
    for node in dag.get_nodes():
        new_dag.add_op_node(node.op, node.qargs, node.cargs)

    return new_dag


def update_sequence(dag, dist_map):
    """Get the sequence of gates from the DAG. Creates a new DAG and removes all gates from it while creating the sequence."""
    working_dag = manual_copy_dag(dag)
    sequence = []
    i = 0
    while True:
        first_gates = get_front_layer(working_dag)
        if not first_gates:
            break
        first_gate_to_excute = find_best_gate(first_gates, dist_map)
        if i == 0:
            first_node = first_gate_to_excute
        i = 1
        remove_node(working_dag, first_gate_to_excute)
        sequence.append(first_gate_to_excute.qindices)
    return sequence, first_node


if __name__ == "__main__":
    # n = 6
    # filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/QASM_files/qft_nativegates_quantinuum_qiskit_opt3_6.qasm"  # "QASM_files/qft_%squbits.qasm" % n
    n = 3
    filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/QASM_files/qft_nativegates_quantinuum_tket_3_short.qasm"  # "QASM_files/qft_%squbits.qasm" % n

    qc = QuantumCircuit.from_qasm_file(filename)

    dag_dep = circuit_to_dagdependency(qc)
    filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/dags_test/dag_dep.pdf"
    dag_dep.draw(filename=filename)

    # dag_test = circuit_to_dag(qc)
    # filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/dags_test/dag_dep_%s.png" % 16
    # dag_test.draw(filename=filename)

    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)

    # remove_node(dag_dep, next_node)
    # sequence, next_node = update_sequence(dag_dep, {0: 2, 1: 1, 2: 3})
    # print(sequence)
    # dag_dep = manual_copy_dag(dag_dep)


# if __name__ == "__main__":
#     n = 6
#     filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/QASM_files/qft_nativegates_quantinuum_qiskit_opt3_6.qasm"  # "QASM_files/qft_%squbits.qasm" % n
#     qc = QuantumCircuit.from_qasm_file(filename)
#     print(qc.to_gate().definition)
#     dag_dep = circuit_to_dagdependency(qc)
#     print([node.node_id for node in dag_dep.get_nodes()])
#     filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/dags_test/dag_dep_%s.png" % 10
#     dag_dep.draw(filename=filename)

#     sequence, next_node = update_sequence(dag_dep, {3: 2, 4: 1, 5: 3})
#     print(sequence)
#     remove_node(dag_dep, next_node)
#     print([node.node_id for node in dag_dep.get_nodes()])
#     filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/dags_test/dag_dep_%s.png" % 11
#     dag_dep.draw(filename=filename)

#     sequence, next_node = update_sequence(dag_dep, {3: 2, 4: 1, 5: 3})
#     print(sequence)

#     remove_node(dag_dep, next_node)
#     sequence, next_node = update_sequence(dag_dep, {3: 2, 4: 1, 5: 3})
#     print(sequence)

#     remove_node(dag_dep, next_node)
#     sequence, next_node = update_sequence(dag_dep, {3: 2, 4: 1, 5: 3})
#     print(sequence)

#     remove_node(dag_dep, next_node)
#     sequence, next_node = update_sequence(dag_dep, {3: 2, 4: 1, 5: 3})
#     print(sequence)


# i = 0
# while True:
#     # plot the dag
#     filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/dags_test/dag_dep_%s.png" % i
#     dag_dep.draw(filename=filename)

#     # create qubit map to access the qubit indices of the DAGDependency
#     front_layer = get_front_layer(dag_dep)
#     for node in front_layer:
#         qbit_indices = node.qindices
#         print("qbit_indices of front layer gates: ", qbit_indices)
#     if not front_layer:
#         break  # Exit if no more gates to execute
#     distance_map = {3: 2 + i, 4: 1 + i, 5: 3 - i}
#     best_gate = find_best_gate(front_layer, distance_map)
#     remove_node(dag_dep, best_gate)
#     i += 1
#     print(f"Executed gate: {best_gate.name}\n")
