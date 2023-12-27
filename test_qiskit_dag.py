from qiskit import QuantumCircuit, converters

# qc = QuantumCircuit.from_qasm_file("/Users/danielschonberger/Desktop/Heuristic_Shuttler/QASM_files/qft_6qubits.qasm")#qft_nativegates_quantinuum_qiskit_opt3_6.qasm")


qc = QuantumCircuit.from_qasm_file(
    "/Users/danielschonberger/Desktop/Heuristic_Shuttler/QASM_files/qft_nativegates_quantinuum_qiskit_opt3_6.qasm"
)
dag_dep = converters.circuit_to_dagdependency(qc)

# # Function to swap commuting gates in a circuit
# def swap_commuting_gates(circuit):
#     modified_circuit = QuantumCircuit(circuit.num_qubits)
#     gates_to_swap = []

#     for i, inst in enumerate(circuit.data):
#         # Example condition for swapping: if Z gate is followed by X gate
#         if inst[0].name == 'x' and i < len(circuit.data) - 1 and circuit.data[i + 1][0].name == 'z':
#             gates_to_swap.append((i, i + 1))

#     skip = False
#     for i, inst in enumerate(circuit.data):
#         if skip:
#             skip = False
#             continue

#         if (i, i + 1) in gates_to_swap:
#             # Swap the gates
#             modified_circuit.append(circuit.data[i + 1][0], circuit.data[i + 1][1])
#             modified_circuit.append(inst[0], inst[1])
#             skip = True
#         else:
#             modified_circuit.append(inst[0], inst[1])

#     return modified_circuit

# # Modify the circuit
# modified_qc = swap_commuting_gates(qc)

# Convert back to DAGDependency
modified_dag_dep = converters.circuit_to_dagdependency(qc)

modified_dag_dep._multi_graph.remove_node(modified_dag_dep.get_node(1).node_id)
modified_dag_dep.draw(filename="/Users/danielschonberger/Desktop/Heuristic_Shuttler/dag_dep_mod3.png")
