import matplotlib.pyplot as plt
import pydot
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dagdependency

filename = "/Users/danielschonberger/Desktop/Heuristic_Shuttler/QASM_files/qft_nativegates_quantinuum_tket_3_short.qasm"  # "QASM_files/qft_%squbits.qasm" % n

qc = QuantumCircuit.from_qasm_file(filename)

# Convert the quantum circuit to a DAG
dag = circuit_to_dagdependency(qc)

nodes = dag.get_nodes()
edges = dag.get_all_edges()


# Create a new graph
graph = pydot.Dot(graph_type="digraph", rankdir="LR")  # 'LR' for horizontal, 'TB' for vertical

# Add nodes to the graph
for node in nodes:
    label = f"{node.node_id}: {node.op.name}" if hasattr(node, "op") else str(node.node_id)
    graph.add_node(pydot.Node(str(node.node_id), label=label))

# Process edges: Add edges to the graph
for edge in edges:
    src = str(edge[0])
    dest = str(edge[1])
    graph.add_edge(pydot.Edge(src, dest))

# Save the graph to a file
graph.write_pdf("dag_visualization_corrected_short.pdf")


plt.show()
