from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, converters

q = QuantumRegister(2, "q")
c = ClassicalRegister(2, "c")
qc = QuantumCircuit(q, c)

qc.h(q[0])
qc.cx(q[0], q[1])
qc.measure(q, c)

# dag = dagcircuit.DAGCircuit(qc)
dag = converters.circuit_to_dag(qc)
dag_dep = converters.circuit_to_dagdependency(qc)

dag.draw(filename="/Users/danielschonberger/Desktop/Heuristic_Shuttler/dag.png")
dag_dep.draw(filename="/Users/danielschonberger/Desktop/Heuristic_Shuttler/dag_dep.png")
