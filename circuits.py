def make_initial_Circuit(backend):
    from pyquil import Program
    from pyquil.gates import H, CNOT, Z, MEASURE, S, X, Y, I, RX, RZ, FENCE
    from pyquil.quilbase import Declare
    qubits = [0,1,2,3]
    num_qubits = len(qubits)
    n = 2
    def trotterLayer(h,J,dt,n):
        trotterLayer = Program(Declare("ro", "BIT", num_qubits))
        for q in qubits:
            trotterLayer += RX(dt*4*h, q)
        for i, j in [(qubits[2*i], qubits[2*i+1]) for i in range(n)]:
            trotterLayer += CNOT(i,j)
        for q in [qubits[2*i+1] for i in range(n)]:
            trotterLayer += RZ(-4*J*dt, q)
        for i, j in [(qubits[2*i], qubits[2*i+1]) for i in range(n)]:
            trotterLayer += CNOT(i,j)
        for i, j in [(qubits[2*i+1], qubits[2*i+2]) for i in range(n-1)]:
            trotterLayer += CNOT(i,j)
        for q in [qubits[2*i+1] for i in range(n-1)]:
            trotterLayer += RZ(-4*J*dt, q)
        for i, j in [(qubits[2*i+1], qubits[2*i+2]) for i in range(n-1)]:
            trotterLayer += CNOT(i,j)
        return trotterLayer

    def maketrotterCircuit(s):
        tL = trotterLayer(h=1, J=-0.15, dt=0.2, n=n)
        trotterCircuit = Program(Declare("ro", "BIT", num_qubits))
        for i in range(s):
            trotterCircuit += tL
            #trotterCircuit += FENCE() #the fence is not currently supported by quilc
        return trotterCircuit

    return [maketrotterCircuit(i) for i in range(1,15)]

def make_initial_Circuit2(backend):
    from pyquil import Program
    from pyquil.gates import H, CNOT, Z, MEASURE, S, X, Y, I, RX, RZ, FENCE
    from pyquil.quilbase import Declare
    prog = Program()
    prog += Declare("ro", "BIT", 2)
    prog += CNOT(0,1)
    #prog += CNOT(0,1)
    return [prog]
