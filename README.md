# Autonomous error mitigation based on probabilistic error reduction

Current quantum computers suffer from a level of noise that prohibits extracting useful results directly from longer computations. The figure of merit in many near-term quantum algorithms is an expectation value measured at the end of the computation, which experiences a bias in the presence of hardware noise. A systematic way to remove such bias is probabilistic error cancellation (PEC). PEC requires the noise to be fully characterized and introduces a sampling overhead which increases exponentially with circuit depth, prohibiting high-depth circuits at realistic noise levels. 
Probabilistic error reduction (PER) is a related quantum error mitigation method which largely reducing the sampling overhead by predictably reducing the noise rather than fully cancelling it. In combination with extrapolation techniques to approximate the expectation value in the zero-noise limit, the accuracy of PER can be comparable to that of PEC. PER is broadly applicable to near-term algorithms, and the autonomous implementation of PER is desirable for facilitating its widespread use. This package is intended to be a set of software tools for autonomously applying error mitigation by parsing a user-specified circuit, obtaining the necessary tomography data through a recently-developed method of sparse Pauli tomography, using this data to generate noise-scaled representations of circuits, and analyzing the results to obtain an error-reduced expectation value.

## Installation:
Python 3.10.10 or higher should be installed. In the main folder create a virtual environment with "python -m venv .venv", activate the environment and then run "pip install -r requirements.txt" in it. It should install everything at the correct version.
### Main Environment
* pyquil        4.0.0
* matplotlib:   3.9.2

## Credits
Project Advisor: [Dr. Peter Orth](https://www.uni-saarland.de/lehrstuhl/orth.html)

A special thanks goes out to Prachi Sharma from the work groupe, who took the time to explain many things to me and guide me througth the scientific process.

## References

* [1] A. Caesura, C. L. Cortes, W. Pol, S. Sim, M. Steudtner, G.-L. R. Anselmetti, M. Degroote, N. Moll, R. Santagati, M. Streif, et al., “Faster quantum hemistry simulations on a quantum computer with improved tensor factorization and active volume compilation,” arXiv preprint arXiv:2501.06165, 2025.
* [2] A. Singh, “Prospect of future computing era using quantum computing,” Available at SSRN 5171765, 2025.
* [3] D. Greenbaum, “Introduction to quantum gate set tomography,” arXiv preprint arXiv:1509.02921, 2015.
* [4] E. Knill, D. Leibfried, R. Reichle, J. Britton, R. B. Blakestad, J. D. Jost, C. Langer, R. Ozeri, S. Seidelin, and D. J. Wineland, “Randomized benchmarking of quantum gates,” Physical Review A—Atomic, Molecular, and Optical Physics, vol. 77, no. 1, p. 012307, 2008.
* [5] E. Van den Berg, Z. Minev, A. Kandala, and K. Temme, “Probabilistic error cancellation with sparse pauli-lindblad models on noisy quantum processors (2022),” arXiv preprint arXiv:2201.09866. 
* [6] Z. Minev, “A tutorial on tailoring quantum noise - twirling 101 (parts i–iv),” 2022. Accessed: 2024-12-17 https://www.zlatko-minev.com/blog/twirling.
* [7] A. Mari, N. Shammah, and W. J. Zeng, “Extending quantum probabilistic error cancellation by noise scaling,” Physical Review A, vol. 104, no. 5, p. 052607, 2021.
* [8] K. Temme, S. Bravyi, and J. M. Gambetta, “Error mitigation for short-depth quantum circuits,” Physical review letters, vol. 119, no. 18, p. 180509, 2017.
* [9] T. Giurgica-Tiron, Y. Hindy, R. LaRose, A. Mari, and W. J. Zeng, “Digital zero noise extrapolation for quantum error mitigation,” in 2020 IEEE International Conference on Quantum Computing and Engineering (QCE), pp. 306–316, IEEE, 2020.
* [10] E. F. Dumitrescu, A. J. McCaskey, G. Hagen, G. R. Jansen, T. D. Morris, T. Papenbrock, R. C. Pooser, D. J. Dean, and P. Lougovski, “Cloud quantum computing of an atomic nucleus,” Physical review letters, vol. 120, no. 21, p. 210501, 2018.
* [11] A. He, B. Nachman, W. A. de Jong, and C. W. Bauer, “Zero-noise extrapolation for quantum-gate error mitigation with identity insertions,” Physical Review A, vol. 102, no. 1, p. 012426, 2020.
* [12] V. R. Pascuzzi, A. He, C. W. Bauer, W. A. De Jong, and B. Nachman, “Computationally efficient zero-noise extrapolation for quantum-gate-error mitigation,” Physical Review A, vol. 105, no. 4, p. 042406, 2022.
* [13] Scipy, “scipy.optimize.curve fit.” Accessed: 2025-03-20 https://docs.scipy.org/doc/scipy/reference/ generated/scipy.optimize.curve_fit.html.
* [14] B. McDonough, A. Mari, N. Shammah, N. T. Stemen, M. Wahl, W. J. Zeng, and P. P. Orth, “Automated quantum error mitigation based on probabilistic error reduction,” in 2022 IEEE/ACM Third International Workshop on Quantum Computing Software (QCS), pp. 83–93, IEEE, 2022.
* [15] P. Zhao, K. Linghu, Z. Li, P. Xu, R. Wang, G. Xue, Y. Jin, and H. Yu, “Quantum crosstalk analysis for simultaneous gate operations on superconducting qubits,” PRX quantum, vol. 3, no. 2, p. 020301, 2022.
