function [Psi ,E_exact] = exact_diagonalization(H)
% Perform exact diagonalization on the Hamiltonian matrix
% Input:
%   H: Hamiltonian matrix
% Output:
%   E: Eigenvalues of the Hamiltonian matrix
%   psi: Eigenvectors of the Hamiltonian matrix

[Psi, E_exact] = eig(H);
% E = diag(E);
% Find the lowest energy (minimum eigenvalue)
E_exact = min(diag(E_exact));
end