flex-general
============

This is a FORTRAN 90-based code that applies the Fluctuation Exchange Approximation (FLEX) to study superconductivity and magnetism in multiband tight-binding models.

To compile the code, please type the following in the src-flex/ subdirectory
% make first-order  (this produces a first-order perturbation theory code)
% make second-order (this produces a second-order perturbation theory code)
% make second-order-mpi (same as above, but with mpi extensions for parallel processing)
% make third-order  (this produes a third-order perturbation theory code)
% make third-order-mpi (same as above, but with mpi extensions for parallel processing)
% make flex (this produces a fluctuation exchange approximation code)
% make flex-mpi (same as above, but with mpi extensions for parallel processing)
