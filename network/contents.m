%Network Toolbox
%
%Version 1.0 16-Dec-2000
%
%Directory: network
%
%COVUPUQ  Plot of covariance function cov(u_P,u_Q) for network with
%	       relative directions, see pages 119--120
%
%EX23	    Example 2.3: Computation of eigensolution for R_N(2,-2,-2),
%	       takes the dimension N of the matrix R_N as argument.
%	       The involved functions y are defined in ex23a, ex23b,
%	       ex23c, and ex23d:
%	           y = 1/tanh(N*x)-3*tanh(x); % ex23a
%	           y = tanh(N*x)-3*tanh(x);   % ex23b
%	           y = cot(N*x)+3*tan(x);     % ex23c
%	           y = tan(N*x)-3*tan(x);     % ex23d
%
%EX33	    Example 3.3: Plot of covariance function G(P,Q) for difference
%	       of abscissae between P and Q the distance r apart and with direction
%	       angle phi
%
%FIG16_18 Figures 1.6, 1.7, and 1.8: Discrete approximation of continuous
%	       1-D leveling with n nodes and combinations of Dirichlet and Neumann
%	       boundary conditions
%
%FIG35	 Figure 3.5: Plot of Green's function for a unit circle, theta = pi/2
%
%FIG35A   Green's fuction for a unit circle with Neumann boundary conditions.
%	       This Green's function is the continuous analogue to the covariance
%         matrix for 2-D leveling
%
%FIG36	 Figure 3.6: Plot of Green's function for a unit circle, theta = 0
%
%FIG37	 Figure 3.7: Confocal ellipses
%
%FIG38	 Figure 3.8: Plot of ellipses tending to fill out an infinite parallel
%         strip
%
%FIG44    Figure 4.4: Plot of the Bessel functions J_0 and J_1 of first kind
%
%FIG51	 Figure 5.1: Plot of 1-D density function for the discrete Laplacian
%
%FIG52	 Figure 5.2: Plot of 2-D density function of the discrete Laplacian
%
%FIG53	 Figure 5.3: Plot of 1-D density function of the discrete bi-harmonic
%         operator 
%
%FIG54	 Figure 5.4: Plots of upper and lower bounds for the spectral distribution
%         function N(lambda) for 1-D density function of the discrete Laplacian 
%
%FIG55	 Figure 5.5: Plots of upper and lower bounds for the spectral distribution
%         function N(lambda) for the 2-D density function of the discrete
%         Laplacian
%
%FIG57    Figure 5.7: Plot of Weyl's asymptotic theorem for eigenvalues
%
%FIG65    Figure 6.5: Difference vectors between REFDK and original GI 
%         cordinates in UTM at common fundamental network points.
%	       Courtesy by National Survey and Cadastre, Denmark
%
%GEONET   This script solves the system of PDEs for geodetic networks (4.23)
%	       with homogeneous boundary conditions (4.29). The components of the
%         right side of (4.23) are all taken to be delta - 1/(area of region).
%         The right-hand side of (4.29) is 0. The solution components (u, v)
%         are plotted as vector field arrows. The components alpha and beta
%         are plotted as graphs, also colored according to function value 
%         (figures 2 and 3, respectively).
%
%	       The script uses functions from the finite element software
%	       FEMLAB that is distributed by COMSOL. Written by Daniel Bertilsson
%
%H	       Generation of an H matrix that describes the geometry of a 2-D 
%         leveling network
%
%LEVRECTA Plot of standard deviations (sigmas) of heights in an m by n
%	       regular, rectangular, free leling network
%
%	       Implementation of formulas to be found in
%	           K. Borre & P. Meissl:
%	           Strength Analysis of Leveling-Type Networks,
%	           An Application of Random Walk Theory, Appendix B
%	           Geod\aetisk Institut, Meddelelse No. 50, K\obenhavn 1974
%
%LW	    Computation of the spectral density function yy for the discrete
%         Laplacian operator in 2-D
%
%MYKRON   Computation of normals for an m by n leveling network by means of
%         Kronecker products
%
%NW1	    Simple least-squares problem. Identical to Example 2 on page 185 in
%	           G. Strang: Introduction to Linear Algebra
%
%NW10	    Solution to Problem 10
%
%NW11	    Green's function for a unit disc with Neumann boundary condition.
%	       The delta function is located at the (complex) zeta position. 
%         The function has created the cover figure
%
%NW12ONE  Spectral density function for 1-D discrete Laplacian operator. The
%         mesh number is n (>100)
%
%NW12TWO  Spectral density function for 2-D discrete Laplacian operator. The
%         number of meshes is m (>200) times n (>200)
%
%NW13	    Computes the structure matrix for each triangle in a network.
%         The structure matrix is a generalized weight matrix derived 
%         according to the continous network theory. The network is described
%         by coordinates of the nodes and a list of the triangles. The
%	       coordinates are given in the file 'nodes' and the triangles are
%         given in the file 'triangles'. Typical call: nw13(triangles,nodes)
%
%NW2	    Simple least-squares problem. Identical to Problem 4.3.1 in
%		           G. Strang: Introduction to Linear Algebra
%
%NW3	    For a free leveling line with n nodes we compute the eigenvalues
%         and -vectors of the normal matrix in alternative ways and plot
%         various results.
%
%NW4	    Computation and plot of eigenvectors of the oscillating n by n
%         tri-diagonal matrix (-1,2,-1)
%
%NW5	    Plots a posteriori standard deviation of each leveling point in
%         a regular m by n network. The normals are generated by means of
%         the Kronecker product. The basic problem (k = 3)is the so-called
%         free network with singular N-matrix. If a boundary is kept fixed
%	       we delete the corresponding columns and rows, and dimensions are
%         deminished accordingly.
%
%NW6	    Propagation of systematic errors in 1-D leveling line with n nodes
%         and fixed terminals
%
%NW7	    Regular n by n tridiagonal matrices for which we modify the (1,1) 
%         and (n,n) components
%
%	       Some of the results can also be found in
%	           W.-D. Schuh (1984) Analyse und Konvergenzbeschleunigung der
%	           Methode der konjugierten Gradienten bei geod\"atischen Netzen. 
%             Mitteilungen der geod\"atischen Institute der Technischen 
%             Universit\"at Graz, Folge 49, page 48
%
%NW8	    Discrete approximation of continuous n node 1-D leveling with Neumann
%         boundary conditions
%
%NW9	    Script that performs all necessary steps in transformation of 
%         observations and computation of covariance matrices, etc. for a
%         single triangle
%
%S2andS4  Computation of stiffness matrix S using both type 2 and type 4 
%         observations for the triangle with vertices (x1,y1),(x2,y2), and (x3,y3)
%	       A typical call:  S = s2ands4(x1, y1, x2, y2, x3, y3, weight2, weight4)
%
%STIFCHEC Check of stiffness matrix S for observation types 2 and 4 for adjacent
%         isoscele triangles. First set with SW-NE diagonal and next with NW-SE
%         diagonal
%
%STIFFNES Computation of stiffness matrix S for the triangle with vertices (x1,y1),
%         (x2,y2), and (x3,y3). A typical call S = stiffnes(x1, y1, x2, y2, x3, y3)
%
%%%%%%%%%%%%%%%%%%%%%%%% end contents.m %%%%%%%%%%%%%%%%%%
