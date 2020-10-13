%% Demo

choice=-1;
 while choice~=0
    fprintf("\n    Select a Matrix \n");
    fprintf("\n    === DENSE ===\n");
    fprintf("\n    1. airy_demo            - An Airy operator. \n");
    fprintf("    2. basor_demo           - Toeplitz matrix (Basor-Morrison). \n");
    fprintf("    3. chebspec_demo        - First order Chebyshev differentiation matrix. \n");
    fprintf("    4. companion_demo       - A companion matrix. \n");
    fprintf("    5. convdiff_demo        - 1-D convection diffusion operator. \n");
    fprintf("    6. davies_demo          - Davies' example. \n");
    fprintf("    7. demmel_demo          - Demmel's matrix. \n");
    fprintf("    8. frank_demo           - The Frank matrix. \n");
    fprintf("    9. gallery3_demo        - The matrix `gallery(3)'. \n");
    fprintf("    10. gallery5_demo        - The matrix `gallery(5)'. \n");
    fprintf("    11. gaussseidel_demo     - Gauss-Seidel iteration matrices. \n");
    fprintf("    12. godunov_demo         - Godunov's matrix. \n");
    fprintf("    13. grcar_demo           - Grcar's matrix. \n");
    fprintf("    14. hatano_demo          - Hatano-Nelson example. \n");
    fprintf("    15. kahan_demo           - The Kahan matrix. \n");
    fprintf("    16. landau_demo          - An application from lasers. \n");
    fprintf("    17. orrsommerfeld_demo   - An Orr-Sommerfeld operator. \n");
    fprintf("    18. random_demo          - A matrix with random entries. \n");
    fprintf("    19. randomtri_demo       - An upper triangular matrix with random entries. \n");
    fprintf("    20. riffle_demo          - The riffle shuffle matrix. \n");
    fprintf("    21. transient_demo       - A matrix with transient behaviour. \n");
    fprintf("\n    === SPARSE ===\n");
    fprintf("      23. convdiff_fd_demo     - A convection diffusion operator (finite-differences). \n");
    fprintf("      24. markov_demo          - A Markov chain transition matrix. \n");
    fprintf("      25. sparserandom_demo    - A sparse random matrix. \n");
    fprintf("      26. skewlap3d_demo       - A non-symmetric 3D Laplacian. \n");
    fprintf("      27. supg_demo            - An SUPG matrix. \n");
    fprintf("\n    === CUSTOM ===\n");
    fprintf("      27. five_lap             - 5-point Laplacian. \n");
    fprintf("      28. ichol_lap            - Incomplete cholseky of 5-point Laplacian. \n");
    fprintf("      29. trefethen            - Toeplitz matrix by Trefethen. \n");
    fprintf("      30. bidiag               - Lower diagonal of -1, diagonal of 1. \n");
    fprintf('      0. exit \n\n');

    fprintf('              select a menu option> ')
    choice = floor(str2num(input('','s'))); 
    while (length(choice) == 0) | (choice<0) | (choice>30)
       fprintf('              select a menu option> ')
       choice = floor(str2num(input('','s'))); 
    end

	n = -1;
    fprintf('              select matrix size from [4,100]> ')
    n = floor(str2num(input('','s'))); 
    while (length(choice) == 0) | (n<4) | (n>100)
       fprintf('              select a matrix size > ')
       n = floor(str2num(input('','s'))); 
    end

    if isnumeric(choice)
       switch choice 
		  case 1, A = airy_demo(n);
		  case 2, A = basor_demo(n);
		  case 3, A = chebspec_demo(n);
		  case 4, A = companion_demo(n);
		  case 5, A = convdiff_demo(n);
		  case 6, A = davies_demo(n);
		  case 7, A = demmel_demo(n);
		  case 8, A = frank_demo(n);
		  case 9, A = gallery3_demo(n);
		  case 10, A = gallery5_demo(n);
		  case 11, A = gaussseidel_demo(n);
		  case 12, A = godunov_demo(n);
		  case 13, A = grcar_demo(n);
		  case 14, A = hatano_demo(n);
		  case 15, A = kahan_demo(n);
		  case 16, A = landau_demo(n);
		  case 17, A = orrsommerfeld_demo(n);
		  case 18, A = random_demo(n);
		  case 19, A = randomtri_demo(n);
		  case 20, A = riffle_demo(n);
		  case 21, A = transient_demo(n);
		  case 22, A = dwave_demo(n);
		  case 23, A = convdiff_fd_demo(n);
		  case 24, A = markov_demo(n);
		  case 25, A = sparserandom_demo(n);
		  case 26, A = skewlap3d_demo(n);
		  case 27, A = supg_demo(n);
		  case 28, A = lap(n);
		  case 29, A = ichol(n);
		  case 30, A = trefethen(n);
		  case 31, A = bidiag(n);
       end
    end

	if(choice > 0) 
		A = full(A); 
	end
	if(choice < 29) 
		A = tril(A); 
	end
	[m,n] = size(A);

	fprintf('\n   We start with the lower triangular portion of the convection diffusion matrix of size %d x %d \n', m,n);
		('\n   First, we plot the diagonal values (in black) and computed eiganluves (in red) via eig() \n');
	input('                - - -   press <return> to continue   - - - \n\n');
		
	% Note that eigs is for sparse matrices
	myfig1=figure('visible','off'); 
	figure(myfig1), clf
		
	ew = eigs(A,round(n/2));
	diag_A = diag(A);
	plot( diag_A , zeros(n), 'kx' ), hold on, axis equal;
	plot( real(ew), imag(ew), 'ro'), hold on;
		
	% ew = eigs(L',n/2-2);
	% plot( real(ew), imag(ew), 'gd'); 
		
	fprintf('\n   Next we will plot the pseudospectra with various epsilon. Notice that the eigenvalues extend up until 10^-10 \n');
	input('                - - -   press <return> to continue   - - - \n\n');
		
	% http://www.cs.ox.ac.uk/pseudospectra/eigtool/documentation/options/index.html
	% opts.npts determines the number of points to sample to compute psuedospectra
	% opts.k determiens the number of eigenvalues to search for
    A = full(A);
	clear opts
	% opts.ax = [diag_L(1) 100 -100 100];
	 opts.print_plot = 1;
	 opts.no_ews = 0;
	 opts.levels = [-16:-12];
	 opts.npts   = 50;
	 opts.k = 10;
	 opts.colourbar = 1;
	 eigtool(A, opts, myfig1);
end
