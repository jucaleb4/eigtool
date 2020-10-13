choice=-1;
 while choice~=0
	print_all_options();

    fprintf('              select a menu option> ')
    choice = floor(str2num(input('','s'))); 
    while (length(choice) == 0) | (choice<0) | (choice>32)
       fprintf('              select a menu option> ')
       choice = floor(str2num(input('','s'))); 
    end

	n = -1;
    fprintf('              select matrix size from [4,...)> ')
    n = floor(str2num(input('','s'))); 
    while (length(choice) == 0) | (n<4) 
       fprintf('              select a matrix size > ')
       n = floor(str2num(input('','s'))); 
    end

    name = "";
    if isnumeric(choice)
		[A,name] = select_matrix(choice, n);
    end

	if(choice > 31) 
		% A = full(A);
        % A = tril(A);
        try
            A = ichol(A);
        catch
            % ref. https://stackoverflow.com/questions/12895228/ichol-as-cholinc-replacement-nonpositive-pivot
            droptol=1e-7;
            [A,~] = ilu(A,struct('type','ilutp','droptol',droptol,'thresh',0));
        end
	end
	if(choice < 29) 
		A = tril(A); 
	end
	[m,n] = size(A);

	%% Experimental Setup
	num_iter = 2; % numiter for gradient descent
	global glb_ct;  % keep track of plots
	glb_ct = 1;
	tol = 1e-6;

	%% Global Parameters for initialization
	% gd_diag_prec 
	display_spy = false;
	if(display_spy)
		figure(glb_ct);
		spy(A)
		glb_ct = glb_ct + 1;
	end

	% backtrack_ls
	global print_bt_scores;
	print_bt_scores = false;
    
	% preconditioner computation
	global plot_gd_score % plots objective value
	plot_gd_score = true;
    compute_diag = true;
    save_diag = true;
    
    if save_diag,
        fname = extractAfter(name, strfind(name, "/"));
    end
    
    if compute_diag     
        fprintf("Computing diagonal preconditioner with %d iterations of gradient descent... \n", num_iter);
        d = gd_diag_prec(A,num_iter);
        if save_diag,
            sname = sprintf("mats/%s_diag%d", fname, num_iter);
            save(sname, 'd');
        end
    else
        fname = extractAfter(name, strfind(name, "/"));
        sd = load(sprintf("mats/%s_diag%d", fname, num_iter));
        d = sd.d;
    end
    
    fprintf("Finished computing diagonal preconditioner. \n");
	
	%% Experiments
	% jacobi_niters
	global plot_jacobi_res
	plot_jacobi_res = true;
	jacobi_analysis(A, d, tol);
	
    % gmres_analysis
	global plot_gmres_res
	plot_gmres_res = true;
	gmres_analysis(A, d);
    
    % pcg_analysis
    global plot_pcg_res
    % plot_pcg_res = true
    % pcg_analysis(A, d);
	
    % preconditioner_analysis
    global eig_analysis;
    eig_analysis = true;
	global plot_diags
	plot_diags = true;
	global print_diags
	print_diags = false;
	global disp_precond_mat
	disp_precond_mat = false;
	global print_precond_off_diag
	print_precond_off_diag = false;
	% preconditioner_analyis(A, d, "argname");

	% FD verification
	% fd_gradient(A);
    input('              === press ENTER to continue === ')
end
	
%% Preconditioner Analysis
function x = eigvec_cond(A)
	[m,n] = size(A);
	[V,D] = eigs(A,n);
	x = cond(V);
end

function preconditioner_analyis(A, d, mat_name)
	% Given a matrix @A named @mat_name, computes condition of its eigenvecs
	% and then after a diagonal preconditioning (derived from the function diag_prec)
    global eig_analysis;
    if eig_analysis,
        condV = eigvec_cond(A);
        condA = cond(A,2);
        fprintf("%s eigenvector cond=%.4e. matrix cond=%.4e\n",mat_name, condV, condA);
    end

	% d = gd_diag_prec(A,num_iter);
	D = sparse(diag(d));
	Dinv = sparse(diag(d.^-1));
	precond_A = D*A*Dinv;
	condV = eigvec_cond(precond_A);
	condA = cond(precond_A,2);
	fprintf("%s preconditioned. eigenvector cond=%.4e. matrix cond=%.4e\n",mat_name, condV, condA);
	global print_diags
	if print_diags
		fprintf("diag=[");
		fprintf("%.2e, ", d);
		fprintf("]\n\n");
	end
	global plot_diags
	if plot_diags
		global glb_ct;
		figure(glb_ct);
		glb_ct = glb_ct+1;
		plot(d)
		title("Diagonals of Preconditioner (Top/Left to Bottom/Right");
	end
	global disp_precond_mat
	if disp_precond_mat
		disp(precond_A);
	end
	global print_precond_off_diag
	if print_precond_off_diag
		fprintf("Upr diag = [");
		fprintf("%.2e, ", diag(precond_A,1));
		fprintf("]\n");

		fprintf("Sub diag = [");
		fprintf("%.2e, ", diag(precond_A,-1));
		fprintf("]\n");
	end
end

%% Preconditioner Calculation
function d = gd_diag_prec(A,num_iter)
	% Computes a diagonal precondition via gradient descent
	% for a matrix @A using @num_iter iterations. 
	% Returns: array @d of the diagonal elements

	% parameters for backtrack
	use_backtrack = true;
	start = 2^10;
	c = 0.5;
	tau = 0.5;

	[m,n]=size(A);
	AA = A.*A;
	trace_A2 = trace(A*A);
	d=ones(n,1);
    % good starting guess
    % d=2.^-(0:n-1);
    % d=d';
    
	first_norm = obj_f(A,d);

	lst=[first_norm];
	pr_list=[ norm(A,'fro') ];
	dff=[];
	alpha=.1;

	for i=1:num_iter
		u=(d).*(AA*(d.^-2)) - (d.^-3).*(AA'*(d.^2));
		v=( (d.^2)'*(AA)*(d.^-2) )^2;
		grad=(4*trace_A2/v)*u;
		if use_backtrack
			alpha = backtrack_ls(A,d,grad,start,c,tau);
		end
		d = d - alpha*grad;
		last_norm = obj_f(A,d);

		lst = [lst, last_norm];
    end
    
	% plot norm of the objective
	global plot_gd_score
	if plot_gd_score
		global glb_ct;
		figure(glb_ct);
		glb_ct = glb_ct+1;
		plot( (0:num_iter), lst );
		title("Gradient Descent Objective Score Over Iterations");
	end
end

function val = obj_f(A,d)
	B = sparse(diag(d))*(A)*sparse(diag(d.^-1));
	val = ( norm(B-B','fro')/norm(B,'fro') )^2;
end

function alpha = backtrack_ls(A,d,grad,start,c,tau)
	% Beings with @start step size and repeatedly decreases step by @tau 
	% until the Armijo-Goldstein conditions are satisifed
	if 0>=c || c>=1
		fprintf("ERROR: Invalid @c. Setting to 1/2\n");
		c = 0.5;
	end
	if 0>=tau || tau>=1
		fprintf("ERROR: Invalid @tau. Setting to 1/2\n");
		tau = 0.5;
	end

	m = dot(grad, grad);
	f_0 = obj_f(A,d);
		
	alpha = start;
	curr_score = obj_f(A,d-alpha*grad);
	scores = [f_0, curr_score];
	% Armijo-Goldstein condition
	diff = [curr_score - f_0 + start*c*m];
	% while curr_score - f_0 > -alpha*c*m
	while ~( f_0 - curr_score >= alpha*c*m )
		alpha = tau*alpha;
		curr_score = obj_f(A,d-alpha*grad);
		scores = [scores, curr_score];
		% diff = [diff, curr_score - f_0 + alpha*c*m];
		diff = [diff, curr_score - f_0];
    end

	global print_bt_scores;
	if print_bt_scores
		fprintf("Original: %.4e\n", obj_f(A,d));
		fprintf("[");
		fprintf("%.4e, ", scores );
		fprintf("]\n\n");
		fprintf("Diff=[");
		fprintf("%.4e, ", diff );
		fprintf("]\n\n");
	end
end

%% Correctness of Gradient
function fd_gradient(A)
	% Uses finite different approximation of the grandient to 
	% verify correctness of the gradient

	[m,n]=size(A);
	AA = A.*A;
	trace_A2 = trace(A*A);

	% have random starting point
	d=rand(n,1);
    
	u=(d).*(AA*(d.^-2)) - (d.^-3).*(AA'*(d.^2));
	v=( (d.^2)'*(AA)*(d.^-2) )^2;
	% gradient at d
	grad=(4*trace_A2/v)*u;

	eps=0.0001;
	% e=@(k,n) [zeros(k-1,1);1;zeros(n-k,1)]
	abs_err = [];
	rel_err = [];

	for i=(1:n)
		% e_i = e(i,n);
		e_i = [zeros(i-1,1);1;zeros(n-i,1)];
		d_fwd = d+(eps*e_i);
		d_bwd = d-(eps*e_i);
		grad_guess = dot(e_i, grad);
		grad_fd_approx = ( obj_f(A,d_fwd) - obj_f(A,d_bwd) )/(2*eps);
		abs_err = [abs_err, abs(grad_fd_approx-grad_guess)];
		rel_err = [rel_err, abs(grad_fd_approx-grad_guess)/grad_guess];
	end

	global glb_ct;
	figure(glb_ct);
	glb_ct = glb_ct+1;
	plot( (1:n), abs_err );
	title("Finite-Approximation Gradient Error Along e_i (absolute err)");
    figure(glb_ct);
	plot( (1:n), rel_err );
	title("Finite-Approximation Gradient Error Along e_i (relative err)");
	glb_ct = glb_ct+1;
end

%% Iterative Methods
function gmres_analysis(A, d)
	[m,n] = size(A);
	b = ones(n,1);
    b = randn(n,1);
    b_norm = norm(b);
    
	restart = 10;
	restol = 1e-6; % default is 1e-6
    maxit=1000;    % to prevent premature termination
	[~,flag,res_wo_pre,~,res_vec_wo_pre] = gmres(A,b,restart,restol,maxit);
	fprintf("GMRES not convg:%d w/o preconditioner w/ rel. res=%.2e in %.d iters\n", ...
        flag, res_wo_pre, numel(res_vec_wo_pre) );
	if( res_vec_wo_pre(numel(res_vec_wo_pre)) == 0. )
		res_vec_wo_pre(numel(res_vec_wo_pre)) = 1e-20;
	end

	% d = gd_diag_prec(A,prec_num_iter);
	D = sparse(diag(d));
	Dinv = sparse(diag(d.^-1));
	precond_A = D*A*Dinv;
	precond_b = D*b;

	[~,flag,res_w_pre,~,res_vec_w_pre] = gmres(precond_A,precond_b,restart,restol,maxit);
	fprintf("GMRES not convg:%d w/ preconditioner w/ rel. res=%.2e in %.d iters\n", ...
        flag, res_w_pre, numel(res_vec_w_pre) );
	if( res_vec_w_pre(numel(res_vec_w_pre)) == 0. )
		res_vec_w_pre(numel(res_vec_w_pre)) = 1e-20;
    end
    
    % convert residual norms to relative
    res_vec_wo_pre = res_vec_wo_pre / b_norm;
    res_vec_w_pre = res_vec_w_pre / b_norm;

	global plot_gmres_res
	if plot_gmres_res
		global glb_ct;
		figure(glb_ct);
		glb_ct = glb_ct+1;
		semilogy( (0:numel(res_vec_wo_pre)-1 ), res_vec_wo_pre, '-r' ); hold on;
		semilogy( (0:numel(res_vec_w_pre)-1 ), res_vec_w_pre, '-b' ); hold off;
		legend({'w/o Preconditioner', 'w/ Preconditioner'}, 'Location', 'southwest');
		title_name = "2-norm of relative residual error of GMRES";
		title(title_name);
		ylabel("Relative Residual Error");
		xlabel("Iteration");
		% ylim([1e-17,1]);
	end
end

function jacobi_analysis(A, d, tol)
	res_vec_wo_pre = jacobi_niters(A, tol, false);
	fprintf("Non-preconditioned Jacobi with tol=%.1e in %.d iters\n", tol, numel(res_vec_wo_pre)-1);

	% d = gd_diag_prec(A,prec_num_iter);
	D = sparse(diag(d));
	Dinv = sparse(diag(d.^-1));
	precond_A = D*A*Dinv;
	res_vec_w_pre = jacobi_niters(precond_A, tol, true);

	fprintf("Preconditioned Jacobi with tol=%.1e in %.d iters\n", tol, numel(res_vec_w_pre)-1);

	global plot_jacobi_res
	if plot_jacobi_res
		global glb_ct;
		figure(glb_ct);
		glb_ct = glb_ct+1;
		semilogy( (0:numel(res_vec_wo_pre)-1 ), res_vec_wo_pre, '-r' ); hold on;
		semilogy( (0:numel(res_vec_w_pre)-1 ), res_vec_w_pre, '-b' ); hold off;
		title_name = "2-norm of relative residual error in Jacobi iteration ";
		legend({'w/o Preconditioner', 'w Preconditioner'}, 'Location', 'southwest');
		ylabel("Relative Residual Error");
		xlabel("Iteration");
		title(title_name);
	end
end

function res_vec = jacobi_niters(A,tol,use_prec)
	[m,n] = size(A);
	b = ones(n,1);
    b = randn(n,1);
	% x_sol = linsolve(full(A),b);
	% x_sol_norm = norm(b-A*x_sol, 2);
	b_norm = norm(b,2);

	d = diag(A);
	Dinv = diag(d.^-1);
	E = -tril(A,-1);
	F = -triu(A,1);
	f = Dinv*b;
	G = Dinv*(E+F);

	x = rand(n,1);
	r = b - A*x;
	rel_res_norm = norm(r,2)/b_norm;
	res_vec = [rel_res_norm];
	while rel_res_norm > tol
		x = G*x + f;
		r = b - A*x;
		rel_res_norm = norm(r,2)/b_norm;
		res_vec = [res_vec, rel_res_norm];
	end
	if( res_vec(numel(res_vec)) == 0. )
		res_vec(numel(res_vec)) = 1e-20;
	end
end

function pcg_analysis(A, d)
    [m,n] = size(A);
    b = randn(n,1);
    b_norm = norm(b);
    
	restol = 1e-6; % default is 1e-6
    maxit=1000;    % to prevent premature termination
	[~,flag,res_wo_pre,~,res_vec_wo_pre] = pcg(A,b,restol,maxit);
	fprintf("PCG not convg:%d w/o preconditioner w/ rel. res=%.2e in %.d iters\n", ...
        flag, res_wo_pre, numel(res_vec_wo_pre) );
	if( res_vec_wo_pre(numel(res_vec_wo_pre)) == 0. )
		res_vec_wo_pre(numel(res_vec_wo_pre)) = 1e-20;
	end

	% d = gd_diag_prec(A,prec_num_iter);
	D = sparse(diag(d));
	Dinv = sparse(diag(d.^-1));
	precond_A = D*A*Dinv;
	precond_b = D*b;

	[~,flag,res_w_pre,~,res_vec_w_pre] = pcg(precond_A,precond_b,restol,maxit);
	fprintf("PCG not convg:%d w/ preconditioner w/ rel. res=%.2e in %.d iters\n", ...
        flag, res_w_pre, numel(res_vec_w_pre) );
	if( res_vec_w_pre(numel(res_vec_w_pre)) == 0. )
		res_vec_w_pre(numel(res_vec_w_pre)) = 1e-20;
    end
    
    % convert residual norms to relative
    res_vec_wo_pre = res_vec_wo_pre / b_norm;
    res_vec_w_pre = res_vec_w_pre / b_norm;

	global plot_pcg_res
	if plot_pcg_res
		global glb_ct;
		figure(glb_ct);
		glb_ct = glb_ct+1;
		semilogy( (0:numel(res_vec_wo_pre)-1 ), res_vec_wo_pre, '-r' ); hold on;
		semilogy( (0:numel(res_vec_w_pre)-1 ), res_vec_w_pre, '-b' ); hold off;
		legend({'w/o Preconditioner', 'w/ Preconditioner'}, 'Location', 'southwest');
		title_name = "2-norm of relative residual error of PCG";
		title(title_name);
		ylabel("Relative Residual Error");
		xlabel("Iteration");
    end
end

%% Matrix Gen
function A = lap(n)
	% Returns discretized 2D Laplacian equation of dimension (n-2)^2
	R = 'S';
	G = numgrid(R,n);
	A = delsq(G);
end

function T = trefethen(n)
	% Returns tridiagonal mat of [1,0,1/4] of dimension n
	T = diag(0*ones(1,n)) + diag(1*ones(1,n-1),1) + diag(1/4*ones(1,n-1),-1);
end

function A = bidiag(n)
	% Returns bidiagonal (lower triangular) mat of [-1,1] of dimension n
	A = diag(1*ones(1,n)) - diag(ones(1,n-1),-1);
end

%% MISC
function print_all_options()
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
    fprintf("      28. five_lap             - 5-point Laplacian. \n");
    fprintf("      29. ichol_lap            - Incomplete cholseky of 5-point Laplacian. \n");
    fprintf("      30. trefethen            - Toeplitz matrix by Trefethen. \n");
    fprintf("      31. bidiag               - Lower diagonal of -1, diagonal of 1. \n");
    fprintf("      32. input                - input matrix \n");
    fprintf('      0. exit \n\n');
end

function [A,name] = select_matrix(choice, n)
    name = "";
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
		case 29, A = ichol(lap(n));
		case 30, A = trefethen(n);
		case 31, A = bidiag(n);
        case 32,  fprintf('              input filename> ')
                  fname = input('','s'); 
                  name = fname;
                 f = load(fname);
                 A = f.Problem.A;
	end
end
