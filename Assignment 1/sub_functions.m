%% sub-functions:
classdef sub_functions
    methods(Static)
        
        function [A]= generate_A(h)
            r = 1 ;
            A = zeros(18,18);
            for m = 1:18
                for n = 1:18
                    a_mn =  sqrt((h + r*sin(((m*pi)/18)) - r*sin(((n*pi)/18))).^2 + (r*cos((m*pi)/18) - r*cos((n*pi)/18)).^2);
                    A(m,n) = 1 ./ (4*pi*a_mn) ;
                end
            end
        end %end - generate_A

        function A = generate_A_q2d(h)
            r = 1;
            A = zeros(18,18);
            for m = 1:18
                for n = 1:18
                    a_mn = (h+r*sin((m*pi)/18) - r*sin((n*pi)/18)).^2+ (r*cos((m*pi)/18) - r*cos((n*pi)/18)).^2;
                    A(m,n) = 1 ./ (4*pi*a_mn) ;
                end
            end
        end %end - generate_A_q2d
        
        function y = Ly_b(L,b) %calc y by the equation Ly=b
            M = length(L);
            y = zeros(1,M);
            y(1) = b(1) / L(1,1);
            sum = 0;
            for j = 2:M
                for i = 1:j-1
                    sum = sum + L(j,i) .* y(i);
                end %for
                y(j) = (b(j) - sum) ./ L(j,j);
                sum = 0;
            end %for
        end %end - Ly_b

        function x = Ux_y(U,y)  %calc x by the equation Ux=y
            M = length(U);
            x = zeros(1,M);
            x(M) = y(M) / U(M,M);
            sum = 0;
            for i=M-1:-1:1
                for j=M:-1:i+1
                    sum = sum + U(i,j).*x(j);
                end %for
                x(i) = (y(i) - sum) ./ U(i,i);
                sum = 0;
            end %for
        end %end - Ux_y
        
        function [x] = solve_linear (L,b,U) %solve by LU decomp
            y = sub_functions.Ly_b(L,b);
            x = sub_functions.Ux_y(U,y);
        end %end - solve_linear
          
        function [Relative_distance,Relative_Error,iter_plot] = Gauss_Seidel(A, v, q, Error_const)
            L = tril(A,-1);
            D = diag(diag(A)); 
            Q = L + D ;
            neg_u = Q - A; %-U
            G = inv(Q) * neg_u  ; % (L+D)^(-1) * (-U)
            C =  inv(Q) * v;   % (L+D)^(-1) * v
            q_k = C;  
            Relative_distance = zeros;
            Relative_Error = zeros;
            i = 1;   
            iter_plot = zeros;
            error = max(abs(q-q_k));
            max_iter = 300; 
            while abs(error) > Error_const && i <= max_iter  
               q_k_prev = q_k; %q_k_prev = q_(k-1)
               q_k = G*(q_k_prev) + C; 
               error = norm(q-q_k, 'inf');
               Relative_distance(i) = norm(q_k - q_k_prev, 'inf') / norm(q_k_prev, 'inf');
               Relative_Error(i) = norm(q_k - q, 'inf') / norm(q, 'inf');
               iter_plot(i) = i;
               i = i + 1;
            end
        end %end - Gauss_Seidel
      
        function [Relative_Distance,Relative_Error,iter_plot] = Jacobi(A, v, q, Error_const)   
            D = diag(diag(A));
            Q = D ;
            I = eye(18);
            inv_Q = inv(Q);
            G = I - (inv_Q * A);
            norm_G = norm(G,'inf'); %inf norm smaller than 1
            C = inv_Q * v;
            q_k = C; % k=1 - q^(1)
            Relative_Distance = zeros;
            Relative_Error = zeros;
            i = 1;   
            iter_plot = zeros;
            error = max(abs(q-q_k));
            max_iter = 300; 
            while abs(error) > Error_const && i <=max_iter  
               q_k_prev = q_k; %q_k_prev = q_(k-1)
               q_k = G*(q_k_prev) + C; 
               error = norm(q-q_k, 'inf');
               Relative_Distance(i) = norm(q_k - q_k_prev, 'inf') / norm(q_k_prev, 'inf');
               Relative_Error(i) = norm(q_k - q, 'inf') / norm(q, 'inf');
               iter_plot(i) = i;
               i = i + 1;
            end
        end %end - Jacobi
          
    end %methods - line 3
end %classdef - line 2