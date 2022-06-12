%% Question 2:
M = 18;
q = [3,1,6,0,9,8,0,5,2,3,1,6,0,9,8,0,5,2]'; %id vector
h_set = [5,2,1];
h1 = pi ./ (5.*M);
h2 = pi ./ (2.*M);
h3 = pi ./ M;
Error_const = 10 ^ (-3);

%% 2.a)
A1 = sub_functions.generate_A(h1);
v1 = A1 * q;
[Relative_Distance_a, Relative_Error_a, Iteration_a] = sub_functions.Gauss_Seidel(A1, v1, q, Error_const); %h1 = pi/5M

    
figure(2);
subplot(2,3,1);
plot_A = semilogy(Iteration_a, Relative_Distance_a, Iteration_a, Relative_Error_a);
title('2.a - h=pi/5M (Gauss Seidel)');
xlabel('Iterations');
legend('Relative Destination of q^(k-1)and q^k' ,'Relative Real Error of  q^(k-1) and q^k','Location','southwest');

%% 2.b)
A2 = sub_functions.generate_A(h2);
v2 = A2 * q;
[Relative_distance_B, Relative_Error_B, Iteration_B] = sub_functions.Gauss_Seidel(A2, v2, q, Error_const); %h2 = pi/2M=

subplot(2,3,2);
plot_B = semilogy(Iteration_B,Relative_distance_B, Iteration_B, Relative_Error_B);
title('2.b - h=pi/2M (Gauss Seidel)');
legend('Relative Destination of q^(k-1)and q^k' ,'Relative Real Error of  q^(k-1) and q^k','Location','southwest');

A2b = sub_functions.generate_A(h3);
v2b = A2b * q;
[Relative_distance_B2, Relatove_Error_B2, Iteration_B2] = sub_functions.Gauss_Seidel(A2b, v2b, q, Error_const); %h=pi/M

subplot(2,3,3);
plot_B2 = semilogy(Iteration_B2, Relative_distance_B2, Iteration_B2, Relatove_Error_B2);
title('2.b - h=pi/M (Gauss Seidel)');
legend('Relative Destination of q^k and q^(k-1)','Relative Real Error of q^k and q^(k-1)');

%% 2.c) 
A3 = sub_functions.generate_A(h1);
v3 = A3 * q;
[Relative_distance_C, Relative_Error_C, Iteration_C] = sub_functions.Jacobi(A3, v3, q, Error_const); %h=pi/5M

subplot(2,3,4);
plot_C = semilogy(Iteration_C, Relative_distance_C, Iteration_C, Relative_Error_C);
plot_C(1).LineWidth = 2.5;
plot_C(2).LineWidth = 2.5;
title('2.c - h=pi/5M (Jacobi)');
legend('Relative Destination of q^k and q^(k-1)','Relative Real Error of q^k and q^(k-1)');

%% 2.d)
A4 = sub_functions.generate_A_q2d(h1);
v4 = A4 * q;
[Relative_distance_D, Relative_Error_D, Iteration_D] = sub_functions.Jacobi(A4, v4, q, Error_const); %h=pi/5M

subplot(2,3,6);
plot_D = semilogy(Iteration_D, Relative_distance_D, Iteration_D, Relative_Error_D);
plot_D(1).LineWidth = 2.5;
plot_D(2).LineWidth = 2.5;
title('2.d - h=pi/5M (Jacobi and new A)');
legend('Relative Destination of q^k and q^(k-1)','Relative Real Error of q^k and q^(k-1)','Location','southeast');
movegui('south');




