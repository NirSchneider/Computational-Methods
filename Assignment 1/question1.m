%% Question 1:
M = 18; %vector length
q = [3,1,6,0,9,8,0,5,2,3,1,6,0,9,8,0,5,2]'; %id vector
r = 1 ; %radius equals 1m
h_set = [1,2,5,10,20,50];
kappa_A = zeros(1,6)'; %k(A) - condition number of A
q2_Relative_Error = zeros(1,6)'; %Relative Error for q' for 1.b
q3_Relative_Error = zeros(1,6)'; %Relative Error for q for 1.c
q4_Relative_Error = zeros(1,6)'; %Relative Error for q for 1.d

%% 1.a)
for i = 1:6   
    h = h_set(i) .* (pi * r ./ 18);
    A = sub_functions.generate_A(h);
    v = A * q;
    [L,U,P] = lu(A);
    kappa_A(i) = norm(A,inf) .* norm(inv(A),inf);
    norm_frobenius_A = norm(A,'fro');
    norm_q = norm(q,2);
    norm_v = norm(v,2);
%% 1.b)
    vector_x = sub_functions.solve_linear(L,P*v,U);
    q2 = vector_x';
    norm_q2 = norm(q2 - q,2);
    q2_Relative_Error(i) = norm_q2 ./ norm_q;
%% 1.c)
    delta_v = norm_v .* 10.^(-3);
    v3 = v + delta_v;
    vector_x = sub_functions.solve_linear(L,P * v3,U);
    q3 = vector_x';
    norm_q3 = norm(q3 - q,2);
    q3_Relative_Error(i) = norm_q3 ./ norm_q;
%% 1.d)
    delta_A = norm_frobenius_A .* 10^(-3);
    A2 = A + delta_A;
    [L,U,P] = lu(A2);
    vector_x = sub_functions.solve_linear(L,P * v,U);
    q4 = vector_x';
    norm_q4 = norm(q4 - q,2);
    q4_Relative_Error(i) = norm_q4 ./ norm_q;
end
%% 1.e)
    figure(1);
    subplot(1,2,1);
    graph_h = h_set .* (pi ./ 18);
    l1 = loglog(graph_h, kappa_A);
    title("1.e");
    ylabel('Condition Number');
    xlabel('h(k)');
    legend('k(A)[h]','Location','southeast');
    grid on

    subplot(1,2,2);
    l2 = loglog(graph_h, q2_Relative_Error, graph_h, q3_Relative_Error, graph_h, q4_Relative_Error);
    title("1.e");
    ylabel('Relative Error ');
    xlabel('h(k)');
    legend('Relative Error','Relative Error_v','Relative Error_A','Location','southeast');
    grid on
    movegui('west');



    
