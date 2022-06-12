%% Question 3:
M = 18; %vector length
q = [3,1,6,0,9,8,0,5,2,3,1,6,0,9,8,0,5,2]'; %id vector
r = 1 ; %radius equals 1m
h_set = [10,5,2,0.5,0.2];
det_A = zeros(5,1);

%% 3.a) + 1.b)
Relative_Error_q = zeros(5,1);
for t=1:5
     h = h_set(t).*((pi * r) ./ M) ;
     A = sub_functions.generate_A(h);
     det_A(t) = abs(det(A));
     v = A * q; 
     A_trans = transpose(A);
     approx_q = inv(A_trans * A) * A_trans * v;
     Relative_Error_q(t) = (norm(approx_q - q)) / norm(q);
end
     
figure(3);
subplot(1,2,1);
graph_h = h_set .* (pi ./ 18);
l1 = loglog(graph_h, det_A);
title('3 - | det(A)(h) |');
xlabel('h');
legend('| det(A)(h) |','Location','southeast');
grid on;

subplot(1,2,2);
l2 = loglog(graph_h, Relative_Error_q);
title('3.a+b');
xlabel('h');
legend('3 - Relative Error of q & approx_q');
grid on;
movegui('east');