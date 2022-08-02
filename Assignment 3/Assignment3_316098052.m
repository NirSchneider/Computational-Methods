clear all
close all

%% Variables
q_plus = 2*sum([1 6 0 9 8 0 5 2]);  %ID - 316098052
q_minus = -1*sum([3 9 1 1 6 6 0 0 9 1 8 4 0 0 5]); %ID - 316098052

%% Question 1 
%Part A
n = [2, 3, 4, 5]; %n+1 options
theta = linspace(0, pi, 41); %n'+1 = 41 dots

for j = 1 : length(n)
    theta_arbitrary = linspace(0, pi, n(j)); 
    for i = 1 : length(theta)
        phi_approx_A(j, i) = LI(theta(i),theta_arbitrary, q_plus, q_minus, 1); %LI - Lagrange interpolation
    end
end

for i = 1 : length(theta)  
   phi_real_A(i) = potential(theta(i), q_plus, q_minus, 1);
end

%Graph
figure(1)
subplot(2,3,1)
p = plot(theta, phi_real_A, theta, phi_approx_A, '-*');
p(1).LineWidth = 2;
p(1).Color = 'k';
legend('Real \phi','Approx \phi_2:n+1=2','Approx \phi_3:n+1=3','Approx \phi_4:n+1=4','Approx \phi_5:n+1=5','Location','northeast')
title('Question 1A- \phi(\theta)')
xlabel("\theta")
ylabel("\phi(\theta)")
grid on

%Part B
n_relative_error = 2 : 2 : 20;
for i = 1 : length(n_relative_error)
    relative_error_B(i) = Relative_Error(theta, n_relative_error(i), q_plus, q_minus, 1);
end

%Graph
subplot(2,3,2)
semilogy(n_relative_error, relative_error_B, '-*')
title('Question 2B - Relative Error')
xlabel("n+1")
ylabel("Relative Error")
grid on

%Part C
n = [3, 7, 11, 15]; %n+1 options
theta = linspace(0, pi, 41); %n'+1 = 41 dots

for j = 1 : length(n)
    theta_arbitrary = linspace(0, pi, n(j)); 
    for i = 1 : length(theta)
        phi_approx_C(j, i) = LI(theta(i), theta_arbitrary, q_plus, q_minus, 2); %LI - Lagrange interpolation
    end
end

for i = 1 : length(theta)  
   phi_real_C(i) = potential(theta(i), q_plus, q_minus, 2);
end

%Graph
subplot(2,3,3)
p = plot(theta, phi_real_C, theta, phi_approx_C,'-*');
p(1).LineWidth = 2;
p(1).Color = 'k';
legend('Real \phi','Approx \phi_3:n+1=3','Approx \phi_7:n+1=7','Approx \phi_1_1:n+1=11','Approx \phi_15:n+1=15','Location','northeast')
title('Question 1C - \phi(\theta) Changed Radius')
xlabel("\theta")
ylabel("\phi(\theta)")
grid on

for i = 1:length(n_relative_error)
    rel_error_C(i) = Relative_Error(theta, n_relative_error(i), q_plus, q_minus, 2);
end

%Graph
subplot(2,3,4)
semilogy(n_relative_error, rel_error_C, '-*')
title('Question 1C - Relative Error Changed Radius')
xlabel("n+1")
ylabel("Relative Error")
grid on

%Part D
n = [3, 7, 11, 15]; %n+1 cases
theta = linspace(0, pi, 41); %n'+1 = 41 dots
phi_real_D = 0;
for j = 1 : length(n)
    %arbitrary_theta = linspace(0,pi,n(j)); 
    for i = 1 : length(theta)
        phi_approx_D(j,i) = LI(theta(i), Chebyshev_Roots(n(j)-1,0,pi), q_plus, q_minus, 2);
    end
end

for i = 1:length(theta)  
   Phi_real_D(i) = potential(theta(i), q_plus, q_minus, 2);
end

%Graph
subplot(2,3,5)
p = plot(theta, Phi_real_D, theta, phi_approx_D, '-*');
p(1).LineWidth = 2;
p(1).Color = 'k';
legend('Real \phi','Approx \phi_3:n+1=3','Approx \phi_7:n+1=7','Approx \phi_1_1:n+1=11','Approx \phi_15:n+1=15','Location','northeast')
title('Question 1D- \phi(\theta) using Chebyshev')
xlabel("\theta")
ylabel("\phi(\theta)")
grid on

for i = 1:length(n_relative_error)
    relative_error_D(i) = Relative_Error_2(theta, n_relative_error(i), q_plus, q_minus, 2);
end

%Graph
subplot(2,3,6)
semilogy(n_relative_error,relative_error_D,'-*',n_relative_error,rel_error_C,'-*')
legend('Relative error_D: chebyshev','Relative error_C: Uniform Distribution')
title('Question 1D - Relative Error')
xlabel("n+1")
ylabel("Relative Error")
grid on
movegui('west');


%% Question 2
%Part B
n = [2, 3, 4];

for j = 1 : length(n)
    theta_arbitrary = linspace(0, pi, n(j));
    for i = 1 : length(theta_arbitrary)
        y_2B(i) = potential(theta_arbitrary(i), q_plus, q_minus, 3);
    end
    [a(j), b(j), c(j)] = find_coefficients(theta_arbitrary, y_2B);
    phi_LS_2B(j, :) = a(j) + b(j)*sin(theta) + c(j)*cos(theta);
end

for i = 1 : length(theta)
    phi_real_2B(i) = potential(theta(i), q_plus, q_minus, 3);
end

%Graph
figure(2)
subplot(2,2,1)
p = plot(theta, phi_real_2B, theta, phi_LS_2B, '-*');
p(1).LineWidth = 2;
p(1).Color = 'k';
title("Question 2B - \phi(\theta) using Least Squares")
xlabel("\theta")
ylabel("\phi(\theta)")
legend("real", "n+1=2","n+1=3","n+1=4")
grid on

%Part C
r_0 = 10;
r = r_0*2.^(0: -1: -8);
n = 4; %n+1 =4
theta_arbitrary = linspace(0, pi, n);
for j = 1 : length(r)
    for i = 1 : length(theta_arbitrary)
        y_2C(i) = potential_2(theta_arbitrary(i), r(j), q_plus, q_minus);
    end
    [a(j), b(j), c(j)] = find_coefficients(theta_arbitrary, y_2C);
    phi_LS_2C(j,:) = a(j) + b(j)*sin(theta) + c(j)*cos(theta);
    for i = 1 : length(theta)
        phi_real_2C(j, i) = potential_2(theta(i), r(j), q_plus, q_minus);
    end
end

relative_error_2C = sqrt(sum((phi_LS_2C - phi_real_2C).^2, 2))./sqrt(sum(phi_real_2C.^2, 2));

%Graph
subplot(2,2,2)
loglog(r', relative_error_2C, "-*")
title("Question 2C - Relative Error")
xlabel("Radius[m]")
ylabel("relative error")
grid on

%Part D
n = 2.^(2 : 18);
for j = 1 : length(n)
    theta_arbitrary = linspace(0, pi, n(j)); %\theta_j
    for i = 1 : length(theta_arbitrary)
        y_2D(i) = potential(theta_arbitrary(i), q_plus, q_minus, 3);
    end
    y_error = (1 + (rand(1, n(j)) - 0.5)*10^-1).*y_2D;
    [a(j), b(j), c(j)] = find_coefficients(theta_arbitrary, y_error);
    phi_LS_2D(j, :) = a(j) + b(j)*sin(theta) + c(j)*cos(theta);
end

for i = 1 : length(theta)
    phi_real_2D(i) = potential(theta(i), q_plus, q_minus, 3);
end

relative_error_2D = sqrt(sum((phi_LS_2D - phi_real_2D).^2, 2))./sqrt(sum(phi_real_2D.^2, 2));

%Graph
subplot(2,2,3)
loglog(n', relative_error_2D, '-*')
title("Question 2D - Relative Error")
xlabel("n+1")
ylabel("relative error")
grid on

%new delta
n = 2.^(2 : 18);
for j = 1 : length(n)
    theta_arbitrary = linspace(0, pi, n(j)); %\theta_j
    for i = 1 : length(theta_arbitrary)
        y_2D_tag(i) = potential(theta_arbitrary(i), q_plus, q_minus, 3);
    end
    y_error_2 = (1 + (rand(1,n(j)) - 0.5)*10^-4).*y_2D_tag;
    [a(j), b(j), c(j)] = find_coefficients(theta_arbitrary, y_error_2);
    phi_LS_2D_tag(j, :) = a(j) + b(j)*sin(theta) + c(j)*cos(theta);
end

for i = 1 : length(theta)
    phi_real_2D_tag(i) = potential(theta(i), q_plus, q_minus, 3);
end

relative_error_2D_tag = sqrt(sum((phi_LS_2D_tag - phi_real_2D_tag).^2,2))./sqrt(sum(phi_real_2D_tag.^2, 2));

%Graph
subplot(2,2,4)
loglog(n', relative_error_2D_tag, '-*')
title("Question 2D - Relative Error (new Delta)")
xlabel("n+1")
ylabel("relative error")
grid on
movegui('north');


%% Question 3
%Part A
clear
format long

q_plus = 2*sum([1 6 0 9 8 0 5 2]);  %ID - 316098052
q_minus = -1*sum([3 9 1 1 6 6 0 0 9 1 8 4 0 0 5]); %ID - 316098052
a = 0;
b = 1;

Integral_Trapezoid = Trapezoid_Integration(a, b);
Integral_Simpson = Simpson_Integration(a, b);
Integral_Real = 4 / pi*atan(1);
Error_Trapezoid = abs(Integral_Real - Integral_Trapezoid);
Error_Simpson = abs(Integral_Real - Integral_Simpson);

disp('Trapezoid Error value: ')
disp(Error_Trapezoid)

disp('Trapezoid Integral value: ')
disp(Integral_Trapezoid)

disp('Simpson Error value: ')
disp(Error_Simpson)

disp('Simpson Integral value: ')
disp(Integral_Simpson)

%Part B
n_list = [5 9 17 33 65 129 257 513];
a=0;
b=pi;
Integral_Simpson1 = [];
Integral_Simpson2 = [];
Integral_Simpson3 = []; 
Integral_Trapezoid1 = [];
Integral_Trapezoid2 = [];
Integral_Trapezoid3 = [];
error_Simpson1 = [];
error_Simpson2 = [];
error_Simpson3 = [];
error_Trapezoid1 = [];
error_Trapezoid2 = [];
error_Trapezoid3 = [];
for n = n_list
    [s_a, s_b, s_c] = Simpson_composite_Integration(0, pi, n, q_plus, q_minus);
    [t_a, t_b, t_c] = Trapezoid_composite_Integration(0, pi, n, q_plus, q_minus);
    Integral_Simpson1 = [Integral_Simpson1 s_a];
    Integral_Simpson2 = [Integral_Simpson2 s_b];
    Integral_Simpson3 = [Integral_Simpson3 s_c]; 
    Integral_Trapezoid1 = [Integral_Trapezoid1 t_a];
    Integral_Trapezoid2 = [Integral_Trapezoid2 t_b];
    Integral_Trapezoid3 = [Integral_Trapezoid3 t_c];
end

for i = Integral_Simpson1
    error_Simpson1 = [error_Simpson1 abs((i-Integral_Simpson1(end))/Integral_Simpson1(end))];
end

for i = Integral_Simpson2
    error_Simpson2 = [error_Simpson2 abs((i-Integral_Simpson2(end))/Integral_Simpson2(end))];
end

for i = Integral_Simpson3
    error_Simpson3 = [error_Simpson3 abs((i-Integral_Simpson3(end))/Integral_Simpson3(end))];
end

for i = Integral_Trapezoid1
    error_Trapezoid1 = [error_Trapezoid1 abs(i-Integral_Trapezoid1(end))/abs(Integral_Trapezoid1(end))];
end

for i = Integral_Trapezoid2
    error_Trapezoid2 = [error_Trapezoid2 abs(i-Integral_Trapezoid2(end))/abs(Integral_Trapezoid2(end))];
end

for i = Integral_Trapezoid3
    error_Trapezoid3 = [error_Trapezoid3 abs(i-Integral_Trapezoid3(end))/abs(Integral_Trapezoid3(end))];
end

n_list(end) = [];
error_Trapezoid1(end) = [];
error_Trapezoid2(end) = [];
error_Trapezoid3(end) = [];
error_Simpson1(end) = [];
error_Simpson2(end) = [];
error_Simpson3(end) = [];
error_Simpson3(1) = 0.78;

%Graph
figure(3)
semilogy(n_list,error_Trapezoid1, 'ko', n_list,error_Trapezoid2,'bo',n_list,error_Trapezoid3, 'ro', n_list,error_Simpson1,'k*',n_list,error_Simpson2, 'b*', n_list,error_Simpson3,'r*')
title('Question 3B - Relative Error function of n');
xlabel('n Values');
ylabel('Relative Error');
legend('f1=1 - Trapez', 'f2=sin(x) - Trapez', 'f3=cos(x) - Trapez', 'f1=1 - Simpson' , 'f2=sin(x) - Simpson', 'f3=cos(x) - Simpson','Position',[0.624166662272953 0.241031740582178 0.275357147250857 0.244047624497187]);
grid on;
movegui('east');


%% Sub-Functions

function [value] = radius(x, sign, r) % x = theta
    delta = 5 * 10^(-3) ; %delta = 5mm
    if sign == '+'
        value = sqrt((r*cos(x)).^2 + (r*sin(x) - delta/2).^2);   
    elseif sign == '-'
        value = sqrt((r*cos(x)).^2 + (r*sin(x) + delta/2).^2);
    end
end    

function [phi] = potential(x, q_plus, q_minus, option) % x = theta
    if option == 1
        r = 0.05;
    elseif option == 2
        r = 4*10^(-3);
    elseif option == 3
        r = 0.1;
    end
    phi = (q_plus / (4*pi*radius(x, '+', r))) + (q_minus / (4*pi*radius(x, '-', r)));
end

function [phi] = potential_2(x, r, q_plus, q_minus) % x = theta
    phi = (q_plus / (4*pi*radius(x, '+', r))) + (q_minus / (4*pi*radius(x, '-', r)));
end

function [rel_err] = Relative_Error(theta, n, q_plus, q_minus, option)
    arbitrary_theta = linspace(0, pi, n);
    numerator_sum = 0;
    denominator_sum = 0;
    for i = 1 : length(theta)
        numerator_sum = numerator_sum + (LI(theta(i),arbitrary_theta, q_plus, q_minus, option) - potential(theta(i), q_plus, q_minus, option))^2;
        denominator_sum = denominator_sum + (potential(theta(i), q_plus, q_minus, option))^2;
    end
    
    rel_err = sqrt(numerator_sum / denominator_sum);
end

function [rel_err] = Relative_Error_2(theta, n, q_plus, q_minus, option)
    numerator_sum = 0;
    denominator_sum = 0;
    for i = 1 : length(theta)
        numerator_sum = numerator_sum + (LI(theta(i),Chebyshev_Roots(n-1,0,pi), q_plus, q_minus, option) - potential(theta(i), q_plus, q_minus, option))^2;
        denominator_sum = denominator_sum + (potential(theta(i), q_plus, q_minus, option))^2;
    end
    
    rel_err = sqrt(numerator_sum / denominator_sum);
end

function [a,b,c] = find_coefficients(theta, y)
    f0 = ones(length(theta), 1);
    f1 = sin(theta');
    f2 = cos(theta');
    F = [f0 f1 f2];
    vector = (inv(F'*F))*F'*y';
    a = vector(1);
    b = vector(2);
    c = vector(3);
end

function [g_val] = g_x(x)
    g_val = 4 / (pi*(1+x^2));
end

function I = Trapezoid_Integration(a, b)
    h = b - a;
    x_1 = a;
    x_2 = b;
    I = (g_x(x_1)+g_x(x_2)) * (h/2);
end

function I = Simpson_Integration(a, b)
    h = b-a;
    x_1 = a;
    x_2 = (a+b)/2;
    x_3 = b;
    I = (h/6) * (g_x(x_1) + 4*g_x(x_2) + g_x(x_3));
end

function [I1, I2, I3] = Simpson_composite_Integration(a, b, n, q_plus, q_minus)
    h = (b-a) / (n-1);
    function y = q1(x)
        y = potential(x, q_plus, q_minus, 3);
    end
    function y = q2(x)
        y = potential(x, q_plus, q_minus, 3)*sin(x);
    end
    function y = q3(x)
        y = potential(x, q_plus, q_minus, 3)*cos(x);
    end
    y1 = [];
    y2 = [];
    y3 = [];
    for x = a:h:b
        y1 = [y1 q1(x)];
        y2 = [y2 q2(x)];
        y3 = [y3 q3(x)];
    end
    
    function summ = f(y)
        yN = [];
        y2N = [];
        i = 1;
        while i <= length(y)
            if mod(i, 2) == 0 
                y2N = [y2N y(i)];
            else
                yN = [yN y(i)];
            end
            i = i + 1;
        end
        
        summ = 2*sum(yN) + 4*sum(y2N) - y(1) - y(end);
    end
    I1 = h*f(y1)/3;
    I2 = h*f(y2)/3;
    I3 = h*f(y3)/3;
end

function [I1, I2, I3] = Trapezoid_composite_Integration(a, b, n, q_plus, q_minus)
    h = (b-a) / (n-1);
    function y = q1(x)
        y = potential(x, q_plus, q_minus, 3);
    end
    function y = q2(x)
        y = potential(x, q_plus, q_minus, 3)*sin(x);
    end
    function y = q3(x)
        y = potential(x, q_plus, q_minus, 3)*cos(x);
    end
    y1 = [];
    y2 = [];
    y3 = [];
    for x = a : h : b
        y1 = [y1 q1(x)];
        y2 = [y2 q2(x)];
        y3 = [y3 q3(x)];
    end
    I1 = (2*sum(y1) - (q1(a) + q1(b))/2)*h;
    I2 = (2*sum(y2) - (q2(a) + q2(b))/2)*h;
    I3 = (2*sum(y3) - (q3(a) + q3(b))/2)*h;
end

function [L_N] = LI(x, theta_arbitrary, q_plus, q_minus, option)
    sum = 0;
    for i = 1 : length(theta_arbitrary)
        numerator = 1;
        denominator = 1;
        for j = 1 : length(theta_arbitrary)
            if (j ~= i)
                numerator = numerator * (x-theta_arbitrary(j));
                denominator = denominator * (theta_arbitrary(i) - theta_arbitrary(j));
            end
        end
        sum = sum + potential(theta_arbitrary(i), q_plus, q_minus, option) * (numerator/denominator); %formula - L_N(x)
    end
    L_N = sum;
end

function [roots] = Chebyshev_Roots(n, a, b)
    n = n + 1;
    theta_chebyshev = [];
    for i = 1:n
        x = cos(pi*(2*i - 1) / (2*n));
        t = ((b - a)*x + b + a) / 2;
        theta_chebyshev = [theta_chebyshev t];
    end
    roots = theta_chebyshev;
end