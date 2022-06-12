%% variables
a = 1;
b = 5;
I1 = 316098052;% my ID
tolerance = 10^(-12);
s = 3^(1/4); %Solution for Question 1 and 2
x_0 = a + (b - a) * (I1 / (I1 + I1)); %Initial guess //I1=I2
format long %long display for output


%% Question 1
[x_n_1, error_n, x_n_diff_1] = Newton_Raphson(x_0, tolerance, s, 1);

%Graph
figure(1)
y_axis = log(error_n(2 : end));  %log(epsilon_n)
x_axis = log(error_n(1 : end-1)); %log(epsilon_n_-1))
plot(x_axis, y_axis, '--o');
title('Question 1 - Newton Raphson');
ylabel('log(epsilon_n)');
xlabel('log(epsilon_n_-_1))'); 
grid on;
movegui('west');

%Required table
n = [1 : length(error_n)]';
X_n = x_n_1(:, 1 : length(x_n_1) - 1)';
X_n_diff = x_n_diff_1';
error = error_n';
T1 = table(n, X_n, X_n_diff, error);

disp(T1);


%% Question 2
x1 = x_0 + (b-x_0)*(I1/(I1+I1)); %Second Initial guess - I1=I2
[x_n_2, error_n, x_n_diff_2] = Secant(x_0, x1, tolerance, s);

%Graph
figure(2)
y_axis = log(error_n(2 : end)); %log(epsilon_n)
x_axis = log(error_n(1 : end-1)); %log(epsilon_n_-1)
plot(x_axis, y_axis, '--o');
title('Question 2 - Secant');
ylabel('log(epsilon_n)');
xlabel('log(epsilon_n_-_1))'); 
grid on;
movegui('north');

%Required table
n = [1:length(error_n)]';
X_n = x_n_2(:,1:length(x_n_2)-1)';
X_n_diff = x_n_diff_2';
error = error_n';
T2 = table(n,X_n, X_n_diff, error);

disp(T2);


%% Question 3
% Part A
x_0 = 5;
s3 = 2;
[x_n_3A, error_n, x_n_diff_3A] = Newton_Raphson_multiple_roots(x_0, 1, tolerance, s3, 1);

%Graph
figure(3)
subplot(2,2,1)
y_axis = log(error_n(2 : end)); %log(epsilon_n)
x_axis = log(error_n(1 : end-1)); %log(epsilon_n_-1)
plot(x_axis, y_axis, '-o');
title('Question 3A - Newton Raphson with multiple roots');
ylabel('log(epsilon_n)');
xlabel('log(epsilon_n_-_1))'); 
grid on;;

%Required table
n = [1 : length(error_n)]';
X_n = x_n_3A(:, 1 : length(x_n_3A) - 1)';
X_n_diff = x_n_diff_3A';
error = error_n';
T3A = table(n, X_n, X_n_diff, error);

disp(T3A); %Display the table in the command window.

% Part B
[X_n_3B, error_n, Xn_diff_3B] = Newton_Raphson_multiple_roots(x_0, 1, tolerance, s3, 2);

%Graph
subplot(2,2,2)
y_axis = log(error_n(2 : end)); %log(epsilon_n)
x_axis = log(error_n(1 : end-1)); %log(epsilon_n_-1)
plot(x_axis, y_axis, '--o');
title('Question 3B - Newton Raphson with multiple roots');
ylabel('log(epsilon_n)');
xlabel('log(epsilon_n_-_1))'); 
grid on;;

%Required table
n = [1 : length(error_n)]';
X_n = X_n_3B(:, 1 : length(X_n_3B) - 1)';
X_n_diff = Xn_diff_3B';
error = error_n';
T3B = table(n, X_n, X_n_diff, error);

disp(T3B);

% Part C
[X_n_3C, error_n, Xn_diff_3C] = Newton_Raphson_multiple_roots(x_0,2.999,tolerance,s3, 1); %the reason for 2.999 explained in the report

%Graph
subplot(2,2,3)
y_axis = log(error_n(2 : end)); %log(epsilon_n)
x_axis = log(error_n(1 : end-1)); %log(epsilon_n_-1)
plot(x_axis, y_axis, '--o');
title('Question 3C - Newton Raphson with multiple roots');
ylabel('log(epsilon_n)');
xlabel('log(epsilon_n_-_1))'); 
grid on;
movegui('east');

%Required table
n = [1:length(error_n)]';
X_n = X_n_3C(:,1:length(X_n_3C)-1)';
X_n_diff = Xn_diff_3C';
error = error_n';
T3C = table(n,X_n, X_n_diff, error);

disp(T3C);


%% Question 4
% Part A
x_0 = pi/2;
s4 = 1.895494267034;
[X_n_4A, error_n, Xn_diff_4A] = Fixed_Point(x_0, tolerance, s4, 1);

%Graph
figure(4)
subplot(2,2,1)
y_axis = log(error_n(2 : end)); %log(epsilon_n)
x_axis = log(error_n(1 : end-1)); %log(epsilon_n_-1)
plot(x_axis, y_axis, '-o');
title('Question 4A - Fixed Point');
ylabel('log(epsilon_n)');
xlabel('log(epsilon_n_-_1))'); 
grid on;

%Required table
n = [1 : length(error_n)]'; 
X_n = X_n_4A(:, 1 : length(X_n_4A) - 1)'; 
X_n_diff = Xn_diff_4A'; 
error = error_n'; 
T4A = table(n, X_n, X_n_diff, error);

disp(T4A);

% Part B
[X_n_4B, error_n, Xn_diff_4B] = Newton_Raphson(x_0, tolerance, s4, 4);

%Graph
subplot(2,2,2)
y_axis = log(error_n(2 : end)); %log(epsilon_n)
x_axis = log(error_n(1 : end-1)); %log(epsilon_n_-1)
plot(x_axis, y_axis, '--o');
title('Question 4B - Newton Raphson');
ylabel('log(epsilon_n)');
xlabel('log(epsilon_n_-_1))'); 
grid on;

%Required table
n = [1 : length(error_n)]'; 
X_n = X_n_4B(:, 1 : length(X_n_4B) - 1)';
X_n_diff = Xn_diff_4B';
error = error_n';
T4B = table(n, X_n, X_n_diff, error);

disp(T4B);

% Part D
s4D = 0;
x_0 = 1/2;
[X_n_4D, error_n, Xn_diff_4D] = Fixed_Point(x_0, tolerance, s4D, 2);

%Graph
subplot(2,2,3)
y_axis = log(error_n(2 : end)); %log(epsilon_n)
x_axis = log(error_n(1 : end-1)); %log(epsilon_n_-1)
plot(x_axis, y_axis, '-o');
title('Question 4D - Fixed Point');
ylabel('log(epsilon_n)');
xlabel('log(epsilon_n_-_1))'); 
grid on;
movegui('south');

%Required table
n = [1:length(error_n)]';
X_n = X_n_4D(:,1:length(X_n_4D)-1)';
X_n_diff = Xn_diff_4D';
error = error_n';
T4D = table(n,X_n, X_n_diff, error);

disp(T4D);


%% Sub-Functions
%Newton-Raphson
function [x_n, error_n, x_n_diff] = Newton_Raphson(x_0, tolerance, s, question)
    iteration = 2;
    error_n = zeros; 
    x_n_diff = zeros;
    x_n(1) = x_0; %initial guess
    if question == 1
        fx_div_fx_tag = (x_0^4 - 3) / (4*x_0^3);
    elseif question == 4
        fx_div_fx_tag = (x_0 - 2*sin(x_0)) / (1 - 2*cos(x_0));
    end
    
    x_n(2) = x_0 - fx_div_fx_tag;
    x_n_diff_abs = abs(x_n(iteration) - x_n(iteration - 1)); %|x_n - x_n_-1)|
    x_n_diff(1) = x_n_diff_abs;
    error_n(1) = abs(x_n(1) - s); %|x_n - s|
    while x_n_diff_abs >= tolerance  
        if question == 1
            fx_div_fx_tag = (x_n(iteration)^4 - 3) / (4*x_n(iteration)^3);
        elseif question == 4
            fx_div_fx_tag = (x_n(iteration) - 2*sin(x_n(iteration))) / (1 - 2*cos(x_n(iteration)));
        end
        x_n(iteration+1) = x_n(iteration) - fx_div_fx_tag;
        error_n(iteration) = abs(x_n(iteration) - s);     
        iteration = iteration + 1;    
        x_n_diff_abs = abs(x_n(iteration) - x_n(iteration - 1));  
        x_n_diff(iteration-1) = x_n_diff_abs;       
    end
end

%Secant
function [x_n, error_n, x_n_diff] = Secant(x0, x1, tolerance, s)
    iteration = 2;
    error_n = zeros;
    x_n_diff = zeros;
    x_n(1) = x0; %initial guess sol.
    x_n(2) = x1; %Second guess sol.
    x_n_diff_abs = abs(x_n(iteration) - x_n(iteration-1)); %|x_n - x_n_-1|
    x_n_diff(1) = x_n_diff_abs;
    error_n(1) = abs(x_n(1) - s); 
    while abs(x_n(iteration) - x_n(iteration-1)) > tolerance
        x_n(iteration + 1) = x_n(iteration) - ((x_n(iteration))^4 - 3) * (x_n(iteration) - x_n(iteration - 1)) / ((x_n(iteration)^4 - 3) - (x_n(iteration - 1)^4 - 3)); %Iteration formula
        error_n(iteration) = abs(x_n(iteration) - s);  % |x_n -s|
        iteration = iteration + 1;    
        x_n_diff_abs = abs(x_n(iteration) - x_n(iteration - 1));  
        x_n_diff(iteration - 1) = x_n_diff_abs;
    end
end

%Newton-Raphson with multiple roots
function [x_n, error_n, x_n_diff] = Newton_Raphson_multiple_roots(x_0, q, tolerance, s, option)
    %fx = x^5 - 6*x^4 + 14*x^3 - 20*x^2 + 24*x - 16;
    %fx_tag = 5*x^4 - 24*x^3 + 42*x^2 - 40*x + 24;
    %fx_double_tag = 20*x^3 - 72*x^2 + 84*x - 40;
    iteration = 2;
    error_n = zeros; 
    x_n_diff = zeros;
    x_n(1) = x_0;
    u_x0 = (x_0^5 - 6*x_0^4 + 14*x_0^3 - 20*x_0^2 + 24*x_0 - 16) / (5*x_0^4 - 24*x_0^3 + 42*x_0^2 - 40*x_0 + 24); %u(x_0) = f(x_0)/f'(x_0)
    y = 0;
    if option == 1
        y = u_x0; % y = u(x_0)
    elseif option == 2
        fx = x_0^5 - 6*x_0^4 + 14*x_0^3 - 20*x_0^2 + 24*x_0 - 16;
        fx_tag = 5*x_0^4 - 24*x_0^3 + 42*x_0^2 - 40*x_0 + 24;
        fx_double_tag = 20*x_0^3 - 72*x_0^2 + 84*x_0 - 40;
        u_x0_tag = 1 - (fx*fx_double_tag) / (fx_tag^2);
        y = (u_x0) / u_x0_tag; % y = u(x_0)/u'(x_0)
    end
    x_n(2) = x_0 - y;
    x_n_diff_abs = abs(x_n(iteration) - x_n(iteration - 1)); %|x_n - x_n_-1|
    x_n_diff(1) = x_n_diff_abs;
    error_n(1) = abs(x_n(1) - s);
    while (x_n_diff_abs > tolerance) && (abs(x_n(iteration) - s) < (abs(x_n(iteration-1) - s))) 
        u_xn = (x_n(iteration)^5 - 6*x_n(iteration)^4 + 14*x_n(iteration)^3 - 20*x_n(iteration)^2 + 24*x_n(iteration) - 16) / (5*x_n(iteration)^4 - 24*x_n(iteration)^3 + 42*x_n(iteration)^2 - 40*x_n(iteration) + 24); %u(x_n) = f(x_n)/f'(x_n)
        if option == 1
            y = u_xn; % y = u(x_n)
        elseif option == 2
            fx = x_n(iteration)^5 - 6*x_n(iteration)^4 + 14*x_n(iteration)^3 - 20*x_n(iteration)^2 + 24*x_n(iteration) - 16;
            fx_tag = 5*x_n(iteration)^4 - 24*x_n(iteration)^3 + 42*x_n(iteration)^2 - 40*x_n(iteration) + 24;
            fx_double_tag = 20*x_n(iteration)^3 - 72*x_n(iteration)^2 + 84*x_n(iteration) - 40;
            u_xn_tag = 1 - (fx*fx_double_tag) / (fx_tag^2);
            y = (u_xn) / u_xn_tag; % y = u(x_n)/u'(x_n)
        end
        x_n(iteration+1) = x_n(iteration) - q * y;
        error_n(iteration) = abs(x_n(iteration) - s);  % |x_n -s|   
        iteration = iteration + 1;    
        x_n_diff_abs = abs(x_n(iteration) - x_n(iteration - 1));  
        x_n_diff(iteration - 1) = x_n_diff_abs; 
    end   
end

%Fixed Point
function [x_n, error_n, x_n_diff] = Fixed_Point(x_0, tolerance, s , option)
    iteration = 2;
    error_n = zeros;
    x_n_diff = zeros;
    x_n(1) = x_0;
    g = 0;
    if option == 1
        g = 2*sin(x_n(1));
    elseif option == 2
        g = asin(x_n(1) / 2);
    end
    x_n(2) = g;
    while abs(x_n(iteration) - x_n(iteration - 1)) >= tolerance
        if option == 1
            g = 2*sin(x_n(iteration));
        elseif option == 2
            g = asin(x_n(iteration) / 2);
        end
        x_n(iteration + 1) = g;
        error_n(iteration) = abs(x_n(iteration) - s);     
        iteration = iteration + 1;    
        Xn_Diff = abs(x_n(iteration)-x_n(iteration-1));  
        x_n_diff(iteration-1) = Xn_Diff;
    end
end

