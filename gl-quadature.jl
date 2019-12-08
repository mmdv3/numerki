using Printf
### Polynomial power base / polynomials
p_0(x) = 1.0;
p_1(x) = x;
p_2(x) = x*x;
p_3(x) = x*x*x;
p_4(x) = x*x*x*x;
p_5(x) = x*x*x*x*x;
p_6(x) = x*x*x*x*x*x;
p_7(x) = x*x*x*x*x*x*x;


### Gauss-Lagrange quadrature over [-1,1] interval
function base_q_1(f)
	return 2.0*f(0.0);
end

function base_q_2(f)
x_1 = -1*sqrt(3.0)/3.0;
x_2 = sqrt(3.0)/3.0;

w_1 = 1;
w_2 = 1;
	return w_1*f(x_1)+w_2*f(x_2);
end

function base_q_3(f)
x_1 = -1*sqrt(15.0)/5.0;
x_2 = 0;
x_3 = sqrt(15.0)/5.0;

w_1 = 5.0/9.0;
w_2 = 8.0/9.0;
w_3 = 5.0/9.0;
	return w_1*f(x_1)+w_2*f(x_2)+w_3*f(x_3);
end

function base_q_4(f)
x_1 = -1*sqrt(525.0 - 70*sqrt(30.0))/35.0;
x_2 = -1*sqrt(525.0 + 70*sqrt(30.0))/35.0;
x_3 = sqrt(525.0 + 70*sqrt(30.0))/35.0;
x_4 = sqrt(525.0 - 70*sqrt(30.0))/35.0;

w_1 = (18.0+sqrt(30.0))/36.0;
w_2 = (18.0-sqrt(30.0))/36.0;
w_3 = (18.0-sqrt(30.0))/36.0;
w_4 = (18.0+sqrt(30.0))/36.0;
	return w_1*f(x_1)+w_2*f(x_2)+w_3*f(x_3)+w_4*f(x_4);
end

### (Shifted) Quadrature over [0,1] interval
function shift(f)
	function out(x)
		return 0.5*f((x+1.0)*0.5);
	end
	return out
end

function shifted_q_1(f)
	return base_q_1(shift(f))
end

function shifted_q_2(f)
	return base_q_2(shift(f))
end

function shifted_q_3(f)
	return base_q_3(shift(f))
end

function shifted_q_4(f)
	return base_q_4(shift(f))
end

### double check
function print_compare(label, polynomial, quadrature, integral_value)
	@printf("\t");
	@printf("polynomial:%s, quadrature value:%.8f, integral value:%.8f\n",
			label, quadrature(polynomial), integral_value);
end

@printf("Interval [-1,1]\n");
	@printf("\t quadrature degree 1:\n");
	print_compare("1", 	p_0, 	base_q_1, 2.0);
	print_compare("x", 	p_1, 	base_q_1, 0.0);

	@printf("\t quadrature degree 2:\n");
	print_compare("1", 	p_0, 	base_q_2, 2.0);
	print_compare("x", 	p_1, 	base_q_2, 0.0);
	print_compare("x^2", p_2, 	base_q_2, 2.0/3.0);
	print_compare("x^3", p_3, 	base_q_2, 0);

	@printf("\t quadrature degree 3:\n");
	print_compare("1", 	p_0, 	base_q_3, 2.0);
	print_compare("x", 	p_1, 	base_q_3, 0.0);
	print_compare("x^2", p_2, 	base_q_3, 2.0/3.0);
	print_compare("x^3", p_3, 	base_q_3, 0);
	print_compare("x^4", p_4, 	base_q_3, 2.0/5.0);
	print_compare("x^5", p_5, 	base_q_3, 0);

	@printf("\t quadrature degree 4:\n");
	print_compare("1", 	p_0, 	base_q_4, 2.0);
	print_compare("x", 	p_1, 	base_q_4, 0.0);
	print_compare("x^2", p_2, 	base_q_4, 2.0/3.0);
	print_compare("x^3", p_3, 	base_q_4, 0);
	print_compare("x^4", p_4, 	base_q_4, 2.0/5.0);
	print_compare("x^5", p_5, 	base_q_4, 0);
	print_compare("x^6", p_6, 	base_q_4, 2.0/7.0);
	print_compare("x^7", p_7, 	base_q_4, 0);

@printf("Interval [0,1]\n");
	@printf("\t quadrature degree 1:\n");
	print_compare("1", 	p_0, 	shifted_q_1, 1.0);
	print_compare("x", 	p_1, 	shifted_q_1, 0.5);

	@printf("\t quadrature degree 2:\n");
	print_compare("1", 	p_0, 	shifted_q_2, 1.0);
	print_compare("x", 	p_1, 	shifted_q_2, 0.5);
	print_compare("x^2", p_2, 	shifted_q_2, 1.0/3.0);
	print_compare("x^3", p_3, 	shifted_q_2, 1.0/4.0);

	@printf("\t quadrature degree 3:\n");
	print_compare("1", 	p_0, 	shifted_q_3, 1.00);
	print_compare("x", 	p_1, 	shifted_q_3, 0.5);
	print_compare("x^2", p_2, 	shifted_q_3, 1.0/3.0);
	print_compare("x^3", p_3, 	shifted_q_3, 1.0/4.0);
	print_compare("x^4", p_4, 	shifted_q_3, 1.0/5.0);
	print_compare("x^5", p_5, 	shifted_q_3, 1.0/6.0);

	@printf("\t quadrature degree 4:\n");
	print_compare("1", 	p_0, 	shifted_q_4, 1.0);
	print_compare("x", 	p_1, 	shifted_q_4, 0.5);
	print_compare("x^2", p_2, 	shifted_q_4, 1.0/3.0);
	print_compare("x^3", p_3, 	shifted_q_4, 1.0/4.0);
	print_compare("x^4", p_4, 	shifted_q_4, 1.0/5.0);
	print_compare("x^5", p_5, 	shifted_q_4, 1.0/6.0);
	print_compare("x^6", p_6, 	shifted_q_4, 1.0/7.0);
	print_compare("x^7", p_7, 	shifted_q_4, 1.0/8.0);


