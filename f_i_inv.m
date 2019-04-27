function [J_173_prime, t] = f_i_inv(J_15, J_173, J_N4)
	%options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
	%%%options = optimoptions(options,'SubproblemAlgorithm','factorization');
	lb = 0;
	ub = 1;
	J_0 = zeros(J_N4, 1);
	for J_i = 1 :4 : J_N4
		J_0(J_i) = 1;
	end
	%solve t from J_173 = t*J_15 + (1-t)*J_0;
	C = J_15 - J_0;
	d = J_173 - J_0;
	options = optimoptions('lsqlin','Algorithm','trust-region-reflective');
	t = lsqlin(C,d,[],[],[],[],lb,ub,[],options);

	%compute J_173_prime with t
	J_173_prime = t*J_15 + (1-t)*J_0;
	for J_i = 1 :4 : J_N4
		q = [J_173_prime(J_i:J_i+3)];
		q = q/norm(q);
		J_173_prime(J_i:J_i+3) = q;
	end
end