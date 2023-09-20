function Mdot = M_to_Mdot(M, q, q_dot)
    n_q = length(q);
    Mdot = zeros(n_q,n_q);
    for i=1:n_q
        Mdot = Mdot + simplify(diff(M,q(i))*q_dot(i),3);
    end
    Mdot = simplify(Mdot,3);
end