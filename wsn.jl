M = 1:3;
J = 1:15;
c = [12.7 22.5 8.9 20.8 13.6 12.4 24.8 19.1 11.5 17.4 24.7 6.8 21.7 14.3 10.5; 19.1 24.8 24.4 23.6 16.1 20.6 15.0 9.5 7.9 11.3 22.6 8.0 21.5 14.7 23.2;  18.6 14.1 22.7 9.9 24.2 24.5 20.8 12.9 17.7 11.9 18.7 10.1 9.1 8.9 7.7; 13.1 16.2 16.8 16.7 9.0 16.9 17.9 12.1 17.5 22.0 19.9 14.6 18.2 19.6 24.2];
w = [61 70 57 82 51 74 98 64 86 80 69 79 60 76 78; 50 57 61 83 81 79 63 99 82 59 83 91 59 99 91;91 81 66 63 59 81 87 90 65 55 57 68 92 91 86; 62 79 73 60 75 66 68 99 69 60 56 100 67 68 54];
Q = [1020 1460 1530];

using JuMP, GLPK;

model = Model(GLPK.Optimizer)
@variable(model, x[m in M, j in J], Bin);
@constraint(model, cov[j in J], sum(x[m, j] for m in M) >= 1);
@constraint(model, knp[m in M], sum(w[m, j] * x[m, j] for j in J) <= Q[m]);
@objective(model, Min, sum(c[m, j] * x[m, j] for m in M, j in J));

optimize!(model);
objective_value(model)

using BlockDecomposition, Coluna;

coluna = optimizer_with_attributes(
    Coluna.Optimizer,
    "params" => Coluna.Params(
        solver = Coluna.Algorithm.TreeSearchAlgorithm() # default branch-cut-and-price
    ),
    "default_optimizer" => GLPK.Optimizer # GLPK for the master & the subproblems
);

@axis(M_axis, M);

M_axis[1]

typeof(M_axis[1])

model = BlockModel(coluna);

@variable(model, x[m in M_axis, j in J], Bin);
@constraint(model, cov[j in J], sum(x[m, j] for m in M_axis) >= 1);
@constraint(model, knp[m in M_axis], sum(w[m, j] * x[m, j] for j in J) <= Q[m]);
@objective(model, Min, sum(c[m, j] * x[m, j] for m in M_axis, j in J));

@dantzig_wolfe_decomposition(model, decomposition, M_axis)

decomposition

master = getmaster(decomposition)
subproblems = getsubproblems(decomposition)

specify!.(subproblems, lower_multiplicity = 0, upper_multiplicity = 1)
getsubproblems(decomposition)

optimize!(model)

value(x[1,3])
