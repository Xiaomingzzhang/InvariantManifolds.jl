f1(x, p, t) = SA[x[2], p[1]*x[1]]

f2(x, p, t) = SA[x[2], -p[2]*x[1]]

hyper(x, p, t) = x[1] - p[3]

dom1(x, p, t) = p[3] < x[1]

dom2(x, p, t) = x[1] > p[3]

initialconditions=(SA[1.0, 5.0], SA[1.0, 1.0, 2.0], 0.0)
@test PiecewiseV((f1, f2), (dom1, dom2), (hyper,)).n ==0
@test PiecewiseV((f1, f2), (dom1, dom2), (hyper,),1)(initialconditions...)==f1(initialconditions...)
@test PiecewiseV((f1, f2), (dom1, dom2), (hyper,),2)(initialconditions...)==f2(initialconditions...)
@test BilliardV(f1, (dom1, dom2), (hyper,hyper))(initialconditions...)==f1(initialconditions...)

