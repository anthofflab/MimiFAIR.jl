using Mimi
using NLsolve

@defcomp fair begin

    # Carbon cycle
    a = Parameter(index=[4])
    τ = Parameter(index=[4])
    E = Parameter(index=[time])
    α = Variable(index=[time])
    R = Variable(index=[time,4])
    C0 = Parameter()
    C = Variable(index=[time])
    Cacc = Variable(index=[time])
    r0 = Parameter()
    rC = Parameter()
    rT = Parameter()

    # RF forcing
    F2x = Parameter()
    Fext = Parameter(index=[time])
    F = Variable(index=[time])

    # Temp
    T = Variable(index=[time])
    Tj = Variable(index=[time,2])
    c = Parameter(index=[2])
    d = Parameter(index=[2])
end

function run_timestep(s::fair, t::Int)
    p = s.Parameters
    v = s.Variables

    if t==1
        # TODO Is this correct?
        for i=1:4
            v.R[t,i] = p.a[i]*p.E[t]
        end
        v.C[t] = p.C0
        v.Cacc[t] = sum(p.E[1:t]) - (v.C[t]-p.C0)
        v.F[t] = p.F2x/log(2)*log(v.C[t]/p.C0) + p.Fext[t]
        for j=1:2
            v.Tj[t,j] = 0.
        end
        v.T[t] = sum(v.Tj[t,:])
        v.α[1]=1.
    else
        # Solve for α
        iIRFT100 = p.r0 + p.rC * v.Cacc[t-1] + p.rT * v.T[t-1]

        function f!(x, fvec)
            α = x[1]
            res = -iIRFT100
            for i=1:4
                temp = α*p.a[i]*p.τ[i] * (1 - exp(-100/α/p.τ[i]))
                res += temp
            end
            fvec[1] = res
        end

        v.α[t] = nlsolve(f!, [0.01], autodiff=true, method=:newton).zero[1]

        # Compute concentrations

        for i=1:4
            v.R[t,i] = v.R[t-1,i] + p.a[i]*p.E[t] - v.R[t-1,i]/v.α[t]/p.τ[i]
        end
        v.C[t] = p.C0 + sum(v.R[t,:])
        # TODO Change to a recursive equation, should be a bit faster
        v.Cacc[t] = sum(p.E[1:t]) - (v.C[t]-p.C0)

        # Compute forcing

        v.F[t] = p.F2x/log(2)*log(v.C[t]/p.C0) + p.Fext[t]

        # Compute temperature

        for j=1:2
            v.Tj[t,j] = v.Tj[t-1,j] + (p.c[j]*v.F[t]-v.T[t-1])/p.d[j]
        end
        v.T[t] = sum(v.Tj[t,:])
    end
end
