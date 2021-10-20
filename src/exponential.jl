mutable struct ExplicitExponentialRungeKuttaSolver{stepsize_T} <: AbstractRungeKuttaSolver
    stepsize::stepsize_T
end

ExplicitExponentialRungeKuttaSolver(h::Real) = ExplicitExponentialRungeKuttaSolver(StepSize(h))
@doc (@doc ExplicitExponentialRungeKuttaSolver) EERK(args...; kwargs...) = ExplicitExponentialRungeKuttaSolver(args...; kwargs...)

############################################################################################
#                                         PRINTING                                         #
############################################################################################

"""
    show(io::IO, solver::ExplicitExponentialRungeKuttaSolver)

prints a full description of `solver` and its contents to a stream `io`.
"""
Base.show(io::IO, solver::ExplicitExponentialRungeKuttaSolver) = NSDEBase._show(io, solver)

"""
    summary(io::IO, solver::ExplicitExponentialRungeKuttaSolver)

prints a brief description of `solver` to a stream `io`.
"""
Base.summary(io::IO, solver::ExplicitExponentialRungeKuttaSolver) = NSDEBase._summary(io, solver)


############################################################################################
#                                          METHODS                                         #
############################################################################################

function (solver::ExplicitExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem)
    solve(problem, solver)
end
