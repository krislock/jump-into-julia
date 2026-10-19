# All code from the handout, in order, plus checks.
# Run with:  julia jump_into_julia.jl

describe(x)          = "something"
describe(x::Number)  = "a number"
describe(x::Integer) = "an integer"
describe(x::String)  = "a string"

describe(1), describe(1.5), describe("hi"), describe([1, 2])
# ("an integer", "a number", "a string", "something")

methods(+)           # + itself has around 200 methods

import Pkg; Pkg.add("Primes")              # one time only; the first add can take a minute
using Primes              # load it in each new session
isprime(7), isprime(91)   # (true, false)   91 = 7 × 13

struct Mod{p} <: Integer
    k::Int
    function Mod{p}(k) where p
        @assert isprime(p)
        new(mod(k, p))
    end
end

import Base: +, -, *, /, inv, show, promote_rule, oneunit, abs, rand

show(io::IO, a::Mod{p}) where p =
    print(io, get(io, :compact, false) ? a.k : "$(a.k) (mod $p)")

+(a::Mod{p}, b::Mod{p}) where p = Mod{p}(a.k + b.k)
-(a::Mod{p}) where p = Mod{p}(-a.k)
-(a::Mod{p}, b::Mod{p}) where p = a + -b
*(a::Mod{p}, b::Mod{p}) where p = Mod{p}(a.k * b.k)
inv(a::Mod{p}) where p = Mod{p}(invmod(a.k, p))
/(a::Mod{p}, b::Mod{p}) where p = a * inv(b)

promote_rule(::Type{Mod{p}}, ::Type{Int}) where p = Mod{p}   # Mod + Int gives Mod
oneunit(::Type{Mod{p}}) where p = Mod{p}(1)
abs(a::Mod{p}) where p = a.k        # used by \ to choose pivots

rand(::Type{Mod{p}}, d::Integer, dims::Integer...) where p =
    Mod{p}.(rand(0:p-1, d, dims...))

a = Mod{7}(3)
b = Mod{7}(5)
a + b        # 1 (mod 7)   3 + 5 = 8, and dividing 8 by 7 leaves a remainder of 1
a * b        # 1 (mod 7)   3 × 5 = 15, and dividing 15 by 7 leaves a remainder of 1
inv(a)       # 5 (mod 7)   since 3 × 5 ≡ 1
a / b        # 2 (mod 7)   since 5 × 2 = 10 ≡ 3
a + 10       # 6 (mod 7)   3 + 10 = 13 ≡ 6; Julia turns 10 into a Mod{7} for us
a^100        # 4 (mod 7)   we never wrote ^ ; Julia's generic power code works

A = [2 1; 1 3]
b = [7, 11]
x = A \ b    # [2.0, 3.0]
A * x        # [7.0, 11.0], which is b
inv(A)       # the inverse of A
inv(A) * b   # also [2.0, 3.0]

A = Mod{7}.([2 1; 1 3])
b = Mod{7}.([1, 2])
x = A \ b        # the solution, mod 7
A * x            # [1, 2], which is b

inv(A)           # the inverse matrix, mod 7
inv(A) * b       # the same solution (but in practice, use A \ b)

struct Dual{T<:Real} <: Number
    p::T     # value
    d::T     # derivative
end
Dual(p, d) = Dual(promote(p, d)...)

import Base: +, -, *, /, inv, sin, cos, exp, show, promote_rule, convert

+(x::Dual, y::Dual) = Dual(x.p + y.p, x.d + y.d)              # sum rule
-(y::Dual)          = Dual(-y.p, -y.d)
-(x::Dual, y::Dual) = x + -y
*(x::Dual, y::Dual) = Dual(x.p * y.p, x.d * y.p + x.p * y.d)  # product rule
inv(y::Dual)        = Dual(1 / y.p, -y.d / y.p^2)
/(x::Dual, y::Dual) = x * inv(y)                               # quotient rule

sin(y::Dual) = Dual(sin(y.p),  cos(y.p) * y.d)                 # chain rule
cos(y::Dual) = Dual(cos(y.p), -sin(y.p) * y.d)
exp(y::Dual) = Dual(exp(y.p),  exp(y.p) * y.d)

promote_rule(::Type{Dual{T}}, ::Type{S}) where {T<:Real, S<:Real} =
    Dual{promote_type(T, S)}
convert(::Type{Dual{T}}, x::Real) where T<:Real = Dual{T}(x, zero(T))
convert(::Type{Dual{T}}, x::Dual) where T<:Real = Dual{T}(x.p, x.d)

show(io::IO, x::Dual) = print(io, x.p, " + ", x.d, "*ε")

ε = Dual(0, 1)

(3 + ε)^2                    # 9 + 6*ε
ε * ε                        # 0 + 0*ε
f(x)  = sin(exp(x + 1))
fp(x) = cos(exp(x + 1)) * exp(x + 1)    # derivative by hand

f(4.0), fp(4.0)              # (-0.6876914116945115, -107.74841574308525)
f(4.0 + ε)                   # -0.6876914116945115 + -107.74841574308525*ε

function babylonian(x; N = 10)
    t = (1 + x) / 2
    for i in 2:N
        t = (t + x / t) / 2
    end
    return t
end

babylonian(2.0 + ε)          # 1.41421356237... + 0.35355339059...*ε
sqrt(2), 1 / (2sqrt(2))      # check: √2 and d/dx √x at x = 2

using LinearAlgebra

Base.log(x::Dual) = Dual(log(x.p), x.d / x.p)
Base.abs(x::Dual) = abs(x.p)    # LU pivots on the value part
Base.oneunit(::Type{Dual{T}}) where T = Dual{T}(oneunit(T), zero(T))

X = [4.0 1 2; 0 3 1; 1 0 5]
D = Dual.(X, 0.0)        # X with zero ε parts
D[1, 2] += ε             # add ε to the (1,2) entry
log(det(D))              # 4.00733... + 0.0181818...*ε
inv(X)[2, 1]             # 0.0181818...  (it is 1/55)

t = 1e-6
E = zeros(3, 3); E[1, 2] = 1
(log(det(X + t*E)) - log(det(X))) / t    # about 0.0181818, close but not exact

fd(f, x, h) = (f(x + h) - f(x)) / h   # f is sin(exp(x + 1)) from section 3

f(4.0 + ε)               # dual numbers: -107.74841574308525 in the ε part
fp(4.0)                  # hand-written derivative: the same number

for k in 1:16
    h = 10.0^(-k)
    err = abs(fd(f, 4.0, h) - fp(4.0))
    println("h = 1e-", k, "    error = ", err)
end
# the error is smallest near h = 1e-9 (about 4e-6), then grows again;
# at h = 1e-16, 4.0 + h == 4.0 in floating point, so the error is the whole
# derivative, about 108

# ---------- checks (added for testing; not on the handout) ----------
println("\nRunning checks on Julia ", VERSION)
@assert describe(1) == "an integer" && describe([1, 2]) == "something"
@assert isprime(7) && !isprime(91)
a, b = Mod{7}(3), Mod{7}(5)
@assert a + b == Mod{7}(1)
@assert a * b == Mod{7}(1)
@assert inv(a) == Mod{7}(5)
@assert a / b == Mod{7}(2)
@assert a + 10 == Mod{7}(6)
@assert a^100 == Mod{7}(4)
@assert [2 1; 1 3] \ [7, 11] == [2.0, 3.0]
@assert inv([2 1; 1 3]) * [7, 11] ≈ [2.0, 3.0]
@assert inv([2 1; 1 3]) ≈ [0.6 -0.2; -0.2 0.4]
A = Mod{7}.([2 1; 1 3]); bb = Mod{7}.([1, 2])
@assert A \ bb == Mod{7}.([3, 2])
@assert inv(A) == Mod{7}.([2 4; 4 6])
@assert inv(A) * bb == Mod{7}.([3, 2])
@assert ε * ε == Dual(0, 0)
@assert (3 + ε)^2 == Dual(9, 6)
@assert Mod{7}(100) == Mod{7}(2)
@assert f(4.0 + ε).d ≈ fp(4.0)
@assert babylonian(2.0 + ε).p ≈ sqrt(2)
@assert babylonian(2.0 + ε).d ≈ 1 / (2sqrt(2))
L = log(det(D))
display(L)
@assert L.p ≈ log(det(X))
@assert L.d ≈ inv(X)[2, 1] ≈ 1/55
@assert abs((log(det(X + t*E)) - log(det(X))) / t - 1/55) < 1e-6
@assert f(4.0 + ε).d ≈ fp(4.0)
@assert abs(fd(f, 4.0, 1e-16) - fp(4.0)) ≈ abs(fp(4.0))
@assert abs(fd(f, 4.0, 1e-9) - fp(4.0)) < 1e-4
println("All checks passed.")
