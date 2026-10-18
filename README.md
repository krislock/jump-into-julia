# Jump into Julia

Handout and code for "Jump into Julia," a Tech Tacos talk by Dr. Nathan Krislock
(NIU, 71-North studio, Tuesday, October 20, 2026).

The talk shows what makes Julia special: types and multiple dispatch. We build two new
kinds of numbers in a few lines each, and Julia's own code works with them.

1. Multiple dispatch in one minute
2. Integers mod p, with linear algebra mod p for free
3. Automatic differentiation with dual numbers, including the derivative of log det X
4. If time allows: why finite differences lose accuracy and dual numbers do not

## Files

- `index.html`: the handout. View it at <https://krislock.github.io/jump-into-julia/>.
- `jump_into_julia.jl`: all of the handout code in order, with checks at the end.

## Running the code

Install Julia from [julialang.org](https://julialang.org/downloads/). Then run:

```
julia jump_into_julia.jl
```

The script installs the [Primes.jl](https://github.com/JuliaMath/Primes.jl) package if you
do not have it. If all goes well, it ends with `All checks passed.`

During the talk we type the code into the Julia REPL instead. To install Primes.jl there,
press `]` at the `julia>` prompt and type `add Primes`.

## Credits

The dual-number code is adapted from an example by Joseph Miller on the Julia Zulip. The
Babylonian example follows Alan Edelman's "Automatic differentiation in 10 minutes with
Julia."

I used Claude, an AI assistant made by Anthropic, to help prepare this material. I tested
all of the code in Julia 1.13.2.
