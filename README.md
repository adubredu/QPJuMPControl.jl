# QPJuMPControl.jl

A re-do of [Twan Koolen's](https://github.com/tkoolen) awesome [QPControl.jl](https://github.com/tkoolen/QPControl.jl) package but without the `Parametron.jl` dependency.

## Side Note
`Parametron.jl` was a cool package back in the day. It's `LazyExpressions` functionality is however broken beyond repair (At least not repairable by me. I tried. For days!). It seg faults when you take the product of `LazyExpressions`. And it doesn't seem to be actively maintained by Twan anymore :(.

As such, I decided to replicate `QPControl.jl` and substitute its `Parametron.jl` features with good old [JuMP.jl](https://jump.dev/) . This implementation may not be as fast as the original functional `QPControl.jl`, but it will at least work, especially with newer Julia versions.

## Acknowledgements
All credits go to [Twan Koolen](https://github.com/tkoolen) for his amazing work with the original [QPControl.jl](https://github.com/tkoolen/QPControl.jl) and on legged locomotion in general.
