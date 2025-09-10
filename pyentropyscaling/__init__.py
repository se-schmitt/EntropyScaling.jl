from juliacall import Main as jl

jl.seval("using EntropyScaling")

for n in jl.names(jl.EntropyScaling):
    s = str(n)
    if not s[0] == '@':
        exec(f'{n} = jl.{n}')
