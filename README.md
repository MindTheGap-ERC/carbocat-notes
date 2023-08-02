# Project README

This repository contains a series of notes and code implementations for the Mind The Gap project. As products materialize from this, they should be moved to their own repository.

## Julia

There is a `MindTheGap` Julia project in here. The best way to interact with code is by starting the Julia REPL.

```shell
julia --project=.
```

or

```shell
make repl
```

Then inside the REPL, first load `Revise`, then `MindTheGap`. `Revise` will keep your running image up-to-date with latests changes in the Julia code:

```julia
using Revise
using MindTheGap
```

You may also want to experiment with the code through Jupyter, see [IJulia](https://github.com/JuliaLang/IJulia.jl) for instructions.

## Entangled

All code in this repository is written using Entangled to synchronize source code with Markdown. To install Entangled run the following:

```shell
pip --user install entangled-cli[rich]
```

Note, that the `[rich]` part is optional, it enables pretty coloured output on the terminal. When editing source files or Markdown, please remember to have the Entangled daemon running:

```shell
entangled watch
```

## Project layout

| directory | contents |
| --- | --- |
| `./docs` | Markdown sources as well as generated figures. |
| `./site` | Generated HTML pages. |
| `./data` | External data, some of it generated. |
| `./src`  | Julia source code, generated from the Markdown in `./docs`. |
| `./test` | Unit tests. |
| `./.entangled` | Additional files needed to generate HTML pages. |

## Building pages

To build the HTML pages, you need `pandoc` installed, see [pandoc.org](https://pandoc.org/). To build the pages run

```shell
make site
```
