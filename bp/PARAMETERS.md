# List of parameters

This page lists all the parameters specific to the PSLP solver.


## Heuristics

**frequency for calling primal heuristic <pslp_init> (-1: never, 0: only at depth freqofs)**

type: int, range: [-1,65534], default: 1

```
heuristics/pslp_init/freq = 1
```


## Reading

**use adjacent conflicts when unspecified?**

type: bool, range: {TRUE,FALSE}, default: TRUE

```
reading/jsonreader/useadjacentconflicts = TRUE
```


## All pricers

These parameters are duplicated for all pricers. Replace `{name}` by the name of the pricer.

**should pricer {name} be activated?**

type: bool, range: {TRUE,FALSE}, default: TRUE

```
pricers/{name}/activate = TRUE
```

**priority of pricer {name}**

type: int, range: [-536870912,536870911], default: 0

```
pricers/{name}/priority = 0
```

**file name to export pricer instances**

type: string, default: ""

```
pricers/{name}/probfilename = ""
```


## CPLEX pricers

Replace `{name}` by the name of the pricer.

**global thread count (0: auto)**

type: int, range: [0,2100000000], default: 0

```
pricers/{name}/threads = 0
```

**optimizer time limit in seconds**

type: real, range: [0,1e+75], default: 1e+75

```
pricers/{name}/timelimit = 1e+75
```

**memory available for working storage**

type: real, range: [0,1e+75], default: 2048

```
pricers/{name}/workmem = 2048
```

**node storage file switch (0: no node file, 1: in memory and compressed, 2: on disk, 3: on disk and compressed)**

type: int, range: [0,3], default: 1

```
pricers/{name}/mip/strategy/file = 1
```

**tree memory limit**

type: real, range: [0,1e+75], default: 1e+75

```
pricers/{name}/mip/limits/treememory = 1e+75
```

**MIP integer solution limit**

type: longint, range: [1,9223372036800000000], default: 9223372036800000000

```
pricers/{name}/mip/limits/solutions = 9223372036800000000
```

**upper objective value limit**

type: real, range: [-1e+75,1e+75], default: 1e+75

```
pricers/{name}/simplex/limits/upperobj = 1e+75
```

**display pricer output log?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/{name}/display/out = FALSE
```

**display pricer warnings?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/{name}/display/warning = FALSE
```

**file name to export lp files**

type: string, default: ""

```
pricers/{name}/lpfilename = ""
```

**can the pricer exceed the global time limit?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/{name}/extratime = FALSE
```

**add initial 2-cycle elimination constraints?**

type: bool, range: {TRUE,FALSE}, default: TRUE

```
pricers/cplex_lazy/add2cyclecons = TRUE
```

**add initial cycle elimination constraints found by DFS?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/cplex_lazy/adddfscyclecons = FALSE
```

**maximal number of cuts added per round (0: unlimited)**

type: int, range: [0,2147483647], default: 0

```
pricers/cplex_lazy/maxcuts = 0
```


## MIP pricers

Replace `{name}` by the name of the pricer.

**file name to export pricer instances**

type: string, default: ""

```
pricers/{name}/probfilename = ""
```

**maximal time in seconds to run**

type: real, range: [0,1e+20], default: 1e+20

```
pricers/{name}/limits/time = 1e+20
```

**maximal memory usage in MB**

type: real, range: [0,8796093022207], default: 8796093022207

```
pricers/{name}/limits/memory = 8796093022207
```

**verbosity level of output**

type: int, range: [0,5], default: 0

```
pricers/{name}/display/verblevel = 0
```

**should the CTRL-C interrupt be caught by the pricer?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/{name}/misc/catchctrlc = FALSE
```

**can the pricer exceed the global time limit?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/{name}/extratime = FALSE
```

**file name to export lp files**

type: string, default: ""

```
pricers/{name}/lpfilename = ""
```

**add initial 2-cycle elimination constraints?**

type: bool, range: {TRUE,FALSE}, default: TRUE

```
pricers/mip_lazy/add2cyclecons = TRUE
```

**add initial cycle elimination constraints found by DFS?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/mip_lazy/adddfscyclecons = FALSE
```

**maximal number of cuts added per round (0: unlimited)**

type: int, range: [0,2147483647], default: 0

```
pricers/mip_lazy/maxcuts = 0
```


## Tree search pricer

**node selection rule (bfs, brfs, cbfs, dfs)**

type: string, default: "dfs"

```
pricers/tree/nodeselection = "dfs"
```

**maximal time in seconds to run**

type: real, range: [0,1e+20], default: 1e+20

```
pricers/tree/timelimit = 1e+20
```

**can the pricer exceed the global time limit?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/tree/extratime = FALSE
```


## Multi pricer

Replace `{name}` by the name of the subpricer.

**use {name} pricer?**

type: bool, range: {TRUE,FALSE}, default: TRUE

```
pricers/multi/{name}/use = TRUE
```

**priority of {name} pricer**

type: int, range: [-536870912,536870911], default: 1000000

```
pricers/multi/{name}/priority = 1000000
```

**should the pricer add variables?**

type: bool, range: {TRUE,FALSE}, default: FALSE

```
pricers/multi/addvars = FALSE
```