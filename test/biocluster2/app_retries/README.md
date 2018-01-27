This README serves to report back on efforts to enable robust  Swift/T app retries. The testing environment is Biocluster2, which is slurm based. Below are the directions for tests thus far, and lessons learnt:

# Slurm resources allocation

The main parameters to allocate resources for a Swift/T run are: 
`PPN`: Processes-per-node
`PROCS`: Number of MPI processes
`NODES`:

These are defined in the script [run-init.zsh  ](https://github.com/swift-lang/swift-t/blob/master/turbine/code/scripts/submit/run-init.zsh), before propogating to the respective slurm submission scripts [turbine-slurm-run.zsh](https://github.com/swift-lang/swift-t/blob/master/turbine/code/scripts/submit/slurm/turbine-slurm-run.zsh) and [turbine-slurm.sh.m4](https://github.com/swift-lang/swift-t/blob/master/turbine/code/scripts/submit/slurm/turbine-slurm.sh.m4i#L46)

The key peace of code that govrens these variables is below:

```
# Round NODE up for extra processes
export NODES=$(( PROCS/PPN ))
(( PROCS % PPN )) && (( NODES++ )) || true
export TURBINE_WORKERS=$(( PROCS - ADLB_SERVERS ))
declare NODES PROCS PPN 
```
But, tests- as per below, don't seem to agree:

## Resources folder:

This is a warm up exercise, just to see how the Swift/T variables `PPN`, `PROCS` and `NODES` affect the actual resources reserved for a typical Swift/T run. Below is a summary table of findings, with actual files from eact test being within the relevant directory within the [Resources directory here](https://github.com/ncsa/Swift-T-Variant-Calling/tree/master/test/biocluster2/app_retries/Resources), including the email_report recieved from the scheduler upon the start of a given job.


|Folder_name|PPN | PROCS| NODES|Resulting_Resources_Allocated|
|-------|----|------|------|------|
|retry_2|3|12|5| 15 MPI ranks, on 5 nodes|
|retry_1|3|12|-| 6 MPI ranks, on 2 nodes\* |
|retry_3|3|-|5|  15 MPI ranks, on 5 nodes|
|retry_4|-|12|5| 5 MPI ranks, on 5 nodes|
|retry_5|-|-|5| 5 MPI ranks, on 5 nodes|
|retry_6|-|5|-| 2 MPI ranks, on 2 nodes\* |
|retry_7|5|-|-| 6 MPI ranks, on 2 nodes+|
\ The default values are `PPN=1`
These tests make me question what the effect of the `PROCS` parameter (defined as the total number of MPI ranks) is. 

# Important note: 
Retries applies to leaf function failures, not to programatic errors. For example, exceeding array limits would cause the program to fail, but no retry apply in this case (this is obviious, but it helps to write out)

# Location directive 

Probably, the main lesson learnt here is that a location object is a data structure with 3 components: `rank`, `strictness`, and `accuracy`; where:
`rank`: is the mpi-rank
`strictness`: is the location specification constratint, which can be either `HARD` or `SOFT`
`Accuracy`: is the choice to run a given task on either a specific `RANK` or `NODE`

- Random number generation in Swift/T is not really random, it always starts from a certain seed (the seed itself being dependent on the MPI rank on which the process is to be evaluated. Control of the seed is achieved via the turbine variable `TURBINE_SRAND`. This also applies to how a random worker maybe selected to carry out a certain function.

- By default, functions are local **local DOES NOT mean they are run in the ADLB server node, or the like- it merely seems to mean that their location can not be controlled. See `t6_function_app.swift` for an example**. To change this and be able to *send a function into any location, it must has the directive: `@dispatch=WORKER`. Otherwise, stc would produce the compilation error: "Tried to call non-targetable function... "* To the contrary, `app functions`, aka leaf functions, can be dispatched without explicitly using the `dispatch directive`

- To check if 2 variables are identical, the `assert` library can be used. In additon to the vanila version, more variants exist with the general body: `assertEqual(var1, var2, test_name_or_message)`.

- App retries do not work in the following cases: 
1) when an `assert` statement fails
2) when there is a coding error (for example, a loop going beyond its boundries)  

- I'm not sure about **the use of the `type` keyword.** In the Swift/T user guide, it should be used to define new structures as in:

```
type person {
	string name,
	int age }
```
There is also the `typedef` keyword which is used to define new types as in:

```
typedef newtype oldtype;
```

In the test file: `t5_types.swift`, which is copied from [the test file](https://github.com/swift-lang/swift-t/blob/master/stc/tests/665-subtype.swift); it is used as a function declaration. **Kindly elaborate**

- As expected, Swift/T is highly typed language. Probably, only int-->float conversion can happen.

- The function `hostmapOneRank` and its alias `hostmapOne` are quite strict about extracting the hostname from the hostmap. The call to this function fails when done repeatedly enough. See the `688-output.txt` for an example error (resulting from running `t8_based_on_688_hard_location.swift`). The easiest way to correct this is to relax the `HARD` strictness of finding the hostmap, via using the `hostmapOneWorkerRank` function.

- There is also the `par` annotation, which has been used in the test: [689](https://github.com/swift-lang/swift-t/blob/7cdece77a97683f8338f29521845a4bbbf6bc635/stc/tests/689-parloc-1.swift). I modified the location to a proper form, and used another function definition, just to get the tast/purpose of this test, but still- it is failing at the compilation stage:

```
ERROR STC internal error: please report this
exm.stc.common.exceptions.STCRuntimeError: Don't support wrapping parallel functions yet
	at exm.stc.frontend.WrapperGen.generateWrapper(WrapperGen.java:198)
	at exm.stc.frontend.WrapperGen.generateWrapper(WrapperGen.java:187)
	at exm.stc.frontend.FunctionCallEvaluator.backendCallWrapped(FunctionCallEvaluator.java:556)
	at exm.stc.frontend.FunctionCallEvaluator.backendFunctionCall(FunctionCallEvaluator.java:534)
	at exm.stc.frontend.FunctionCallEvaluator.evalFunctionCallInner(FunctionCallEvaluator.java:148)
	at exm.stc.frontend.FunctionCallEvaluator.evalFunctionCall(FunctionCallEvaluator.java:107)
	at exm.stc.frontend.ExprWalker.evalToVars(ExprWalker.java:99)
	at exm.stc.frontend.ASTWalker.exprStatement(ASTWalker.java:1583)
	at exm.stc.frontend.ASTWalker.walkStatement(ASTWalker.java:654)
	at exm.stc.frontend.ASTWalker.walkTopLevelCompileStatements(ASTWalker.java:389)
	at exm.stc.frontend.ASTWalker.walkTopLevel(ASTWalker.java:283)
	at exm.stc.frontend.ASTWalker.walkFile(ASTWalker.java:256)
	at exm.stc.frontend.ASTWalker.compileTopLevel(ASTWalker.java:216)
	at exm.stc.frontend.ASTWalker.walk(ASTWalker.java:180)
	at exm.stc.ui.STCompiler.compileOnce(STCompiler.java:109)
	at exm.stc.ui.STCompiler.compile(STCompiler.java:71)
	at exm.stc.ui.Main.main(Main.java:103)
```
In the above example, I also note that I couldn't dispatch a normal swift/t function as shown. Is there a rule here? Ans. yes:
**You can only dispatch: 1)app function, 2) tcl functions**
**For pure Swift/t functions, you can apply `par` annotation**




