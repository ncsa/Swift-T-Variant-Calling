# Desription

Testing on BW lacked organization, and after moving on for months, it was hard to come back and analyze all of the results.

Now, Jacob is trying to replicate the errors and successes of our previous experimentation, and he will keep the records here.

## Background on running on BW

On Blue waters, on each test, I use a starter sh script to make the call to Swift/T, which generates a qsub. But this qsub is missing some info that must be added manually. I kill the original submitted qsub and add the following

## 2 Node Test

**Note that even though the test is named 2 Node Test, it actually uses 3 nodes, as one extra is needed for Swift/T itself**

### 2 Node Test (2 samples per node)

/scratch/sciteam/jacobrh/purge_exempt/Swift_testing/newTesting/2Nodes

#### Background

This test appears to have always failed before. Going to replicate the error if possible.

In the 2Node.2samplesPerNode.sh file:

```
PPN=2 # this will make Swift/T run with 2 samples on each node
NODES=3
PROCS=$((PPN * NODES)) # This is the value directly passed to Swift/T's -n flag
```

This test is running on 4 samples

#### Result
