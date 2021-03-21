# Sun Grid Engine (SGE)

* [Frequently used commands](#Frequently-used-commands)
    * [Submit jobs](#Submit-jobs)
    * [Parallel environment](#Parallel-environment)
    * [Queue configuration](#Queue-configuration)
    * [Queue status](#Queue-status)
* [Command not found error](#Command-not-found-error)

## Frequently used commands <a name="Frequently-used-commands"></a>

### Submit jobs <a name="Submit-jobs"></a>

To request a specific node:

```
qsub -l h=node_name example.sh
```

To request node A or node B:

```
qsub -l h='node_name_A|node_name_B' example.sh
```

To request 20 threads (cores) within a specific node using the parallel environment:

```
$ qsub -l h=node_name -pe pe_name 20 example.sh
```

To delete all jobs from a user:

```
qdel -u user_name
```

To delete a specific job:

```
qdel job_id
```

To print error message from a job:

```
qstat -j job_id | grep "error"
```

### Parallel environment <a name="Parallel-environment"></a>

To list all parallel environments:

```
qconf -spl
```

To print the configuration of a parallel environment:

```
qconf -sp pe_name
```

### Queue configuration <a name="Queue-configuration"></a>

To list all queues:

```
qconf -sql
```

To print the configuration of a queue:

```
qconf -sq queue_name
```

To list all administrative hosts (i.e. nodes for submitting jobs):

```
qconf -sh
```

To list all execution hosts (i.e. nodes for running jobs):

```
qconf -sel
```

### Queue status <a name="Queue-status"></a>

To print the status of all queues:

```
qstat -g c
```

To print the availability of all queues:

```
qstat -f
```

To print the availability of a queue:

```
qstat -f -q queue_name
```

To print all jobs currently occupying a queue:

```
qstat -u "*" | grep "queue_name"
```

To print the status of a host:

```
qhost -h host_name
```


## Command not found error <a name="Command-not-found-error"></a>

In some servers, even when a user submits a simple script to SGE, as simple as defining an environment variable, it returns an error complaining that command could not be found. However, when the user runs the same script locally or on a different cluster, it runs just fine. According to this Stack Overflow [post](https://stackoverflow.com/questions/17271931/sge-command-not-found-undefined-variable), the issue is "most likely the queues on your cluster are set to posix_compliant mode with a default shell of /bin/csh. The posix_compliant setting means your #! line is ignored. You can either change the queues to unix_behavior or specify the required shell using qsub's -S option":

```
#$ -S /bin/sh
``` 
