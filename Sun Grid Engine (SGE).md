# Sun Grid Engine (SGE)

* [Frequently used commands](#Frequently-used-commands)
    * [Submit jobs](#Submit-jobs)

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





## FUCs for Sun Grid Engine <a name="FUCs-for-Sun-Grid-Engine"></a>

```
# -- Submitting jobs ---------------------------------------------------------

Delete existing jobs from a user (e.g. sbslee).
$ qdel -u <user_name>

Delete a specific job (e.g. 35345).
$ qdel <job_id>

Print error message from a job (e.g. 35345).
$ qstat -j <job_id> | grep "error"

# -- Parallel environment ----------------------------------------------------

List all parallel environments.
$ qconf -spl

Print the configuration of a parallel environment (e.g. make).
$ qconf -sp <pe_name>

# -- Queue configuration -----------------------------------------------------

List all queues.
$ qconf -sql

Print the configuration of a queue (e.g. biall.q).
$ qconf -sq <queue_name>

List all administrative hosts (i.e. nodes for submitting jobs).
$ qconf -sh

List all execution hosts (i.e. nodes for running jobs).
$ qconf -sel

# -- Queue status ------------------------------------------------------------

Print the status of all queues.
$ qstat -g c

Print the availability of all queues.
$ qstat -f

Print the availability of a queue (e.g. biall.q).
$ qstat -f -q <queue_name>

Print all jobs currently occupying a queue (e.g. biall.q).
$ qstat -u "*" | grep "<queue_name>"

Print the status of a host (e.g. grc145).
$ qhost -h <host_name>
```










## Command not found error <a name="Command-not-found-error"></a>

A Stack Overflow user mentioned:

> Most likely the queues on your cluster are set to posix_compliant mode with a default shell of /bin/csh.

Related posts:

* https://stackoverflow.com/questions/17271931/sge-command-not-found-undefined-variable
