Bash
****

Frequently used commands for Bash
=================================

Checksum
--------

Determine SHA-512 checksum:

.. code-block:: text

    $ shasum -a 256 example.txt

Determine MD5 checksum:

.. code-block:: text

    $ md5sum example.txt

Determine MD5 checksum in macOS:

.. code-block:: text

    $ md5 example.txt

List things
-----------

List files that don't contain a given string:

.. code-block:: text

    $ ls -1 | grep -v "pattern"

List one file per line:

.. code-block:: text

    $ ls -1 dir_name

List files in a space-separated view:

.. code-block:: text

    $ ls dir_name | tr '\n' ' '
    $ tr '\n' ' ' < list.txt

List all environment variables:

.. code-block:: text

    $ set

List all currently running processes:

.. code-block:: text

    $ ps aux

Here, the ``aux`` option means:

* ``a`` - show processes for all users
* ``u`` - show the process's owner
* ``x`` - show processes not attached to a terminal

List computer resources:

.. code-block:: text

    $ lscpu

List files with 0 bytes:

.. code-block:: text

    find . -size 0 -type f

Zipped files
------------

Create a .tar.gz file:

.. code-block:: text

    $ tar -czvf dir_name.tar.gz dir_name

Unzip a .tar.gz file:

.. code-block:: text

    $ tar -xf dir_name.tar.gz

Count things
------------

Count unique lines in a file:

.. code-block:: text

    $ sort example.txt | uniq -c | sort -bgr

Count files in a directory:

.. code-block:: text

    $ find dir_name | wc -l

Estimate size
-------------

Estimate storage size:

.. code-block:: text

    $ df -h

Estimate directory size:

.. code-block:: text

    $ du -sh dir_name
    $ du -sh dir_name --apparent-size

.. note::
    According to this `post <https://unix.stackexchange.com/questions/173947/du-s-apparent-size-vs-du-s>`__:
        The "apparent size" of a file is how much valid data is actually in the file. It is the actual amount of data that can be read from the file. Block-oriented devices can only store in terms of blocks, not bytes. As a result, the disk usage is always rounded up to the next highest block. A "block" in this case may not equate to a physical block on the storage device, either, depending on how the file system allocates space.

Estimate a server's memory usage

.. code-block:: text

    $ free -m

Comparison
----------

Find difference between two files:

.. code-block:: text

    $ grep -Fxvf file1 file2

Find difference between two directories:

.. code-block:: text

    $ diff -qr dir_name1 dir_name2

Check things
------------

Check if the system has a proxy server:

.. code-block:: text

    $ env | grep "proxy"

Check if the system is running a 32-bit or 64-bit OS:

.. code-block:: text

    $ uname -a

Check whether a file exists or not:

.. code-block:: text

    if test -f example.txt
    then
      echo "Found"
    else
      echo "Not found"
    fi

Check whether a variable exists or not:

.. code-block:: text

    if [ -z ${LC_ALL+x} ]
    then
      echo "LC_ALL is unset"
    else
      echo "LC_ALL is set to '$LC_ALL'"
    fi

Check Linux stuff for a Debian-based system:

.. code-block:: text

    cat /etc/*_version

Check Linux stuff for a Red Hat or CentOS-based system:

.. code-block:: text

    cat /etc/*-release

Module
------

List currently loaded modules:

.. code-block:: text

    $ module list

Load the latest version of a tool:

.. code-block:: text

    $ module load tool_name/latest

List available modules:

.. code-block:: text

    $ module avail

Load module or specify which dependencies have not been loaded:

.. code-block:: text

    $ module load modulefile

Internet
--------

Get IP address:

.. code-block:: text

    $ ip addr

Configure IP networking:

.. code-block:: text

    $ vi /etc/sysconfig/network-scripts/ifconfig-em1
    $ vi /etc/sysconfig/network-scripts/ifcfg-em1

File transfer
-------------

From local to server:

.. code-block:: text

    $ scp file.txt user_name@host_name:/path/to/destination
    $ scp file1.txt file2.txt user_name@host_name:/path/to/destination

From server to local:

.. code-block:: text

    $ scp user_name@host_name:/path/to/server/file.txt /path/to/destination
    $ scp -T user_name@host_name:"file1.txt file2.txt" "/path/to/destination"

Copy all files in a directory from server to local:

.. code-block:: text

    $ wget -r -c --no-parent --retry-connrefused /path/to/server/dir_name/

Here, the options mean:

    * ``-r`` - turn on recursive retrieving
    * ``-c`` - continue getting a partially-downloaded file
    * ``--no-parent`` - do not ever ascend to the parent directory when retrieving recursively
    * ``--retry-connrefused`` - consider "connection refused" a transient error and try again

Copy entire directory:

.. code-block:: text

  $ rsync -avzP source destination

Here, the ``-avzP`` option means:

* ``a`` - use archive mode
* ``v`` - be verbose
* ``z`` - compress file data during the transfer
* ``P`` - display progress and preserve partial files

* To only move files, and not directories, within the current directory to another:

    .. code-block:: text

        $ find . -maxdepth 1 -type f -exec mv {} dir_name \;

* To access a server and copy files:

    .. code-block:: text

        $ lftp sftp://user_id:user_pw.@host_name:port_number
        $ mirror -c target_dir destination_dir

awk
===

* To list columns by header name for a tab-delimited file:

    .. code-block:: text

        awk '
        NR==1 {
            for (i=1; i<=NF; i++) {
                f[$i] = i
            }
        }
        { print $(f["foo"]), $(f["baz"]) }
        ' example.txt

* To list columns by header name for a .csv file:

    .. code-block:: text

        awk -F "\"*,\"*" '
        NR==1 {
            for (i=1; i<=NF; i++) {
                f[$i] = i
            }
        }
        { print $(f["foo"]), $(f["baz"]) }
        ' example.csv

* To print lines that are both in file1.txt and file2.txt (intersection):

    .. code-block:: text

        $ awk 'NR == FNR{a[$0];next} $0 in a' file1.txt file2.txt

* To print lines that are only in file1.txt and not in file2.txt:

    .. code-block:: text

        $ awk 'NR == FNR{a[$0];next} !($0 in a)' file2.txt file1.txt

sed
===

* To search and replace a specific word from a line:

    .. code-block:: text

        $ echo "exampleword" | sed 's/word/new/g'


* To search and remove a specific word from a line:

    .. code-block:: text

        $ echo "exampleword" | sed 's/word//g'

vi and vim
==========

Frequently used commands for vi and vim
---------------------------------------

Search a pattern:

    * Press ``/``.
    * Type the search pattern.
    * Press ``Enter`` to perform the search.
    * Press ``n`` to find the next occurrence or ``N`` to find the previous occurrence.

Search and replace in the entire file:

.. code-block:: text

    :%s/foo/bar/g

Search and replace a pattern involving the ``/`` character:

.. code-block:: text

    :%s#/foo#/bar#g

Remove Windows line endings (type ``Ctrl-v`` and ``Ctrl-m`` to input ``^M``):

.. code-block:: text

    :%s/^M//g

Move the cursor to end of the file:

    Press the ``Esc`` key and then press the ``Shift`` and ``G`` keys together.

Bash configuration
==================

The ``.bashrc`` file is used to provide a place where you can set up variables, functions and aliases, define your (PS1) prompt and define other settings that you want to use every time you open a new terminal window. The following command will activate the configuration:

.. code-block:: text

    $ source .bashrc

There is also the ``.bash_profile`` file, which is executed for login shells, while ``.bashrc`` is executed for interactive non-login shells. When an installed program cannot be called from the command line, add the line ``export PATH=~/.local/bin:$PATH`` to the ``.bash_profile`` file.

System permission
=================

+-------+----------------+
| User  | rwx            |
+=======+================+
| Owner | ``-rwx------`` |
+-------+----------------+
| Group | ``----rwx---`` |
+-------+----------------+
| Other | ``-------rwx`` |
+-------+----------------+

+---+-------------------------+---------+--------+
| # | Permission              | rwx     | Binary |
+===+=========================+=========+========+
| 7 | read, write and execute | ``rwx`` | 111    |
+---+-------------------------+---------+--------+
| 6 | read and write          | ``rw-`` | 110    |
+---+-------------------------+---------+--------+
| 5 | read and execute        | ``r-x`` | 101    |
+---+-------------------------+---------+--------+
| 4 | read only               | ``r--`` | 100    |
+---+-------------------------+---------+--------+
| 3 | write and execute       | ``-wx`` | 011    |
+---+-------------------------+---------+--------+
| 2 | write only              | ``-w-`` | 010    |
+---+-------------------------+---------+--------+
| 1 | execute only            | ``--x`` | 001    |
+---+-------------------------+---------+--------+
| 0 | none                    | ``---`` | 000    |
+---+-------------------------+---------+--------+

For example, to give read, write, and execute permissions for everyone:

.. code-block:: text

    $ chmod 777 dir_name

Give permissions for all files inside the directory:

.. code-block:: text

    $ chmod 777 -R dir_name

Restore default permissions for the home directory:

.. code-block:: text

    $ chmod 755 home_dir

Make files read-only:

.. code-block:: text

    $ find . -type f -exec chmod a-w '{}' \;

Restore writing permission:

.. code-block:: text

    $ find . -type f -exec chmod +w '{}' \;

Sun Grid Engine (SGE)
=====================

`This <https://bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html>`__ website has good overview on how to submit a job using qsub.

Another good resource: https://info.hpc.sussex.ac.uk/hpc-guide/resource.html#memory-requests

Frequently used commands for SGE
--------------------------------

Submit jobs
^^^^^^^^^^^

Request more physical memory (default is 2 GB):

.. code-block:: text

    $ qsub -l m_mem_free=4G example.sh

Request more virtual memory (default i 2.5 GB):

.. code-block:: text

    $ qsub -l m_hvmem=4G example.sh

Request a specific node:

.. code-block:: text

    $ qsub -l h=node_name example.sh

Request node A or node B:

.. code-block:: text

    $ qsub -l h='node_A|node_B' example.sh

Request nodes in a specific queue:

.. code-block:: text

    $ qsub -q queue_name example.sh

Request 20 slots within a specific node using the parallel environment:

.. code-block:: text

    $ qsub -l h=node_name -pe pe_name 20 example.sh

Delete specific jobs:

.. code-block:: text

    $ qstat | grep "PATTERN" | awk '{print $1}' | xargs qdel

Delete all jobs from a user:

.. code-block:: text

    $ qdel -u user_name

Delete a specific job:

.. code-block:: text

    $ qdel job_id

Monitor resuource usage:

.. code-block:: text

    $ qstat -j job_id

Print error message from a running job:

.. code-block:: text

    $ qstat -j job_id | grep "error"

Print error message from a finished job:

.. code-block:: text

    $ qacct -j job_id | grep "error"

Queue configuration
^^^^^^^^^^^^^^^^^^^

To list all queues:

.. code-block:: console

    $ qconf -sql

To print the configuration of a queue:

.. code-block:: console

    $ qconf -sq queue_name

To list all administrative hosts (i.e. nodes for submitting jobs):

.. code-block:: console

    $ qconf -sh

To list all execution hosts (i.e. nodes for running jobs):

.. code-block:: console

    $ qconf -sel

To show the resources available for each node:

.. code-block:: console

    $ qhost -F
    $ qhost -F -h node_name

Queue status
^^^^^^^^^^^^

To print the status of all queues:

.. code-block:: console

    $ qstat -g c

To print the availability of all queues:

.. code-block:: console

    $ qstat -f

To print the availability of a queue:

.. code-block:: console

    $ qstat -f -q queue_name

To print all jobs currently occupying a queue:

.. code-block:: console

    $ qstat -u "*" | grep "queue_name"

To print the status of a host:

.. code-block:: console

    $ qhost -h host_name

Queue states
------------

Queue states or combinations of states can be:

+-------------+--------------------------------------------+
| Code        | Description                                |
+=============+============================================+
| a(larm)     | Alarm by reaching the limit of load.       |
+-------------+--------------------------------------------+
| A(larm)     | Alarm by reaching the limit of suspension. |
+-------------+--------------------------------------------+
| s(uspended) | Suspended.                                 |
+-------------+--------------------------------------------+
| S(uspended) | Suspended by subordination.                |
+-------------+--------------------------------------------+
| C           | Suspended by a timetable.                  |
+-------------+--------------------------------------------+
| d(isable)   | Execution host excluded from scheduling.   |
+-------------+--------------------------------------------+
| D(isable)   | Queue excluded from scheduling.            |
+-------------+--------------------------------------------+
| E(rror)     | Error state.                               |
+-------------+--------------------------------------------+
| U           | Unreachable                                |
+-------------+--------------------------------------------+
| au          | SGE is not running on the compute node.    |
+-------------+--------------------------------------------+

According to this `post <http://gridengine.org/pipermail/users/2011-September/001651.html>`__:

> The default is np_load_avg=1.75 with is more or less useless nowadays. Problem is, that also processes in state "D" (uninterruptible kernel task" which points to "waiting for disk" are there). So, a load higher than the installed cores times 1.75 can still be fine. [Originally it was the length of the process chain, i.e. number of process which are eligible to get some cpu cycles. As long as this number is lower than the number of installed cores, all processes are running at full speed (despite any set nice values), as there is noone to be nice to. Only with more processes than cores there is something to share)]

Command not found error
-----------------------

In some servers, even when a user submits a simple script to SGE, as simple as defining an environment variable, it returns an error complaining that command could not be found. However, when the user runs the same script locally or on a different cluster, it runs just fine. According to this Stack Overflow `post <https://stackoverflow.com/questions/17271931/sge-command-not-found-undefined-variable>`__, the issue is most likely the queues on your cluster are set to ``posix_compliant`` mode with a default shell of ``/bin/csh``. The ``posix_compliant`` setting means your ``#!`` line is ignored. You can either change the queues to ``unix_behavior`` or specify the required shell using the ``qsub -S`` option:

.. code-block:: console

    #$ -S /bin/sh

Usage monitoring
----------------

.. code-block:: console

    $ qstat -j job_id

+---------+----------------------------------------------+----------------+
| Field   | Description                                  | Examples       |
+=========+==============================================+================+
| cpu     | Total processing time                        | '00:00:08'     |
+---------+----------------------------------------------+----------------+
| mem     | Accumulated RAM of the job in Gbytes seconds | '80.89515 GBs' |
+---------+----------------------------------------------+----------------+
| io      | Total I/O usage                              | '0.07841'      |
+---------+----------------------------------------------+----------------+
| vmem    | Currently available RAM for the job          | '9.633G'       |
+---------+----------------------------------------------+----------------+
| maxvmem | Maximum RAM of the job when it was running   | '9.633G'       |
+---------+----------------------------------------------+----------------+

.. code-block:: console

    $ qhost -h node_name

+----------+-----------------------+-------------+
| Field    | Description           | Examples    |
+==========+=======================+=============+
| HOSTNAME | Node name             | 'cm463'     |
+----------+-----------------------+-------------+
| ARCH     | Architecture type     | 'linux-x64' |
+----------+-----------------------+-------------+
| NCPU     | Number of CPUs        | '56'        |
+----------+-----------------------+-------------+
| LOAD     | Current memory load   | '4.78'      |
+----------+-----------------------+-------------+
| MEMTOT   | Total memory          | '251.4G'    |
+----------+-----------------------+-------------+
| MEMUSE   | Currently used memory | '14.7G'     |
+----------+-----------------------+-------------+
| SWAPTO   | Swap space            | '8.0G'      |
+----------+-----------------------+-------------+
| SWAPUS   | Swap usage            | '2.3M'      |
+----------+-----------------------+-------------+

Parallel environment (PE)
-------------------------

List avaialble PEs:

.. code-block:: text

    $ qconf -spl

Print the configuration of a PE:

.. code-block:: text

    $ qconf -sp pe_name

Allocation rules
^^^^^^^^^^^^^^^^

Users can set ``allocation_rule`` to define how to assign slots to a job:

.. code-block:: text

    $ qconf -sp openmp
    pe_name            openmp
    slots              99999
    user_lists         NONE
    xuser_lists        NONE
    start_proc_args    NONE
    stop_proc_args     NONE
    allocation_rule    $fill_up
    control_slaves     FALSE
    job_is_first_task  TRUE
    urgency_slots      min
    accounting_summary FALSE

SGE provides three rules:

- ``$fill_up``: Use all of the job slots on a given host before moving to the next host.
- ``$round_robin``: Select one slot from each host in a round-robin fashion until all job slots are assigned. This setting can result in more than one job slot per host.
- ``$pe_slots``: Place all the job slots on a single machine. Grid Engine will only schedule such a job to a machine that can host the maximum number of slots requested by the job.

For example, if we want to assign 100 slots to two separate nodes ("bdcm02" and "bdcm04"):

.. code-block:: text

    $ qlogin -l h="bdcm02|bdcm04" -pe openmp 100

Ascp
====

`Ascp <https://www.ibm.com/docs/en/aci/3.9.2?topic=macos-ascp-transferring-from-command-line-ascp>`__ is a scriptable FASP transfer binary that enables you to transfer to and from Aspera transfer servers to which you have authentication credentials. Transfer settings are customizable and can include file manipulation on the source or destination, filtering of the source content, and client-side encryption-at-rest.

.. code-block:: text

    /mnt/garnet/Users/sbslee/.aspera/connect/bin/ascp -QT -l 100m -P33001 -k 1 -i /mnt/garnet/Users/sbslee/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/run/ERR323/ERR3239276/NA06985.final.cram .


Run commands in background
==========================

To run a command in the background, add the ampersand symbol (&) at the end
of the command:

.. code-block:: text

    $ command &

The job ID will be printed:

.. code-block:: text

    [1] 25177

Display the status of all background jobs:

.. code-block:: text

    $ jobs -l

Terminate a background job:

.. code-block:: text

    $ kill -9 25177

References:

   - `How to Run Linux Commands in Background <https://linuxize.com/post/how-to-run-linux-commands-in-background/#:~:text=A%20background%20process%20is%20a,the%20background%20processes%20is%20Linux.>`__

Connecting to a proxy server
============================

.. code-block:: text

    $ export http_proxy=http://129.156.243.243:3128
    $ export https_proxy=http://129.156.243.243:3128

Quick links:

- `Overriding system-repository Proxies by Using https_proxy and http_proxy <https://docs.oracle.com/cd/E36784_01/html/E37628/gmgas.html>`__

Miscellaneous
=============

A shebang is the character sequence consisting of the characters number sign and exclamation mark (#!) at the beginning of a script (``#!/bin/bash``).

Rename part of a filename:

.. code-block:: text

    $ for file in *.txt ; do mv $file ${file//ABC/XYZ} ; done

Remove file extension (e.g. ``.gz``):

.. code-block:: text

    $ mv -- "$file" "${file%%.gz}"

Access hard drives:

.. code-block:: text

    $ cd /
    $ cd Volumes
    $ cd ls

Extract lines repeated at least three times:

.. code-block:: text

    $ awk '++a[$0] == 3 { print $0 }' example.txt

Print every fifth line:

.. code-block:: text

    $ awk 'NR % 5 == 0' example.txt

Skip the first two lines of a file:

.. code-block:: text

    $ tail -n +3 example.txt

Concatenate a string to each line of the ``ls`` command output:

.. code-block:: text

    $ ls | xargs -i echo "Hello World {}"

Combine arrays as columns:

.. code-block:: text

    $ a=(A B C)
    $ b=(1 2 3)
    $ paste <(printf "%s\n" "${a[@]}") <(printf "%s\n" "${b[@]}")

Echo tab characters:

.. code-block:: text

    $ echo Hello$'\t'World

Read file names in the current directory into an array:

.. code-block:: text

    $ a=(*)

Redirect stdout and stderr:

.. code-block:: text

    $ some_command > out_file 2>error_file

Create a hard link or a symbolic link to an existing file or directory:

.. code-block:: text

    $ ln -s original_file new_file

Change group ownership:

.. code-block:: text

    $ chgrp -R group_name *

Get file basename:

.. code-block:: text

    $ basename /path/to/foo.txt
    foo.txt
    $ basename /path/to/foo.txt .txt
    foo

According to this `post <https://phoenixnap.com/kb/ulimit-linux-command#:~:text=ulimit%20is%20a%20built%2Din,users%20and%20system%20performance%20issues.>`__, ``ulimit`` is a built-in Linux shell command that allows viewing or limiting system resource amounts that individual users consume. Limiting resource usage is valuable in environments with multiple users and system performance issues.

Show the maximum number of file descriptors:

.. code-block:: text

    $ ulimit -n
    256

Reset the maximum number of file descriptors:

.. code-block:: text

    $ ulimit -n 300
    $ ulimit -n
    300
