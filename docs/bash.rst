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

Estimate a server's memory usage

.. code-block:: text

    $ free -m

Comparison
----------

Find difference between two directories:

.. code-block:: console

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

.. code-block:: console

    $ module list

Load the latest version of a tool:

.. code-block:: console

    $ module load tool_name/latest

List available modules:

.. code-block:: console

    $ module avail

Load module or specify which dependencies have not been loaded:

.. code-block:: console

    $ module load modulefile

File transfer
-------------

From local to server:

.. code-block:: console

    $ scp file.txt user_name@host_name:/path/to/destination
    $ scp file1.txt file2.txt user_name@host_name:/path/to/destination

From server to local:

.. code-block:: console

    $ scp user_name@host_name:/path/to/server/file.txt /path/to/destination
    $ scp -T user_name@host_name:"file1.txt file2.txt" "/path/to/destination"

Copy all files in a directory from server to local:

.. code-block:: console

    $ wget -r -c --no-parent --retry-connrefused /path/to/server/dir_name/

Here, the options mean:

    * ``-r`` - turn on recursive retrieving
    * ``-c`` - continue getting a partially-downloaded file
    * ``--no-parent`` - do not ever ascend to the parent directory when retrieving recursively
    * ``--retry-connrefused`` - consider "connection refused" a transient error and try again

Copy entire directory:

.. code-block:: console

  $ rsync -avzP source destination

Here, the ``-avzP`` option means:

* ``a`` - use archive mode
* ``v`` - be verbose
* ``z`` - compress file data during the transfer
* ``P`` - display progress and preserve partial files

* To only move files, and not directories, within the current directory to another:

    .. code-block:: console

        $ find . -maxdepth 1 -type f -exec mv {} dir_name \;

* To access a server and copy files:

    .. code-block:: console

        $ lftp sftp://user_id:user_pw.@host_name:port_number
        $ mirror -c target_dir destination_dir

Miscellaneous
-------------

To access hard drives:

.. code-block:: console

    $ cd /
    $ cd Volumes
    $ cd ls

To extract lines repeated at least three times:

.. code-block:: console

    $ awk '++a[$0] == 3 { print $0 }' example.txt

To print every fifth line:

.. code-block:: console

    $ awk 'NR % 5 == 0' example.txt

To skip the first two lines of a file:

.. code-block:: console

    $ tail -n +3 example.txt

To concatenate a string to each line of the ``ls`` command output:

.. code-block:: console

    $ ls | xargs -i echo "Hello World {}"

To combine arrays as columns:

.. code-block:: console

    $ a=(A B C)
    $ b=(1 2 3)
    $ paste <(printf "%s\n" "${a[@]}") <(printf "%s\n" "${b[@]}")

To echo tab characters:

.. code-block:: console

    $ echo Hello$'\t'World

To read file names in the current directory into an array:

.. code-block:: console

    $ a=(*)

To redirect stdout and stderr:

.. code-block:: console

    $ some_command > out_file 2>error_file

To create a hard link or a symbolic link to an existing file or directory:

.. code-block:: console

    $ ln -s original_file new_file

To change group ownership:

.. code-block:: console

    $ chgrp -R group_name *

awk
===

* To list columns by header name for a tab-delimited file:

    .. code-block:: console

        awk '
        NR==1 {
            for (i=1; i<=NF; i++) {
                f[$i] = i
            }
        }
        { print $(f["foo"]), $(f["baz"]) }
        ' example.txt

* To list columns by header name for a .csv file:

    .. code-block:: console

        awk -F "\"*,\"*" '
        NR==1 {
            for (i=1; i<=NF; i++) {
                f[$i] = i
            }
        }
        { print $(f["foo"]), $(f["baz"]) }
        ' example.csv

* To print lines that are both in file1.txt and file2.txt (intersection):

    .. code-block:: console

        $ awk 'NR == FNR{a[$0];next} $0 in a' file1.txt file2.txt

* To print lines that are only in file1.txt and not in file2.txt:

    .. code-block:: console

        $ awk 'NR == FNR{a[$0];next} !($0 in a)' file2.txt file1.txt

sed
===

* To search and replace a specific word from a line:

    .. code-block:: console

        $ echo "exampleword" | sed 's/word/new/g'


* To search and remove a specific word from a line:

    .. code-block:: console

        $ echo "exampleword" | sed 's/word//g'

vi and vim
==========

Frequently used commands for vi and vim
---------------------------------------

* To search a pattern:

    * Press ``/``.
    * Type the search pattern.
    * Press ``Enter`` to perform the search.
    * Press ``n`` to find the next occurrence or ``N`` to find the previous occurrence.

* To search and replace in the entire file:

    .. code-block:: console

        :%s/foo/bar/g

* To search and replace a pattern involving the ``/`` character:

    .. code-block:: console

        :%s#/foo#/bar#g

* To move the cursor to end of the file:

    Press the ``Esc`` key and then press the ``Shift`` and ``G`` keys together.

Looping through
===============

Print every line of a file:

.. code-block:: console

    for x in `cat example.txt`
    do
      echo "$x"
    done

Print only the second column:

.. code-block:: console

    for x in `awk '{print $2}' example.txt`
    do
      echo "$x"
    done

Print multiple columns:

.. code-block:: console

    while read -r first second
    do 
      echo "$first and $second"
    done < input.txt

Arrays
======

* To create an array:

    .. code-block:: console

        $ a=(1 2 3)
        $ a=(A B C)
        $ a=('A 1' 'B 2' 'C 3')

* To print an array:

    .. code-block:: console

        $ echo "${a[@]}"

* To print elements on separate lines:

    .. code-block:: console

        $ printf '%s\n' "${a[@]}"

* To loop through an array:

    .. code-block:: console

        $ cat example.sh
        a=(1 2 3)
        for x in ${a[@]}
        do
          echo $x
        done
        $ sh example.sh
        1
        2
        3

Bash configuration
==================

The ``.bashrc`` file is used to provide a place where you can set up variables, functions and aliases, define your (PS1) prompt and define other settings that you want to use every time you open a new terminal window. The following command will activate the configuration:

.. code-block:: console

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

.. code-block:: console

    $ chmod 777 dir_name

To give permissions for all files inside the directory:

.. code-block:: console

    $ chmod 777 -R dir_name

OpenSSH
=======

Frequently used commands for OpenSSH
------------------------------------

* To remove all keys belonging to a host name:

    .. code-block:: console

        $ ssh-keygen -R host_name

* To delete a select key from the authentication agent:

    .. code-block:: console

        $ ssh-add -d ~/.ssh/host_id_rsa.pub
        $ rm ~/.ssh/host_id_rsa
        $ rm ~/.ssh/host_id_rsa.pub

Creating a channel with password
--------------------------------

First, open your SSH configuration file:

.. code-block:: console

    $ vi ~/.ssh/config

Next, add the following:

.. code-block:: console

    Host host_id
        HostName host_name
        User user_name

Here, ``host_id`` is the nickname that will be used for the ``ssh`` command and ``host_name`` can be an IP address or an actual host name in the server. Lastly, ``user_name`` is your user ID for the server. After the configuration file is saved, you can access the server by (you still need to enter your password):

.. code-block:: console

    $ ssh host_id

Creating a channel without password
-----------------------------------

First, set up a channel with password as described above. Then, run the following:

.. code-block:: console

    $ ssh-keygen -t rsa -b 4096 -C "host_id"

Save the private key as ``host_id_rsa`` and the public key as ``host_id_rsa.pub``. Add the private key to the authentication agent:

.. code-block:: console

    $ ssh-add ~/.ssh/host_id_rsa

Check whether the addition was successful:

.. code-block:: console

    $ ssh-add -L

Add the public key to the server:

.. code-block:: console

    $ cat ~/.ssh/host_id_rsa.pub | ssh host_id 'cat >> ~/.ssh/authorized_keys'

Finally, update the configuration:

.. code-block:: console

    Host host_id
        HostName host_name
        User user_name
        IdentityFile ~/.ssh/host_id_rsa

Now, you shouldn't need to enter the password when logging in.

Channeling through multiple servers
-----------------------------------

Imagine the server you work on everyday (server C) can only be accessed through another server (server B). Inconveniently, server B can only be accessed through server A. So, your task is to set up a channel that looks like this: local > server A > server B > server C. To do this, you need to set up the SSH configuration as follows:

.. code-block:: console

    Host host_id_A
        HostName host_name_A
        User user_name_A
        IdentityFile ~/.ssh/host_id_A_rsa

    Host host_id_B
        HostName host_name_B
        User user_name_B
        ProxyCommand ssh host_id_A nc %h %p 2> /dev/null
        IdentityFile ~/.ssh/host_id_B_rsa

    Host host_id_C
        HostName host_name_C
        User user_name_C
        ProxyCommand ssh host_id_B nc %h %p 2> /dev/null
        IdentityFile ~/.ssh/host_id_C_rsa

You can now access server C directly by:

.. code-block:: console

    $ ssh host_id_C

Sun Grid Engine (SGE)
=====================

`This <https://bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html>`__ website has good overview on how to submit a job using qsub.

Frequently used commands for SGE
--------------------------------

Submit jobs
^^^^^^^^^^^

To request a specific node:

.. code-block:: console

    $ qsub -l h=node_name example.sh

To request node A or node B:

.. code-block:: console

    $ qsub -l h='node_A|node_B' example.sh

To request nodes in a specific queue:

.. code-block:: console

    $ qsub -q queue_name example.sh

To request 20 slots within a specific node using the parallel environment:

.. code-block:: console

    $ qsub -l h=node_name -pe pe_name 20 example.sh

To delete all jobs from a user:

.. code-block:: console

    $ qdel -u user_name

To delete a specific job:

.. code-block:: console

    $ qdel job_id

To monitor resuource usage:

.. code-block:: console

    $ qstat -j job_id

To print error message from a running job:

.. code-block:: console

    $ qstat -j job_id | grep "error"

To print error message from a finished job:

.. code-block:: console

    $ qacct -j job_id | grep "error"

Parallel environment
^^^^^^^^^^^^^^^^^^^^

To list all parallel environments:

.. code-block:: console

    $ qconf -spl

To print the configuration of a parallel environment:

.. code-block:: console

    $ qconf -sp pe_name

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
