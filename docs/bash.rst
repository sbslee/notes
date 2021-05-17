Bash
****

Arrays
======

To create an array:

.. code-block:: console

    $ a=(1 2 3)
    $ a=(A B C)
    $ a=('A 1' 'B 2' 'C 3')

To print an array:

.. code-block:: console

    $ echo "${a[@]}"

To print elements on separate lines:

.. code-block:: console

    $ printf '%s\n' "${a[@]}"

To loop through an array:

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

To remove all keys belonging to a host name:

.. code-block:: console

    $ ssh-keygen -R host_name

To delete a select key from the authentication agent:

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

.. parsed-literal::

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

.. parsed-literal::

    Host host_id
        HostName host_name
        User user_name
        IdentityFile ~/.ssh/host_id_rsa

Now, you shouldn't need to enter the password when logging in.

Channeling through multiple servers
-----------------------------------

Imagine the server you work on everyday (server C) can only be accessed through another server (server B). Inconveniently, server B can only be accessed through server A. So, your task is to set up a channel that looks like this: local > server A > server B > server C. To do this, you need to set up the SSH configuration as follows:

.. parsed-literal::

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

Frequently used commands for SGE
--------------------------------

Submit jobs
^^^^^^^^^^^

To request a specific node:

.. code-block:: console

    $ qsub -l h=node_name example.sh

To request node A or node B:

.. code-block:: console

    $ qsub -l h='node_name_A|node_name_B' example.sh

To request 20 threads (cores) within a specific node using the parallel environment:

.. code-block:: console

    $ qsub -l h=node_name -pe pe_name 20 example.sh

To delete all jobs from a user:

.. code-block:: console

    $ qdel -u user_name

To delete a specific job:

.. code-block:: console

    $ qdel job_id

To print error message from a job:

.. code-block:: console

    $ qstat -j job_id | grep "error"

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

.. parsed-literal::

    #$ -S /bin/sh
