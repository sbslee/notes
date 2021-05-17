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
