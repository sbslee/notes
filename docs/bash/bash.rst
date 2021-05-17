Bash
****

Bash configuration
==================

The ``.bashrc`` file is used to provide a place where you can set up variables, functions and aliases, define your (PS1) prompt and define other settings that you want to use every time you open a new terminal window. The following command will activate the configuration:

.. code-block:: console

   $ source .bashrc

There is also the ``.bash_profile`` file, which is executed for login shells, while ``.bashrc`` is executed for interactive non-login shells. When an installed program cannot be called from the command line, add the line ``export PATH=~/.local/bin:$PATH`` to the ``.bash_profile`` file.
