Jupyter Notebook
****************

Installation
============

To install with ``conda``:

.. code-block:: console

    $ conda install -c conda-forge notebook

To install with ``pip``:

.. code-block:: console

    $ pip install notebook

Running in server
=================

In the remote server, activate a conda environment with Jupyter Notebook and then run the following command.

.. code-block:: console

    $ jupyter notebook --no-browser --port=8080

Next, from the local machine, run the following command (cm401 is the name of the remote server).

.. code-block:: console

    $ ssh -L 8080:localhost:8080 cm401

Finally, navigate to http://localhost:8080/.

In case you get an error like this:

.. code-block:: console

    ...
    channel 3: open failed: connect failed: Connection refused
    ...

Run the following:

.. code-block:: console

    $ ssh -L 8080:127.0.0.1:8080 cm401

In case you need to run Jupyter Notebook in an interactive session (i.e. ``qlogin``) because you are in a login node as opposed to a computing node, you have to do the port forwarding twice. To do this, follow these instructions:

| 1. Start an interactive session in the remote server. Take a note of the node name (e.g. grc141) because you will need this later.
| 2. Activate a conda environment with Jupyter Notebook.
| 3. Start a Jupyter Notebook server (e.g. `jupyter notebook --no-browser --port=8080`).
| 4. Open another Terminal window, log into the server, and then do the first port forwarding: `ssh -L 8080:localhost:8080 grc141`. This step establishes the connection between your login node and computing node.
| 5. Open yet another Terminal window and do your second port forwarding: `ssh -L 8080:localhost:8080 <login_node>`. This step establises the connection between your login node and local computer. So here's the final port forwarding: local > login node > computing node.
| 6. Finally, navigate to http://localhost:8080/.

References:

- `SSH -L connection successful, but localhost port forwarding not working “channel 3: open failed: connect failed: Connection refused <https://stackoverflow.com/questions/18705453/ssh-l-connection-successful-but-localhost-port-forwarding-not-working-channel>`__
- `Jupyter Notebook on UIowa's HPCs: An Example of Using Argon <https://zhiyzuo.github.io/Argon-Jupyter/>`__

Creating notebook programmatically
==================================

We can create a notebook via Python scripting as follows.

.. code:: python3

    import nbformat as nbf

    nb = nbf.v4.new_notebook()
    text = """\
    # My first automatic Jupyter Notebook
    This is an auto-generated notebook."""

    code = """\
    %pylab inline
    hist(normal(size=2000), bins=50);"""

    nb['cells'] = [nbf.v4.new_markdown_cell(text),
                   nbf.v4.new_code_cell(code)]
    fname = 'test.ipynb'

    with open(fname, 'w') as f:
        nbf.write(nb, f)

References:

- `How to create/modify a jupyter notebook from code (python)? <https://stackoverflow.com/questions/38193878/how-to-create-modify-a-jupyter-notebook-from-code-python>`__
