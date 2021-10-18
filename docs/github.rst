GitHub
******

Frequently used commands for GitHub
===================================

Clone a repository:

.. code-block:: console

    $ git clone repo_url

Add all changes to the staging area:

.. code-block:: console

    $ git add -A

Commit changes:

.. code-block:: console

    $ git commit -m "put message here"

Show current branch:

.. code-block:: console

    $ git branch

Create a new branch:

.. code-block:: console

    $ git branch branch_name

Switch to an existing branch:

.. code-block:: console

    $ git checkout branch_name

Create and switch to a new branch simultaneously:

.. code-block:: console

    $ git checkout -b branch_name

List all existing branches:

.. code-block:: console

    $ git branch -a

Clone a specific branch:

.. code-block:: console

    $ git clone repo_url
    $ git branch -a
    $ git checkout branch_name

Push commits to master:

.. code-block:: console

    $ git push origin master

Push commits to a specific branch:

.. code-block:: console

    $ git push origin branch_name

Update the local repo:

.. code-block:: console

    $ git pull

Merge a branch into master:

.. code-block:: console

    $ git checkout master
    $ git merge branch_name

Delete a local branche:

.. code-block:: console

    $ git branch -d branch_name

Delete a remote branch:

.. code-block:: console

    $ git push origin --delete branch_name

Update email address of user:

.. code-block:: console

    $ git config user.email example@gmail.com

Delete a remote Git tag:

.. code-block:: console

    $ git push --delete origin tag_name

Problem with the SSL CA cert
============================

I encountered below error in the company's server:

.. code-block:: text

    $ git clone https://github.com/sbslee/dokdo
    Cloning into 'dokdo'...
    fatal: unable to access 'https://github.com/sbslee/dokdo': Problem with the SSL CA cert (path? access rights?)

I solved it with:

.. code-block:: text

    $ git config --global http.sslVerify false
    $ git clone https://github.com/sbslee/dokdo
    
I then restored the verification setting:
    
.. code-block:: text

    $ git config --global http.sslVerify true
