GitHub
******

Frequently used commands for GitHub
===================================

* To clone a repository:

    .. code-block:: console

        $ git clone repo_url

* To add changes to the staging area:

    .. code-block:: console

        $ git add -A

* To commit the changes:

    .. code-block:: console

        $ git commit -m "message"

* To show the current branch:

    .. code-block:: console

        $ git branch

* To create a new branch:

    .. code-block:: console

        $ git branch branch_name

* To switch to a different branch:

    .. code-block:: console

        $ git checkout branch_name

* To list all existing branches:

    .. code-block:: console

        $ git branch -a

* To clone a specific branch:

    .. code-block:: console

        $ git clone repo_url
        $ git branch -a
        $ git checkout branch_name

* To push to master:

    .. code-block:: console

        $ git push origin master

* To push to a specific branch:

    .. code-block:: console

        $ git push origin branch_name

* To update the local repo:

    .. code-block:: console

        $ git pull

* To merge a branch into master:

    .. code-block:: console

        $ git checkout master
        $ git merge branch_name

* To delete a local branche:

    .. code-block:: console

        $ git branch -d branch_name

* To delete a remote branch:

    .. code-block:: console

        $ git push origin --delete branch_name

* To update the email address of user:

    .. code-block:: console

        $ git config user.email example@gmail.com

* To delete a remote Git tag:

    .. code-block:: console

        $ git push --delete origin tag_name
