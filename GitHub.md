# GitHub

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)

## Frequently used commands <a name="Frequently-used-commands"></a>

To clone a repository:

```
git clone repo_url
```

To add changes to the staging area:

```
git add -A
```

To commit the changes:

```
git commit -m "message"
```

To show the current branch:

```
git branch
```

To create a new branch:

```
git branch branch_name
```

To switch to a different branch:

```
git checkout branch_name
```

To list all existing branches:

```
git branch -a
```

To clone a specific branch:

```
git clone repo_url
git branch -a
git checkout branch_name
```






```
Push to master.
$ git push origin master

Push to a specific branch.
$ git push origin <name_of_branch>

Update the local repo.
$ git pull

Merge a branch into master.
$ git checkout master
$ git merge <name_of_branch>

Delete a local branche.
$ git branch -d <name_of_branch>

Deleting a remote branche.
$ git push origin --delete <name_of_branch>

$ git add requirements_dev.txt
$ git config user.email <sbstevenlee@gmail.com>




# Delete a remote Git tag.
$ git push --delete origin <tag_name>







$ pip install -e git+https://github.com/sbslee/stargazer#egg=stargazer
$ python -m pip install git+https://github.com/sbslee/stargazer
$ python -m pip install git+https://github.com/sbslee/stargazer.git@1307e7094251fc8b0335ef716b4fc2be7b041658
$ pip install sphinx
$ pip install sphinx_rtd_theme
$ make html
$ make clean
```
