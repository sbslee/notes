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

To push to master:

```
git push origin master
```

To push to a specific branch:

```
git push origin branch_name
```

To update the local repo:

```
git pull
```

To merge a branch into master:

```
git checkout master
git merge branch_name
```

To delete a local branche:

```
git branch -d branch_name
```

To deleting a remote branch:

```
git push origin --delete branch_name
```

To update the email address of user:

```
git config user.email example@gmail.com
```

To delete a remote Git tag:

```
git push --delete origin tag_name
```
