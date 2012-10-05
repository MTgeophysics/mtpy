# Git tips etc.

## Set up local git repo 

First I created a fork on GitHub of https://github.com/geophysics/MTpy,
which was created at https://github.com/nietky/MTpy

Now clone that locally:

    C:\Users\nietky\github> git clone https://github.com/nietky/MTpy.git
    Cloning into 'MTpy'...
    remote: Counting objects: 295, done.
    remote: Compressing objects: 100% (235/235), done.
    remote: Total 295 (delta 151), reused 196 (delta 52)
    Receiving objects: 100% (295/295), 519.98 KiB | 758 KiB/s, done.
    Resolving deltas: 100% (151/151), done.
    C:\Users\nietky\github> cd mtpy

In the following the current branch is shown in the command prompt. This is
very handy ([Powershell](http://stackoverflow.com/a/10106906) recipe; lots of
methods for [bash](http://lmgtfy.com/?q=git+branch+in+bash+prompt)).

    C:\Users\nietky\github\mtpy [master]> ls
        Directory: C:\Users\nietky\github\mtpy
    Mode                LastWriteTime     Length Name
    ----                -------------     ------ ----
    d----         5/10/2012   4:10 PM            core
    d----         5/10/2012   4:10 PM            doc
    d----         5/10/2012   4:10 PM            imaging
    d----         5/10/2012   4:10 PM            utils
    -a---         5/10/2012   4:10 PM         72 README.md
    -a---         5/10/2012   4:10 PM          0 __init__.py
    C:\Users\nietky\github\mtpy [master]> git st
    # On branch master
    nothing to commit (working directory clean)
    C:\Users\nietky\github\mtpy [master]> git branch
    * master

Only the master branch was pulled down, but there are more branches on the
Github repo I just cloned, they are called "remote branches". Show them with:

    C:\Users\nietky\github\mtpy [master]> git branch -a
    * master
      remotes/origin/HEAD -> origin/master
      remotes/origin/Peacock
      remotes/origin/gitignore
      remotes/origin/lk
      remotes/origin/master
      remotes/origin/reconfig_as_pkgs

If I want to play with one of those, I need to track the remote branch:

    C:\Users\nietky\github\mtpy\doc [master]> git checkout --track origin/lk
    Branch lk set up to track remote branch lk from origin.
    Switched to a new branch 'lk'
    C:\Users\nietky\github\mtpy\doc [lk]> git branch
    * lk
      master

Now I have a local copy, and can pull any future changes by using:

    $ git pull origin lk

## Add a new remote repository

Add the central repo:

    $ git remote add upstream https://github.com/geophysics/MTpy.git

## Pushing a local feature branch

First I should check which remote repo I want:

C:\Users\nietky\github\mtpy\doc [git_tips]> git remote -v
origin  https://github.com/nietky/MTpy.git (fetch)
origin  https://github.com/nietky/MTpy.git (push)
upstream        https://github.com/geophysics/MTpy.git (fetch)
upstream        https://github.com/geophysics/MTpy.git (push)

I want them on "origin" (my github fork), and the changes I want to push are on
the "git_tips" feature branch:

    $ git push -u origin git_tips
    C:\Users\nietky\github\mtpy\doc [git_tips]> git push -u origin git_tips
    Counting objects: 6, done.
    Delta compression using up to 4 threads.
    Compressing objects: 100% (4/4), done.
    Writing objects: 100% (4/4), 1.12 KiB, done.
    Total 4 (delta 1), reused 0 (delta 0)
    To https://github.com/nietky/MTpy.git
     * [new branch]      git_tips -> git_tips
    Branch git_tips set up to track remote branch git_tips from origin.

After making some changes, pushing is very similar:

    C:\Users\nietky\github\mtpy\doc [git_tips]> git push origin git_tips

## Testing a pull request/feature branch

Just track the branch associated with the pull request, as above, and test
away/make changes/push them, etc.  In this way, more than one person can work
on a single pull request. The discussion thread on Github updates automatically
Just make sure you fetch and merge the most recent changes before starting work:

    $ git fetch


