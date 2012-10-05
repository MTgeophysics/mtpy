Example of setting up a local copy of a GitHub repository:

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

Showing remote branches:

    C:\Users\nietky\github\mtpy [master]> git branch -a
    * master
      remotes/origin/HEAD -> origin/master
      remotes/origin/Peacock
      remotes/origin/gitignore
      remotes/origin/lk
      remotes/origin/master
      remotes/origin/reconfig_as_pkgs

Track a remote branch locally:

    C:\Users\nietky\github\mtpy\doc [git]> git checkout -b lk origin/lk
    Branch lk set up to track remote branch lk from origin.
    Switched to a new branch 'lk'
