Quick guide how to develop the documentation, which can be automatically published online

To build in local machine:
cd docs
make html

after the successful sphinx-build, view html files in \docs\build\html\

To modify/add auto docs, edit the files:
docs/source/conf.py
docs/source/index.rst
docs/source/core.rst, analysis.rst, model.rst, imaging.rst, etc

make html
And view it in your local dev environment

Once you are comfortable with the new changes in docs, commit and push into github repo (thedevelop branch as default)
Then the system will automatically build html docs and publish it in https://mtpy2.readthedocs.io/en/develop/
(Just wait a few minutes to see the new changes take effect)


For all these to happen automatically behind the scene,
a developer (FeiZhang) has to set up hooksin both Github and
readthedocs service at: https://readthedocs.org/


See also
http://www.sphinx-doc.org/en/master/contents.html
Web Link: https://readthedocs.org/projects/mtpy2/builds/7813943/
click view raw
https://readthedocs.org/api/v2/build/7813943.txt
