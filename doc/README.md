# SUNDIALS Documentation

The SUNDIALS documentation is written using reStructuredText and
[Sphinx](https://www.sphinx-doc.org/). We host the generated HTML documentation
at https://sundials.readthedocs.io.

To build the documentation with Sphinx you will need Python 3.9+. Sphinx and the
necessary extensions can be installed using the requirements file i.e., `pip
install -r requirements.txt`. Additionally, building the developer
documentation requires [Graphviz](https://graphviz.org/) for generating
flowcharts.

Once you have the dependencies installed, you can choose what you want to build.
To build the so-called "superbuild" of HTML docs that includes everything
(like our readthedocs) then do

```bash
cd superbuild
make -j4
```

If you want to build the docs separately for each of the individual SUNDIALS packages,
then you should run

```bash
make -j4 html|latexpdf # build HTML docs or PDF (using latex)
```

Finally, if you just want to build the docs for one SUNDIALS package, then you can
run

```bash
cd <package>
make -j4 html|latexpdf # build HTML docs or PDF (using latex)
```
