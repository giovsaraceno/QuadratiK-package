# Contributing to QuadratiK package

*Version: 1.1.2*

First of all, thanks for considering contributing to `QuadratiK`! 👍

`QuadratiK` is a comprehensive statistical analysis package providing a comprehensive set of goodness-of-fit tests and a clustering technique for spherical data using kernel-based quadratic distances.

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

### Share the love ❤️

Think `QuadratiK` is useful? Let others discover it, by telling them in person, via Twitter, ResearchGate or a blog post.

Using `QuadratiK` for a paper you are writing? [Cite it](https://arxiv.org/abs/2402.02290).

### Ask a question

Using `QuadratiK` and got stuck? Browse the [documentation](https://github.com/giovsaraceno/QuadratiK-package) to see if you can find a solution. Post your questions and requests as an [issue on GitHub](https://github.com/giovsaraceno/QuadratiK-package/issues/new). We'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [gsaracen\@buffalo.edu](mailto:gsaracen@buffalo.edu).

### Propose an idea 💡

Have an idea for a new `QuadratiK` feature? Take a look at the [documentation](https://github.com/giovsaraceno/QuadratiK-package) and [issue list](https://github.com/giovsaraceno/QuadratiK-package/issues) to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub](https://github.com/giovsaraceno/QuadratiK-package/issues/new). Please

-   Explain in detail how it would work.
-   Keep the scope as narrow as possible.
-   Functions are described as comments near their code and translated to documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). Write documentation using roxygen2 package.
-   Test your changes with goodpractice before submitting.

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1.  Fork [this repo](https://github.com/giovsaraceno/QuadratiK-package) and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2.  If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3.  Open the RStudio project file (`.Rproj`).
4.  Make your changes:
    -   Write your code.
    -   Test your code (bonus points for adding unit tests).
    -   Document your code (see function documentation above).
    -   Check your code with `devtools::check()` and aim for 0 errors and warnings.
5.  Commit and push your changes.
6.  Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).

### Report a bug

Using `QuadratiK` and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub](https://github.com/giovsaraceno/QuadratiK-package/issues/new) so we can fix it. A good bug report makes it easier for us to do so, so please include:

-   Your operating system name and version (e.g. Mac OS 10.13.6).
-   Any details about your local setup that might be helpful in troubleshooting.
-   Detailed steps to reproduce the bug.

#### The website

[This website](https://github.com/giovsaraceno/QuadratiK-package) is generated with [`pkgdown`](http://pkgdown.r-lib.org/). That means we don't have to write any html: content is pulled together from documentation in the code, vignettes, [Markdown](https://guides.github.com/features/mastering-markdown/) files, the package `DESCRIPTION` and `_pkgdown.yml` settings. If you know your way around `pkgdown`, you can [propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to improve documentation. If not, [report an issue](https://github.com/giovsaraceno/QuadratiK-package/issues/new) and we can point you in the right direction.

## Workflow

The `QuadratiK` package has been engineered with a combination of R and C++ code to balance usability with high-speed computation, catering to complex data analysis needs. The package has undergone extensive testing with various examples. Continuous integration practices have been employed to maintain high standards of reliability and functionality.

`QuadratiK` was first released on CRAN on February 23, 2024. Regular updates and bug fixes are planned to continually enhance the package's functionality and user experience. We are actively planning to include additional methods based on the kernel-based quadratic-distance. One of our primary goals is to make `QuadratiK` increasingly user-friendly, with improvements to the user experience and the layout of the outputs. User feedback is highly valued and will be a key driver of future development.

This Life Cycle Statement is subject to periodic review and will be updated to reflect the evolving nature of `QuadratiK`.
