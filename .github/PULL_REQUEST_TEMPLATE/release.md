This is a list of tasks that need to be done before a SUNDIALS release.

## Before Merging

- [ ] Create release branch from develop in the Git repository

- [ ] If this is a major release, search the SUNDIALS code for 'DEPRECATION NOTICE' and 'SUNDIALS_DEPRECATED'. All deprecated functions should be removed (unless this is the first version that they are deprecated).

- [ ] Regenerate the Fortran 2003 interfaces. It is possible nothing will be updated.

- [ ] Update the "Changes in ..." sections in all user guides. The changes should be sorted so that major new features are above bug fixes.

- [ ] Update version numbers and release date information using the `updateVerson.sh` script. This will update the following files:

   * `CITATIONS.md`
   * `CHANGELOG.md`
   * `CMakeLists.txt`
   * `README.md`
   * `src/arkode/README`
   * `src/cvode/README`
   * `src/cvodes/README`
   * `src/ida/README`
   * `src/idas/README`
   * `src/kinsol/README`
   * `doc/arkode/examples/source/conf.py`
   * `doc/shared/versions.py`
   * `doc/shared/History.rst`
   * `doc/shared/sundials.bib`
   * `doc/sundials/biblio.bib`
   * `scripts/tarscript.sh`

   The following files are no longer maintaianed:

   * `html/main.html` (This is no longer maintained as of at least 2016)
   * `sundialsTB/install_STB.m` (This is no longer maintained as of 2016)

- [ ] Update version numbers of third party libraries in the Install Guide in doc directory.

## After Merging

- [ ] Tag the release on main

- [ ] Sync develop with main

- [ ] Update internal web pages for SUNDIALS:
   https://computing-staging.llnl.gov/user

- [ ] After web changes are pushed to production, verify content and functionality of https://computing.llnl.gov/projects/sundials.

- [ ] Add the new release to the SUNDIALS Spack package