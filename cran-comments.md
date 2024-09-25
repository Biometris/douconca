## Patch release

Patch release addressing check error on several CRAN machines. After previous submission no longer check full output objects to avoid issues with different signs in eigenvalue decompositions for different configurations.

## Test environments

* local Windows 11 install, R 4.4.1
* winbuilder (develop)
* macbuilder (release)
* Ubuntu (via github actions, devel and release)
* macOS (via github actions, release)
* mkl (via github actions)

## R CMD check results

0 errors | 0 warnings | 0 notes
