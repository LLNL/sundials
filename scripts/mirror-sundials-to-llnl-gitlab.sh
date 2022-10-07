#!/usr/bin/bash

git clone --bare ssh://git@mybitbucket.llnl.gov:7999/sundials/sunrepo.git $HOME/.sundials-mirroring/sunrepo
cd $HOME/.sundials-mirroring/sunrepo
git push --mirror ssh://git@czgitlab.llnl.gov:7999/sundials/sundials.git
cd ..
rm -rf $HOME/.sundials-mirroring/sunrepo
