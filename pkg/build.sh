#!/bin/bash
cd ..
rm -rf quadmod-release
cp -r pkg quadmod-release
##grep -v quadmodData quadmod/DESCRIPTION | grep -v Remotes > quadmod-release/DESCRIPTION
##rm quadmod-release/tests/testthat/test-quadmodData.R
PKG_TGZ=$(R CMD build quadmod-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
