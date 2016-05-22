#!/bin/bash

# ========================================================================
# SUMMARY
# ========================================================================
#
# "pre-commit-git-diff-docx.sh": Small git (https://git-scm.com/)
# hook. It works in combination with another hook,
# "pre-commit-git-diff-docx.sh".
#
# Together, they keep a Markdown (.md) copy of .docx files so that git
# diffs of the .md files show the changes in the document (as .docx
# files are binaries, they produce no diffs that can be checked in
# emails or in the repository's commit page).
#
# ========================================================================
# DEPENDENCIES
# ========================================================================
#
# pre-commit-git-diff-docx.sh
# pandoc (http://pandoc.org/)
#
# ========================================================================
# INSTALLATION
# ========================================================================
#
# See pre-commit-git-diff-docx.sh
#
# ========================================================================
# DETAILS:
# ========================================================================
#
# This script checks whether file .commit-amend-markdown exists. If it
# exists, it amends the previous commit adding the names of .md files
# inside it.

# Author: Ramon Casero <rcasero@gmail.com>
# Version: 0.1.0
# Copyright © 2016 University of Oxford
# 
# University of Oxford means the Chancellor, Masters and Scholars of
# the University of Oxford, having an administrative office at
# Wellington Square, Oxford OX1 2JD, UK. 
#
# This file is part of Gerardus.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. The offer of this
# program under the terms of the License is subject to the License
# being interpreted in accordance with English Law and subject to any
# action against the University of Oxford being under the jurisdiction
# of the English Courts.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# go to the top directory of this project, because filenames are
# referred to that location
cd `git rev-parse --show-toplevel`

# check whether the commit included .docx files that were converted to
# markdown (.md) format
if [ -a .commit-amend-markdown ]
then

    # add Mardown versions (.md) of the .docx files to amend the
    # commit
    cat .commit-amend-markdown | xargs git add || {
    	echo "Git cannot add Markdown files to amend the commit";
    	exit 1;
    }

    # delete the file with the list of Markdown files to avoid an
    # infinite loop
    rm .commit-amend-markdown

    # add the .md file by amending the last commit
    ## --no-verify: prevent infinite loop, don't go into the pre-commit
    ##              hook again
    echo Amend last commit adding .md files
    git commit --amend -C HEAD --no-verify || {
    	echo "Git cannot amend the commit";
    	exit 1;
    }

fi
exit
