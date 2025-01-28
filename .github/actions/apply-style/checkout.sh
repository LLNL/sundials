#!/bin/bash

###
# Attempt to find the branch of the PR from the detached head state
##

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Attempting to find branch that matches commit..."

# Get the current commit SHA
current_commit_sha=$(git rev-parse HEAD)
# List all branches containing the current commit SHA
branches=$(git branch -r --contains $current_commit_sha)

echo "all branches:"
echo "$(git branch -a)"
echo "branches with SHA $current_commit_sha:"
echo "$branches"

# Loop over the string split by whitespace
branch=""
num_branches_found=0
for _possible_branch in $branches; do
  # Skip items that start with "pull/"
  if [[ $_possible_branch == pull/* ]]; then
    continue
  fi
  if [[ $_possible_branch == origin/* ]]; then
    _possible_branch=$(echo "$_possible_branch" | sed 's/origin\///')
  fi
  echo "Possible Branch: $_possible_branch"
  branch=$_possible_branch
  num_branches_found=$((num_branches_found+1))
done

if [ "$num_branches_found" -ne 1 ]; then
  echo "Error: Unable to find a single branch that matched git sha $current_commit_sha"
  exit 1
fi

echo "Found branch: $branch"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

git checkout $branch
git submodule update --init --recursive
