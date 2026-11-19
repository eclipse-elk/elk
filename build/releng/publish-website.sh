#!/bin/bash
#
# This script clones the website repository, copies the (previously built) site
# over to the repository, and, if necessary, pushes any changes. Execute this
# script in the directory it's placed in.

echo "Website Publisher"
echo "================="

# Change to the directory the ELK repository is a subfolder of
cd ../../..

# Check if the website exists in the first place and abort otherwise
if [ ! -d elk/docs/public ]
then
  echo "docs/public/ not found"
  exit 1
fi

# Produce a shallow clone of the repository
git clone --depth 1 git@github.com:eclipse-elk/elk-website.git website

# Synchronize the website over to the repository
rsync -crv --delete --exclude=".*/" --exclude=".*" --exclude="CODE_OF_CONDUCT.md" --exclude="CONTRIBUTING.md" --exclude="LICENSE.md" --exclude="NOTICE.md" --exclude="README.md" --exclude="SECURITY.md" elk/docs/public/ website


# Change into the git repository folder
cd website

# If nothing has changed, abort
git add -A
if ! git diff --cached --exit-code; then
  echo "Changes have been detected, publishing to repository."

  git config --global user.email "elk-bot@eclipse.org"
  git config --global user.name "ELK Bot"
  git commit -m "Nightly website build job"
  git log --graph --abbrev-commit --date=relative -n 5
  git push origin master

else
  echo "No changes have been detected since last build, nothing to publish"
fi
