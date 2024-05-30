#!/bin/sh

# Try to remove wiki
rm -rf Providentia.wiki

# Cloning wiki to get Markdown files
echo "Cloning wiki..."
git clone https://earth.bsc.es/gitlab/ac/Providentia.wiki.git
cd Providentia.wiki

echo "Copying md files..."

# Move introduction
mv * ../source/

# Removing wiki
cd ..
rm -rf Providentia.wiki

# Removing unused files
cd source
rm _sidebar.md
rm Tips-and-tricks-for-developers.md
rm Connection-setup.md
rm FAQ.md

# Generating docs
cd ..
echo "Generating docs..."
make clean
make html