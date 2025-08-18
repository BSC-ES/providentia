#!/bin/sh

# Try to remove wiki, uploads and previous docs
rm -rf providentia.wiki
rm -rf source/uploads
rm -rf source/*.md

# Cloning wiki to get Markdown files
echo "Cloning wiki..."
git clone https://earth.bsc.es/gitlab/ac/providentia.wiki.git
cd providentia.wiki

echo "Copying md files..."

# Move introduction
mv * ../source/

# Removing wiki
cd ..
rm -rf providentia.wiki

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