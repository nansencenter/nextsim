#! /bin/sh

# remove the links to the data
find . -type l -exec unlink {} \;
