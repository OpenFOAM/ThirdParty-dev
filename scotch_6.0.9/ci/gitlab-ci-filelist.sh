#!/bin/sh

SRCDIR_TO_ANALYZE="src"

echo $PWD
rm -f filelist.txt
for dir in ${SRCDIR_TO_ANALYZE}
do
    find $dir -name '*\.[ch]' >> filelist.txt
done

