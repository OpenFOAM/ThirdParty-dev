# Set external libraries

#export LD_LIBRARY_PATH=...:$LD_LIBRARY_PATH
#export PKG_CONFIG_PATH=.../lib/pkgconfig:$PKG_CONFIG_PATH
#export PATH=.../bin:$PATH

# ENV var used during the analysis
export GCOV_PREFIX=$PWD/src/coverage/src
export GCOV_PREFIX_STRIP=`echo $PWD/src | grep -o "/" | wc | awk '{print $1}'`
