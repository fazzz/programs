#!/bin/sh

if [ -z $2 ]; then
    echo "USAGE: '$0' sourcename libname"
    exit
fi

source=$1
lib=$2

makelib=${PROG}/compile/makeLIB.sh
echo "#!/bin/sh" > ${HOME}/mylib/makelib_${lib}.sh
echo "chmod +x" ${makelib} >> ${HOME}/mylib/makelib_${lib}.sh
echo "cd " $(pwd) >> ${HOME}/mylib/makelib_${lib}.sh
echo "sh " ${makelib}  $@ >> ${HOME}/mylib/makelib_${lib}.sh

if [ -f ${source}.o ]; then
#if [ -e ${source}.o ]; then
    rm ${source}.o
fi
if [ -e lib${lib}.a ]; then
    rm lib${lib}.a
fi
if [ -f ~/myinclude/${lib}.h ]; then
    rm ~/myinclude/${lib}.h
fi
if [ -f ~/mylib/lib${lib}.a ]; then
    rm lib${lib}.a
fi

gcc -O4 -c ${source}.c -L${HOME}/mylib -I${HOME}/myinclude -L${HOME}/lib -I${HOME}/include -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0
ar cru lib${lib}.a ${source}.o
ranlib lib${lib}.a

mv lib${lib}.a ~/mylib/
cp ${lib}.h ~/myinclude/

cat ${HOME}/mylib/makelib_${libname}.sh


