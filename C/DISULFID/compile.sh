
if [ -z $3 ]; then
    echo "USAGE: ./'$0' o/g/b  programname sourcefilename1 sourcefilename2 ......"
    exit
fi

gcc -o ${programname}g -pg -g ${sourcefilename[@]} -L${HOME}/lib -I${HOME}/include -L${HOME}/mylib -I${HOME}/myinclude -Xlinker -rpath -Xlinker ${HOME}/mylib -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0 -I${HOME}/lib/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -L/lib64  -lEF ;
mv ${programname}g ${HOME}/mybin

gcc -o ${programname} -O4 ${sourcefilename[@]} -L${HOME}/lib -I${HOME}/include -L${HOME}/mylib -I${HOME}/myinclude -Xlinker -rpath -Xlinker ${HOME}/mylib -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0 -I${HOME}/lib/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -L/lib64  -lEF ;
mv ${programname} ${HOME}/mybin
