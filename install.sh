#!/bin/sh

INSTALL_DIR=`pwd`

#if [ ! -z "${1}" ]; then
#   INSTALL_DIR=${1}
#fi

     
case "$(uname -s)" in

   Darwin)

     git clone https://github.com/simongog/sdsl-lite.git
     cd sdsl-lite
     ./install.sh ${INSTALL_DIR}
     cd ..

     wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2
     tar -xf htslib-1.3.2.tar.bz2
     cd htslib-1.3.2
     ./configure --libdir=${INSTALL_DIR}/lib CC=clang
     make 
     make prefix=${INSTALL_DIR} install
     cd ..
     ;;

   Linux)
   
     git clone https://github.com/simongog/sdsl-lite.git
     cd sdsl-lite
     ./install.sh ${INSTALL_DIR}
     cd ..

     wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2
     tar -xf htslib-1.3.2.tar.bz2
     cd htslib-1.3.2
     ./configure --libdir=${INSTALL_DIR}/lib  
     make
     make prefix=${INSTALL_DIR} install
     cd ..
     ;;

   # Add here more strings to compare
   # See correspondence table at the bottom of this answer

   *)
     echo 'unsupported OS' 
     ;;
esac
