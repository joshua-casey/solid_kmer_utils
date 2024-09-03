mkdir -p libs
cd libs
wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
tar -xf htslib-1.20.tar.bz2
cd htslib-1.20
./configure --prefix=$(pwd)/install
make
make install
