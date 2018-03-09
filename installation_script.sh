#gcc5
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-5 g++-5
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 1


#cmake 3
sudo apt-get install build-essential
wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz
tar xf cmake-3.2.2.tar.gz
cd cmake-3.2.2
./configure
make
sudo make install
sudo apt-get install cmake-curses-gui
##FFTW
wget http://www.fftw.org/fftw-3.3.6-pl2.tar.gz 
tar xf fftw-3.3.6-pl2.tar.gz 
cd fftw-3.3.6-pl2
./configure
make
sudo make install


#tfhe starts
git clone https://github.com/tfhe/tfhe.git
cd tfhe
mkdir build
cd build
ccmake ../src
TFHE_PREFIX=/usr/local #the prefix where you installed tfhe
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$TFHE_PREFIX/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$TFHE_PREFIX/include
export LIBRARY_PATH=$LIBRARY_PATH:$TFHE_PREFIX/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TFHE_PREFIX/lib
##install R
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt-get update
sudo apt-get install r-base


###for NFLlib
##install GMP
wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.xz
tar -xf gmp-6.1.2.tar.xz
cd gmp-6.1.2/
sudo apt-get install m4
./configure
make
sudo make install


##install MPFR
wget http://www.mpfr.org/mpfr-current/mpfr-3.1.5.tar.xz
tar -xf mpfr-3.1.5.tar.xz
cd mpfr-3.1.5
./configure
make
sudo make install


##NFL
sudo apt-get install libgmp-dev
git clone https://github.com/quarkslab/NFLlib.git
Cd NFLlib
mkdir _build
cd _build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/nfllib
make
make test
sudo make install
