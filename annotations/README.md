

Install GSL (https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/):

wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz
tar -zxvf gsl-2.5.tar.gz
cd gsl-2.5
mkdir /mnt/seagate10/gsl
./configure --prefix=/mnt/seagate10/gsl
make
make check
make install

Update .bashrc file adding path to GSL lib:
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/seagate10/gsl/lib
export LD_LIBRARY_PATH

