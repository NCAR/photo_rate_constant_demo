FROM fedora:32

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        gcc-c++ \
        netcdf-fortran-devel \
        cmake \
        make \
        wget \
#        python \
#        python3 \
        git \
#        ncview \
    && dnf clean all

# copy the MusicBox code
COPY . /photo-demo/

# install nc4fortran
RUN git clone https://github.com/geospace-code/nc4fortran \
    && cd nc4fortran \
    && cmake -DNetCDF_ROOT=/usr/lib64/libnetcdff.so -B build \
    && cd build \
    && make install

# install json-fortran
RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.2.0.tar.gz \
    && tar -zxvf 8.2.0.tar.gz \
    && cd json-fortran-8.2.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && make install \
    && echo '/usr/local/jsonfortran-gnu-8.2.0/lib' >> /etc/ld.so.conf.d/local.conf \
    && ldconfig

RUN cd /photo-demo/src/build \
    && export JSONFORTRAN_LIB=/usr/local/jsonfortran-gnu-8.2.0/lib \
    && export JSONFORTRAN_INC=/usr/local/jsonfortran-gnu-8.2.0/lib\
    && export NC4FORTRAN_LIB=/usr/local/lib \
    && export NC4FORTRAN_INC=/usr/local/include \
    && export NETCDF_DIR=/usr/lib64 \
    && make


