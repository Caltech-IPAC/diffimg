# Ensure VPN is NOT running as it interferes with the ubi9 download.
FROM registry.access.redhat.com/ubi9/ubi

# Install commands and tools.

RUN yum -y install gcc
RUN yum -y install g++
RUN yum -y install gfortran
RUN yum -y install autoconf
RUN yum -y install automake
RUN yum -y install libtool
RUN yum -y install zlib-devel
RUN yum -y install atlas
RUN yum -y install atlas-devel
RUN yum -y install postgresql
RUN yum -y install postgresql-devel
RUN yum -y install unzip
RUN yum -y install git
RUN yum -y install sqlite


# The following failed to install:
#RUN yum module install -y llvm-toolset
# The following failed to install:
#RUN yum -y install epel-release
# The following is required to pip install gdal??? It failed anyway.
#RUN yum -y install gdal
#RUN yum -y install gdal-devel





#RUN dnf -y install wget
#RUN wget http://s3.amazonaws.com/etc-data.koordinates.com/gdal-travisci/install-libkml-r864-64bit.tar.gz
#RUN tar xzf install-libkml-r864-64bit.tar.gz
#RUN cp -r install-libkml/include/* /usr/local/include
#RUN cp -r install-libkml/lib/* /usr/local/lib
#RUN ldconfig
#RUN wget http://download.osgeo.org/gdal/3.11.3/gdal-3.11.3.tar.gz
#RUN tar xzf gdal-3.11.3.tar.gz
#RUN cd gdal-3.11.3
#RUN ./configure --with-libkml
#RUN make
#RUN make install

#Dockerfile:43
#--------------------
#  41 |     RUN tar xzf gdal-3.11.3.tar.gz
#  42 |     RUN cd gdal-3.11.3
#  43 | >>> RUN ./configure --with-libkml
#  44 |     RUN make
#  45 |     RUN make install
#--------------------
#ERROR: failed to solve: process "/bin/sh -c ./configure --with-libkml" did not complete successfully: exit code: 127
#[1]  + exit 1     docker build --no-cache --file ~/git/diffimg/docker/Dockerfile --tag  . >&





# Install packages in Python 3.9 (/usr/bin/python3)
# Cannot install under Python 3.9: photutils

RUN yum install -y pip
RUN pip install --upgrade pip
RUN pip install --upgrade numpy
RUN pip install --upgrade astropy
RUN pip install --upgrade psycopg2-binary
RUN pip install --upgrade healpy
RUN pip install --upgrade boto3
RUN pip install --upgrade pyfftw
RUN pip install --upgrade sympy
RUN pip install --upgrade sip_tpv
RUN pip install --upgrade matplotlib
RUN pip install --upgrade sfft
RUN pip install --upgrade scipy
RUN pip install --upgrade pyvo


# Install Python 3.12 with packages (/usr/bin/python3.12)
# Cannot install under Python 3.12: sfft,photutils(?)

RUN dnf -y install python3.12
RUN python3.12 -m ensurepip
RUN python3.12 -m pip install --upgrade pip
RUN python3.12 -m pip install --upgrade numpy
RUN python3.12 -m pip install --upgrade astropy
RUN python3.12 -m pip install --upgrade psycopg2-binary
RUN python3.12 -m pip install --upgrade healpy
RUN python3.12 -m pip install --upgrade boto3
RUN python3.12 -m pip install --upgrade pyfftw
RUN python3.12 -m pip install --upgrade sympy
RUN python3.12 -m pip install --upgrade sip_tpv
RUN python3.12 -m pip install --upgrade matplotlib
RUN python3.12 -m pip install --upgrade scipy
RUN python3.12 -m pip install --upgrade pyvo
RUN python3.12 -m pip install --upgrade Cython
# The following failed...
# On Russ's Mac: ERROR: failed to solve: process "/bin/sh -c python3.12 -m pip install --upgrade photutils" did not complete successfully: exit code: 1
#RUN python3.12 -m pip install --upgrade gdal
#RUN python3.12 -m pip install --upgrade "photutils[all]"


# Install AWS CLI

RUN cd /tmp && \
    curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    cd /tmp && \
    rm -rf /tmp/*

# Install ps command.

RUN yum -y install procps

# Checkout RAPID git repo.

RUN git clone https://github.com/Caltech-IPAC/rapid code
WORKDIR /code

# Build C executables in RAPID software stack.

ARG RAPID_SW=/code
RUN cd /code/c/builds && \
    ./build_inside_container.sh

# Configure run-time environment.

ENV PATH=/code/c/bin:/root/.local/bin:/root/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV LD_LIBRARY_PATH=/code/c/lib
ENV PYTHONPATH=/code
ENV PYTHONUNBUFFERED=1

WORKDIR /tmp

# Make the entrypoint bash so arbitrary commands can override the default.

ENTRYPOINT ["bash"]

#CMD ["/code/pipeline/awsBatchSubmitJobs_runSingleSciencePipeline.sh"]
