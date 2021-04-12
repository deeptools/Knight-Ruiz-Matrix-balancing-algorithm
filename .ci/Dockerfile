FROM quay.io/pypa/manylinux2010_x86_64
MAINTAINER Dilawar Singh <dilawar.s.rajput@gmail.com>

ENV PATH=$PATH:/usr/local/bin

# Read PYPI_PASSWORD  
ARG PYPI_PASSWORD
ENV PYPI_PASSWORD=$PYPI_PASSWORD

RUN git config --global user.name "Dilawar Singh"
RUN git config --global user.email "dilawar@subcom.tech"
WORKDIR /root
ADD . /root/
RUN ls -ltrh
CMD cd .ci && PYPI_PASSWORD=$PYPI_PASSWORD ./build_wheels_linux.sh 
