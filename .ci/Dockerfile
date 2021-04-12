FROM quay.io/pypa/manylinux2010_x86_64
MAINTAINER Dilawar Singh <dilawar.s.rajput@gmail.com>

ENV PATH=$PATH:/usr/local/bin
ENV PYDIR39=/opt/python/cp39-cp39/
RUN $PYDIR39/bin/python3 -m pip install conan

RUN ln -s $PYDIR39/bin/conan /usr/local/bin/conan
RUN which conan && conan search poco

# Read PYPI_PASSWORD  
ARG PYPI_PASSWORD
ENV PYPI_PASSWORD=$PYPI_PASSWORD

RUN git config --global user.name "Dilawar Singh"
RUN git config --global user.email "dilawar@subcom.tech"
WORKDIR /root
ADD . /root/
RUN ls -ltrh
CMD cd .ci && PYPI_PASSWORD=$PYPI_PASSWORD ./build_wheels_linux.sh 
