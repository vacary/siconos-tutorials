# Pull image from Siconos project registry (with all required deps for Siconos)
FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04
WORKDIR /home
# Get siconos options file
COPY siconos-options.cmake /home/siconos-options.cmake
# Clone Siconos
RUN git clone https://github.com/siconos/siconos.git
RUN cd siconos ; git rev-parse --short HEAD > $HOME/siconos-commit-number.txt
RUN mkdir -p /home/build
# Configure (default), build and install (system wide)
RUN cd /home/build ; cmake /home/siconos -DUSER_OPTIONS_FILE=/home/siconos-options.cmake -DCMAKE_INSTALL_PREFIX=/home/install-siconos
RUN cd /home/build ; make -j 4 ; make install
ENV SICONOS_INSTALL_DIR /home/install-siconos
RUN apt update  && apt install -y -qq libboost-timer-dev
RUN apt clean && apt autoremove
