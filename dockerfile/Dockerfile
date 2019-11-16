FROM ubuntu
MAINTAINER yibinwang '738339463@qq.com'

RUN apt-get update
RUN apt-get -y install build-essential
RUN apt-get -y install manpages-dev
RUN apt-get -y install zlib1g-dev
RUN apt-get -y install parallel
RUN apt-get -y install git
RUN apt-get -y install ncbi-blast+
RUN apt-get -y install bioperl
RUN apt-get -y install vim
RUN apt-get -y install libmeep-openmpi-dev
RUN apt-get -y install abyss

RUN adduser --disabled-password --gecos '' stu_wangyibin
USER stu_wangyibin
WORKDIR /home/stu_wangyibin
RUN mkdir software && cd software
RUN git clone https://github.com/hyattpd/Prodigal.git
RUN cd Prodigal && make 
RUN echo 'export PATH=/home/stu_wangyibin/software/Prodigal:$PATH' >> ~/.bashrc

WORKDIR /home/stu_wangyibin
RUN git clone https://github.com/liaochenlanruo/Bt_toxin_scanner2.git
RUN cd Bt_toxin_scanner2
RUN echo 'export PATH=/home/stu_wangyibin/Bt_toxin_scanner2/:$PATH' >> ~/.bashrc


CMD /bin/bash

