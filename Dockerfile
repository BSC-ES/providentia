FROM ubuntu:20.04

MAINTAINER Alba Vilanova Cortezon <alba.vilanova@bsc.es>

RUN apt-get update \
&& apt-get install wget -y \
&& apt-get install ffmpeg libsm6 libxext6 -y \
&& apt-get install ghostscript -y \
&& apt-get install libx11-dev -y \
&& apt-get install libqt5gui5 -y

RUN wget https://repo.continuum.io/archive/Anaconda3-2024.02-1-Linux-x86_64.sh \
&& bash Anaconda3-2024.02-1-Linux-x86_64.sh -b \
&& rm Anaconda3-2024.02-1-Linux-x86_64.sh

ENV PATH /root/anaconda3/bin:$PATH
RUN conda update --all -y

RUN conda install -c conda-forge cartopy -y \
&& conda install -c conda-forge jupyterlab -y

RUN pip install -r /tmp/Providentia/requirements.txt

WORKDIR /tmp/
