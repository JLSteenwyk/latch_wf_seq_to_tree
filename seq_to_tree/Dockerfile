FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:02ab-main

# Its easy to build binaries from source that you can later reference as
# subprocesses within your workflow.
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

RUN curl -L https://mafft.cbrc.jp/alignment/software/mafft-7.490-with-extensions-src.tgz -o mafft-7.490-src.tgz &&\
    gunzip -cd mafft-7.490-src.tgz | tar xfv - &&\
    cd mafft-7.490-with-extensions/core/ &&\
    make clean &&\
    make &&\
    su &&\
    make install &&\
    rm -rf mafft-7.490-with-extensions mafft-7.490-src.tgz 

RUN curl -L https://github.com/iqtree/iqtree2/releases/download/v2.2.0/iqtree-2.2.0-Linux.tar.gz -o iqtree-2.2.0-Linux.tar.gz &&\
    gunzip -cd iqtree-2.2.0-Linux.tar.gz | tar xfv - &&\
    mv iqtree-2.2.0-Linux/bin/iqtree2 /usr/bin &&\
    rm -rf iqtree-2.2.0-Linux iqtree-2.2.0-Linux.tar.gz

#create local directory to store output of cmd later
RUN mkdir /root/seq_to_tree/

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
WORKDIR /root
