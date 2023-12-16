FROM scilus/scilus:1.6.0

LABEL version="ChildBrainFlow-1.0.0"

RUN wget -O FS_BN_GL_SF_utils.tar.gz "https://www.dropbox.com/scl/fi/6s1tc4eanf2sutejw7fkd/FS_BN_GL_SF_utils.tar.gz?rlkey=3gvhvpepv7ldkqef3go10cb5e&dl=0" && \
    tar -xzvf FS_BN_GL_SF_utils.tar.gz && \
    rm FS_BN_GL_SF_utils.tar.gz

# Installing dependencies.
RUN apt-get update && \
    apt-get install -y csh tcsh && \
    apt-get install -y libglu1-mesa libxt6 libxmu6 libgl1 freeglut3-dev && \
    echo "/usr/local/lib" >> /etc/ld.so.conf && \
    ldconfig
RUN wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2 \
    tar xjf parallel-latest.tar.bz2 \
    cd parallel-* && ./configure && make && make install \
    echo 'will cite' | parallel --citation 1> /dev/null 2> /dev/null &

# Setup Conda for FastSurfer
ARG CONDA_FILE=Miniconda3-py38_4.11.0-Linux-x86_64.sh
RUN wget --no-check-certificate -qO ~/miniconda.sh https://repo.continuum.io/miniconda/${CONDA_FILE} && \
     chmod +x ~/miniconda.sh && \
     ~/miniconda.sh -b -p /opt/conda && \
     rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH 

# Setup FastSurfer
WORKDIR /
RUN wget -O FastSurfer.tar.gz "https://github.com/Deep-MI/FastSurfer/archive/refs/tags/v2.1.2.tar.gz" && \
    tar -xzvf FastSurfer.tar.gz && \
    mv FastSurfer-2.1.2 FastSurfer && \
    rm FastSurfer.tar.gz && \
    conda env create -f FastSurfer/fastsurfer_env_cpu.yml

# Install conda-pack:
RUN conda install -c conda-forge conda-pack

# Use conda-pack to create a standalone env in /venv:
RUN conda-pack -n fastsurfer_cpu -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

# Fix paths.
RUN /venv/bin/conda-unpack
ENV PATH /venv/bin:$PATH

# Installing freesurfer on top of scilus:1.6.0
WORKDIR /
RUN wget "https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.3.2/freesurfer-linux-centos7_x86_64-7.3.2.tar.gz" -O fs.tar.gz && \
    tar --no-same-owner -xzvf fs.tar.gz && \
    rm fs.tar.gz
RUN wget -O freesurfer/license.txt "https://www.dropbox.com/scl/fi/0s8lp6lydyd0rxawxb4jm/license.txt?rlkey=hz54oc0d4sor69avqphtrjvgn&dl=0"

# Setup freesurfer env
ENV OS=Linux
ENV PATH=${PATH}:/FastSurfer:/freesurfer/bin:/freesurfer/fsfast/bin:/freesurfer/tktools:/freesurfer/mni/bin:/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV FREESURFER_HOME=/freesurfer
ENV FREESURFER=/freesurfer
ENV SUBJECTS_DIR=/freesurfer/subjects
ENV LOCAL_DIR=/freesurfer/local
ENV FSFAST_HOME=/freesurfer/fsfast
ENV FMRI_ANALYSIS_DIR=/freesurfer/fsfast
ENV FUNCTIONALS_DIR=/freesurfer/sessions
ENV FS_LICENSE=/freesurfer/license.txt

# set default fs options
ENV FS_OVERRIDE=0
ENV FIX_VERTEX_AREA=""
ENV FSF_OUTPUT_FORMAT=nii.gz

# mni env requirements
ENV MINC_BIN_DIR=/freesurfer/mni/bin
ENV MINC_LIB_DIR=/freesurfer/mni/lib
ENV MNI_DIR=/freesurfer/mni
ENV MNI_DATAPATH=/freesurfer/mni/data
ENV MNI_PERL5LIB=/freesurfer/mni/share/perl5
ENV PERL5LIB=/freesurfer/mni/share/perl5