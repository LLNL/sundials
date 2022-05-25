FROM ecpe4s/ubuntu20.04-runner-x86_64:2022-02-01 as bootstrap

ENV SPACK_ROOT=/opt/spack \
    CURRENTLY_BUILDING_DOCKER_IMAGE=1 \
    container=docker

ENV DEBIAN_FRONTEND=noninteractive   \
    LANGUAGE=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN apt-get -yqq update \
 && apt-get -yqq install --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        file \
        g++ \
        gcc \
        gfortran \
        git \
        gnupg2 \
        iproute2 \
        locales \
        lua-posix \
        make \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-wheel \
        python3-venv \
        unzip \
        csh \
 && locale-gen en_US.UTF-8 \
 && pip3 install boto3 \
 && rm -rf /var/lib/apt/lists/*

RUN mkdir $SPACK_ROOT && cd $SPACK_ROOT && \
    git clone https://github.com/spack/spack.git . && git checkout e4s-22.02  && \
    mkdir -p $SPACK_ROOT/opt/spack

RUN ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/docker-shell \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/interactive-shell \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/spack-env

RUN mkdir -p /root/.spack \
 && cp $SPACK_ROOT/share/spack/docker/modules.yaml \
        /root/.spack/modules.yaml \
 && rm -rf /root/*.* /run/nologin $SPACK_ROOT/.git

# [WORKAROUND]
# https://superuser.com/questions/1241548/
#     xubuntu-16-04-ttyname-failed-inappropriate-ioctl-for-device#1253889
RUN [ -f ~/.profile ]                                               \
 && sed -i 's/mesg n/( tty -s \\&\\& mesg n || true )/g' ~/.profile \
 || true

WORKDIR /root
SHELL ["docker-shell"]

# Creates the package cache
RUN spack spec hdf5+mpi

ENTRYPOINT ["/bin/bash", "/opt/spack/share/spack/docker/entrypoint.bash"]
CMD ["interactive-shell"]

# Build stage with Spack pre-installed and ready to be used
FROM bootstrap as builder

# Install python environment
RUN python3 -m venv /opt/python-venv/sundocs && \
    . /opt/python-venv/sundocs/bin/activate && \
    pip install wheel && \
    pip install sphinx >=4.0.0 && \
    pip install sphinx-fortran && \
    pip install sphinx_rtd_theme && \
    pip install sphinxcontrib.bibtex && \
    pip install sphinx-copybutton

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  packages:" \
&&   echo "    all:" \
&&   echo "      providers:" \
&&   echo "        blas:" \
&&   echo "        - openblas+ilp64" \
&&   echo "        mpi:" \
&&   echo "        - openmpi" \
&&   echo "      compiler:" \
&&   echo "      - gcc" \
&&   echo "  specs:" \
&&   echo "  - cmake" \
&&   echo "  - petsc~double+int64 ^openblas+ilp64" \
&&   echo "  - openmpi" \
&&   echo "  - openblas+ilp64" \
&&   echo "  - superlu-mt+int64~blas" \
&&   echo "  - trilinos+float+tpetra gotype=long_long ^openblas+ilp64" \
&&   echo "  config:" \
&&   echo "    install_tree: /opt/software" \
&&   echo "  mirrors:" \
&&   echo "    local: /var/cache/spack-mirror" \
&&   echo "    E4S: https://cache.e4s.io" \
&&   echo "  repos: []" \
&&   echo "  upstreams: {}" \
&&   echo "  concretization: together" \
&&   echo "  view: /opt/view") > /opt/spack-environment/spack.yaml

RUN --mount=type=cache,target=/var/cache/spack-mirror                   \
    cd /opt/spack-environment                                           \
 && spack env activate .                                                \
 && spack buildcache keys -it                                           \
 && cd /opt/spack-environment                                           \
 && status=0                                                            \
 && spack install --reuse --no-check-signature                          \
    || spack install --reuse --no-check-signature                       \
    || spack install --reuse --no-check-signature                       \
    || status=$?                                                        \
 && mkdir -p /var/cache/spack-mirror/build_cache                        \
 && spack buildcache update-index -d /var/cache/spack-mirror            \
 && spack buildcache list --allarch --very-long                         \
    | sed '/^$/d;/^--/d;s/@.\+//;s/\([a-z0-9]*\) \(.*\)/\2\/\1/'        \
    | sort > tmp.buildcache.txt                                         \
 && spack find --format {name}/{hash} | sort                            \
    | comm -23 - tmp.buildcache.txt                                     \
    | xargs --no-run-if-empty                                           \
      spack buildcache create --allow-root --only package --unsigned    \
                              --directory /var/cache/spack-mirror       \
                              --rebuild-index                           \
 && spack clean -a                                                      \
 && exit $status

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

# Add spack installed binaries to path
ENV PATH="/opt/view/bin:${PATH}"

# Allow mpirun as root
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
# Suppress HWLOC errors
ENV HWLOC_HIDE_ERRORS=2

ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l"]
