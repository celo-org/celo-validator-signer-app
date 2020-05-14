FROM ubuntu:18.04
LABEL maintainer="Marek Olszewski <999594+marekolszewski@users.noreply.github.com>"

VOLUME ["/code"]
WORKDIR /code

ENV BOLOS_ENV /opt/bolos-env
ENV BOLOS_SDK /opt/bolos-sdk

RUN apt-get update --fix-missing

# Ledger SDK 
RUN apt-get install -y --no-install-recommends libc6-dev-i386 python3 python-pil curl ca-certificates bzip2 xz-utils git make

RUN echo "Create install directories" && \
  mkdir ${BOLOS_ENV} ${BOLOS_SDK}

RUN echo "Install custom gcc" && \
  curl -L https://developer.arm.com/-/media/Files/downloads/gnu-rm/6-2017q2/gcc-arm-none-eabi-6-2017-q2-update-linux.tar.bz2 --output /tmp/gcc.tar.bz2 && \
  echo "e68e4b2fe348ecb567c27985355dff75b65319a0f6595d44a18a8c5e05887cc3  /tmp/gcc.tar.bz2" | sha256sum -c && \
  tar -xvf /tmp/gcc.tar.bz2 -C ${BOLOS_ENV} && \
  mv ${BOLOS_ENV}/gcc-arm-none-eabi-6-2017-q2-update ${BOLOS_ENV}/gcc-arm-none-eabi && \
  rm /tmp/gcc.tar.bz2

RUN echo "Install custom clang" && \
  curl -L https://releases.llvm.org/8.0.0/clang+llvm-8.0.0-x86_64-linux-gnu-ubuntu-18.04.tar.xz --output /tmp/clang.tar.xz && \
  echo "0f5c314f375ebd5c35b8c1d5e5b161d9efaeff0523bac287f8b4e5b751272f51  /tmp/clang.tar.xz" | sha256sum -c && \
  tar -xvf /tmp/clang.tar.xz -C ${BOLOS_ENV} && \
  mv ${BOLOS_ENV}/clang+llvm-8.0.0-x86_64-linux-gnu-ubuntu-18.04 ${BOLOS_ENV}/clang-arm-fropi && \
  rm /tmp/clang.tar.xz

RUN echo "Install Ledger Nano X SDK"
COPY sdk-nanox-1.2.4-1.5 ${BOLOS_SDK}

# Rust setup
RUN apt-get update --fix-missing
RUN apt-get install -y build-essential
RUN curl https://sh.rustup.rs -sSf | sh -s -- -v -y

ENV PATH="${PATH}:/root/.cargo/bin"

RUN rustup install nightly
RUN rustup override set nightly
RUN rustup target add thumbv6m-none-eabi

# Ledgerblue tooling
RUN apt-get install -y software-properties-common
RUN add-apt-repository universe
RUN apt-get install -y python3-pip
RUN apt-get install -y libudev-dev libusb-dev libusb-1.0-0-dev
RUN pip3 install ledgerblue

COPY ./bin/init /usr/local/bin/init
ENTRYPOINT ["init"]
