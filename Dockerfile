FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

ENV TERM linux

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN apt-get update && \
	apt-get -y --no-install-recommends install sudo apt-utils git ca-certificates wget make build-essential zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev perl perl-base libpng-dev && \
	rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/ramonamela/PyAnnotation.git && \
	bash -c /PyAnnotation/dependencies/install_dependencies.sh
