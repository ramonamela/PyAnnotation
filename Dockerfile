FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

ENV TERM linux

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN apt-get update \
    && apt-get install -y \
        vim

RUN sed -i '66d;65d;64d;63d;62d;61d;60d' ~/.bashrc && \
	sed -i '56d;55d;54d;53d;52d' ~/.bashrc && \
	sed -i '50d;49d;48d;47d;46d;45d;44d;43d;42d;41d' ~/.bashrc && \
	sed -i '29d;28d;27d' ~/.bashrc && \
	sed -i '6d' ~/.bashrc

RUN apt-get update && \
	apt-get -y --no-install-recommends install sudo apt-utils git ca-certificates wget make build-essential zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev perl perl-base libpng-dev && \
	rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/ramonamela/PyAnnotation.git && \
	bash -c /PyAnnotation/dependencies/install_dependencies.sh && \
	bash -c /PyAnnotation/useful_commands/download_cache_files.sh
