FROM rocker/binder:latest
MAINTAINER ross@ecohealthalliance.org

## Copies your repo files into the Docker Container
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
RUN install2.r deSolve

## Become normal user again
USER ${NB_USER}



## Run an install.R script, if it exists.

RUN if [ -f DESCRIPTION ]; then R --quiet -e "devtools::install(dep=TRUE)"; fi
