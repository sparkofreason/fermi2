FROM jupyter/scipy-notebook

USER root

RUN apt-get update

RUN apt-get install r-base libgeos-dev graphviz -y

RUN R -e 'install.packages("skellam")'

# Apt doesn't install a new enough version
RUN conda install -c conda-forge proj

RUN pip install rpy2 astropy cartopy gprof2dot

# NB_USER is set by the notbook
USER $NB_USER
