FROM ubuntu:latest

# TODO: More build libraries needed (see cyan-waterbody)?
RUN apt-get update && \
	apt-get install -y \
		wget \
		gcc \
		build-essential

# Manual installation of Mamba:
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh"
RUN /bin/bash Mambaforge-Linux-x86_64.sh -b -p /opt/mamba
RUN rm Mambaforge-Linux-x86_64.sh

ENV PATH /opt/mamba/bin:$PATH

COPY . /src

RUN chmod 755 /src/start_flask.sh

COPY uwsgi.ini /etc/uwsgi/uwsgi.ini

RUN /opt/mamba/bin/mamba init

RUN mamba env create -f /src/pkasolver/devtools/conda-envs/test_env.yaml

RUN conda run --no-capture-output -n pka_prediction python /src/pkasolver/setup.py install

RUN conda run --no-capture-output -n pka_prediction pip install -r /src/requirements.txt

CMD ["sh", "start_flask.sh"]