FROM condaforge/mambaforge:4.10.3-1

COPY ["environment.yml", "./"]

RUN apt update && \
    apt install --assume-yes wget nano

RUN mamba env update -f environment.yml

CMD echo "This is the ngs_pipelines image with versions:" $(mamba --version | tr "\n" "\ ")
            
