FROM ubuntu:20.04

ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install necessary system dependencies
RUN apt-get update && \
    apt-get install -y git cmake g++ curl python3 python3-pip r-base && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*



# Install R packages
RUN R -e "install.packages(c('data.table','magrittr','dplyr'), repos='http://cran.rstudio.com/', dependencies = TRUE)"


# Clone and setup your Python project
RUN git clone https://github.com/IEDB/TCRMatch.git && \
    cd /TCRMatch/ && \
    make && \
    cd /TCRMatch/scripts && \
    ./update.sh && \
    pip3 install pyrodigal biopython pandas

# Set the default CMD to run the tcrmatch script
WORKDIR /app
RUN mv /TCRMatch/* /app/


CMD ["bash"]
