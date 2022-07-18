FROM sagemath/sagemath:9.5

RUN sudo apt-get update \
    && sudo apt-get install -y make \
    && sudo rm -r /var/lib/apt/lists/*

WORKDIR "/home/sage/mq_estimator"
CMD ["/bin/bash"]
