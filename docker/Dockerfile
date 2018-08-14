FROM nextsim-reqs

ENV NEXTSIMDIR=/nextsim-src
WORKDIR /nextsim-src
COPY Makefile ./
COPY contrib ./contrib/
COPY core ./core/
COPY model ./model/
SHELL ["/bin/bash", "-c"]
RUN source /root/nextsim_libs.src && make -j7
RUN source /root/nextsim_libs.src && make All -j7

CMD [ "/bin/bash" ]
