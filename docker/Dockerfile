FROM test/test1

ENV NEXTSIMDIR=/nextsim-src
WORKDIR /nextsim-src
COPY Makefile ./
COPY contrib ./contrib/
COPY core ./core/
COPY model ./model/
SHELL ["/bin/bash", "-c"]
#RUN source /root/nextsim_libs.src && make -j7

CMD [ "/bin/bash" ]
