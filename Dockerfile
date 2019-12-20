FROM perl:5.22

WORKDIR /app

ADD . /app

RUN cpanm 

EXPOSE 5000

CMD ["perl", "IgServer.psgi"]
