FROM tiangolo/uwsgi-nginx-flask:python3.7

WORKDIR /app

COPY ./app /app

RUN pip install matplotlib sympy numpy flask-restful flask-cors

