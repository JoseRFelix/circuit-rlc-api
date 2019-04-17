FROM tiangolo/uwsgi-nginx-flask:python3.7

WORKDIR /app

COPY ./app /app

EXPOSE 443

RUN pip install matplotlib sympy numpy flask-restful flask-cors seaborn pyopenssl

