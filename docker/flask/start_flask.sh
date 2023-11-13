#!/bin/bash
cd /src
exec conda run -n pka_prediction uwsgi --ini /etc/uwsgi/uwsgi.ini