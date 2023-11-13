FROM nginx:stable-alpine

# Remove default configuration from Nginx
RUN rm /etc/nginx/conf.d/default.conf

ARG config=nginx.conf

# Overwrite the NGINX conf
COPY ./docker/nginx/$config /etc/nginx/conf.d/
