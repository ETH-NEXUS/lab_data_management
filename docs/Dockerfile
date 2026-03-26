FROM squidfunk/mkdocs-material:9.5.18

RUN apk add --no-cache bash git nodejs npm \
  && npm install -g git-release-notes

WORKDIR /docs

COPY . /docs
