default:
	@echo "USAGE: make up"

env_var: export GIT_VERSION=$(git describe --always)
env_var: export GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
env_var: export GIT_LASTCOMMITDATE=$(git log -1 --format=%cI)
env_var: export GIT_COMMITHASH=$(git rev-parse HEAD)

build: env_var
	@docker-compose -f docker-compose.yml -f docker-compose.dev.yml build

up: env_var	
	@docker-compose -f docker-compose.yml -f docker-compose.dev.yml up -d