default:
	@echo "USAGE: make up"
up: export GIT_VERSION=$(git describe --always)
up: export GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
up: export GIT_LASTCOMMITDATE=$(git log -1 --format=%cI)
up: export GIT_COMMITHASH=$(git rev-parse HEAD)
up:	
	@docker-compose -f docker-compose.yml -f docker-compose.dev.yml up -d