default:
	@echo "USAGE: make up"

env_var: 
export GIT_VERSION=$(git describe --always)
export GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
export GIT_LASTCOMMITDATE=$(git log -1 --format=%cI)
export GIT_COMMITHASH=$(git rev-parse HEAD)

build: env_var
	@docker-compose -f docker-compose.yml -f docker-compose.dev.yml build

up: env_var	
	@docker-compose -f docker-compose.yml -f docker-compose.dev.yml up -d

format:
	# @find ./api/app -type f -name '*.py' ! -path '*/migrations/*' -exec pycodestyle --first '{}' \;
	# @find ./api/app -type f -name '*.py' ! -path '*/migrations/*' -exec autopep8 --in-place --aggressive --aggressive '{}' \;
	@find ./api/app -type f -name '*.py' ! -path '*/migrations/*' -exec black '{}' \;
	@find ./api/app -type f -name '*.py' ! -path '*/migrations/*' -exec autoflake --in-place --remove-unused-variables --remove-all-unused-imports '{}' \;
	@flake8 ./api/app --exclude */migrations/*

mainton:
	@docker-compose -f docker-compose.yml -f docker-compose.prod.yml exec ws sh -c 'touch /web_root/.maintenance'

maintoff:
	@docker-compose -f docker-compose.yml -f docker-compose.prod.yml exec ws sh -c 'rm -f /web_root/.maintenance'

redeploy: env_var
	@git pull
	@docker-compose -f docker-compose.yml -f docker-compose.prod.yml up --build -d


ps:
	@docker ps --format "$(FORMAT)"