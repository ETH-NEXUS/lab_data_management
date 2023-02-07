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
	@find ./api/app -type f -name '*.py' ! -path '*/migrations/*' -exec autoflake --in-place --remove-unused-variables --remove-all-unused-imports '{}' \;
	@find ./api/app -type f -name '*.py' ! -path '*/migrations/*' -exec pycodestyle --first '{}' \;
	@find ./api/app -type f -name '*.py' ! -path '*/migrations/*' -exec autopep8 --in-place --aggressive --aggressive '{}' \;
	@flake8 ./api/app --exclude */migrations/*