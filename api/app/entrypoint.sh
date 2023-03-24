#!/usr/bin/env bash

PORT=5000

python manage.py collectstatic --noinput

python manage.py wait_for_db
python manage.py makemigrations
python manage.py migrate
python manage.py initadmin
python manage.py db init

cat README.md > docs/docs/index.md

cd docs || exit
[ -x ./scripts/generate_changelog.sh ] || chmod +x ./scripts/generate_changelog.sh
git config --global --add safe.directory /app/docs
./scripts/generate_changelog.sh
mkdocs build --clean
cd ..

# if [ "$ENABLE_JUPYTER" == "True" ]; then
#   jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root --config notebook/ipython_config.py &
# fi

if [ "$ENABLE_JUPYTER" == "True" ]; then
  python manage.py shell_plus --notebook &
fi

if [ "$DJANGO_DEBUG" == "True" ]; then
  python manage.py runserver 0.0.0.0:${PORT}
else
  gunicorn ldm.asgi:application \
    --log-file - \
    --workers 16 \
    --worker-class uvicorn.workers.UvicornWorker \
    --timeout 300 \
    --bind 0.0.0.0:${PORT}
fi