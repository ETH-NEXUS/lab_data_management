#!/usr/bin/env bash

PORT=5000

python manage.py collectstatic --noinput
python manage.py wait_for_db
python manage.py makemigrations
python manage.py migrate
python manage.py initadmin
python manage.py db init

if [ "$ENABLE_JUPYTER" == "True" ]; then
  jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root --config notebook/ipython_config.py &
fi

if [ "$DEV" == "True" ]; then
  python manage.py runserver 0.0.0.0:${PORT}
else
  gunicorn labmgmt.asgi:application \
    --log-file - \
    --workers 16 \
    --worker-class uvicorn.workers.UvicornWorker \
    --timeout 300 \
    --reload \
    --bind 0.0.0.0:${PORT}
fi

