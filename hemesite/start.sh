conda activate djangoenv
celery -A hemesite worker --detach -l debug --concurrency=1 -Q wesa,celery
celery -A hemesite beat --detach -l info --scheduler django_celery_beat.schedulers:DatabaseScheduler
python manage.py runserver 131.220.127.81:8000
