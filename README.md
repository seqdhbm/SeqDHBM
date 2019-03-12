# AKImhof_heme

## To get the website running, it is necessary to 
- Install biopython, prettytables and django
- Install RabbitMQ, cellery and django-celery-beat:
  - sudo apt-get install rabbitmq-server
  - sudo rabbitmqctl add_user myuser mypassword
  - sudo rabbitmqctl add_vhost myvhost
  - sudo rabbitmqctl set_permissions -p myvhost myuser ".*" ".*" ".*"
  - to start: sudo rabbitmq-server -detached
  - to stop: sudo rabbitmqctl stop
  - pip install celery 
  - pip install django_celery_beat
  - pip install Flask-SQLAlchemy
- If there is a problem, install django celery results...
  - pip install django-celery-results
  - python manage.py migrate django_celery_results
- ... and update the database model
  - python manage.py makemigrations
  - python manage.py migrate
- Set the environment variable for the seqdhbm package;
- Run celery using hemesitefolder$ celery -A hemesite worker -l info;
- Also create a file to protect the password at hemesite/pass.env. An example was created at ./example_pass.env.
