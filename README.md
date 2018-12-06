# AKImhof_heme

## To get the website running, it is necessary to 
- Install biopython, prettytables, django, RabbitMQ and cellery;
- Set the environment variable for the seqdhbm package;
- Run celery using hemesitefolder$ celery -A hemesite worker -l info;
- Also create a file to protect the password at hemesite/pass.env. An example was created at ./example_pass.env.
